#!/usr/bin/env bash
# ==========================================================
# Script: run_emu_amplicons.sh
# Purpose: Install/check Emu + DBs, split FASTQ by length bins, run Emu per bin
# Notes:
#   - Requires conda/mamba (Miniforge/Conda) + seqkit
#   - Input FASTQ files live in data/
#   - Outputs to results/, logs to logs/
# ==========================================================

set -euo pipefail

# ---------------- User-tunable defaults -------------------
THREADS="${THREADS:-8}"

# Databases:
# 1) Emu official 16S (bacteria+archaea) — will be downloaded automatically if missing
EMU_DB16S_DIR="${EMU_DB16S_DIR:-$PWD/refdb/emu_16S}"   # auto-fetched

# 2) Custom Emu-format DB paths for ITS / LSU (must be provided by user if we want Emu on these)
# Each directory must contain: species_taxid.fasta  and  taxonomy.tsv
EMU_DB_ITS_DIR="${EMU_DB_ITS_DIR:-}"   # e.g., /path/to/emu_db_unite_its
EMU_DB_LSU_DIR="${EMU_DB_LSU_DIR:-}"   # e.g., /path/to/emu_db_silva_lsu

# Length bins (non-overlapping; edit if needed)
ITS_MIN="${ITS_MIN:-200}"; ITS_MAX="${ITS_MAX:-650}"
LSU_MIN="${LSU_MIN:-651}"; LSU_MAX="${LSU_MAX:-1199}"
S16_MIN="${S16_MIN:-1200}"; S16_MAX="${S16_MAX:-1900}"
AMB_MIN="${AMB_MIN:-1901}"; AMB_MAX="${AMB_MAX:-3300}"

# Env name
ENV_NAME="${ENV_NAME:-emu-env}"

# Type of sequencing
EMU_TYPE="${EMU_TYPE:-map-ont}"

# ----------------------------------------------------------
# Paths & global state
# ----------------------------------------------------------
LOGDIR="logs"
OUTDIR="results"
METADIR="metadata"
mkdir -p "$OUTDIR" "$METADIR" refdb "$LOGDIR"

RUN_LOG="${LOGDIR}/run_emu_amplicons.$(date +%Y%m%d_%H%M%S).log"
touch "$RUN_LOG"

#Root for bin outputs inside filtered/
BINS_ROOT="${BINS_ROOT:-results/filtered/bin_splits}"
mkdir -p "$BINS_ROOT"

# Marker enable flags (computed after DB checks)
ENABLE_16S=1
ENABLE_ITS=0
ENABLE_LSU=0

# ----------------------------------------------------------
# Helper: time_function
# ----------------------------------------------------------
time_function() {
    local fn="$1"
    local start=$(date +%s)
    echo ">>> Running $fn ..." | tee -a "$RUN_LOG"
    local status=0

    set +e
    $fn
    status=$?
    set -e

    local end=$(date +%s)
    local dur=$((end - start))
    mkdir -p "$LOGDIR"
    echo -e "${fn}\t${dur}" >> "${LOGDIR}/.timing.tsv"
    echo ">>> $fn completed in ${dur}s" | tee -a "$RUN_LOG"
    return $status
}

# Log helper
log()  { echo -e ">>> $*" | tee -a "$RUN_LOG" >&2; }
warn() { echo -e "!!! $*" | tee -a "$RUN_LOG" >&2; }
die()  { echo -e "xxx $*" | tee -a "$RUN_LOG" >&2; exit 1; }

# ----------------------------------------------------------
# Function: ensure_channels
#   Reset channels order and strict priority; clear index cache
# ----------------------------------------------------------
ensure_channels() {
    echo ">>> Ensuring Conda channels (conda-forge, bioconda, defaults) with strict priority..." | tee -a "$RUN_LOG"
    conda config --remove-key channels >/dev/null 2>&1 || true
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --add channels defaults
    conda config --set channel_priority strict
    echo ">>> Clearing index cache..." | tee -a "$RUN_LOG"
    conda clean --index-cache -y >/dev/null 2>&1 || true

    echo ">>> Current channels:" | tee -a "$RUN_LOG"
    conda config --show channels | tee -a "$RUN_LOG"
}

# ----------------------------------------------------------
# Create/activate env
#   - If exists: activate
#   - Else: ensure_channels + create (retry once on solver failure)
# ----------------------------------------------------------
create_env_emu() {
  log "Checking conda/mamba..."
  if ! command -v conda >/dev/null 2>&1; then
    die "conda not found. Install Miniforge/Conda first."
  fi
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"

  if conda env list | grep -qE "^${ENV_NAME}\s"; then
    log "Environment '${ENV_NAME}' already exists. Activating it..."
    conda activate "${ENV_NAME}"
    log "Environment '${ENV_NAME}' is now active."
  else
    # Only enforce channels & clean index when creating
    ensure_channels

    log "Environment '${ENV_NAME}' not found. Creating it now (strict priority)..."
    set +e
    if command -v mamba >/dev/null 2>&1; then
      mamba create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        python=3.11 emu minimap2 "seqkit>=2.6" \
        "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" -y
    else
      conda create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        python=3.11 emu minimap2 "seqkit>=2.6" \
        "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" -y
    fi
    status=$?
    set -e

    if [ $status -ne 0 ]; then
      echo "!!! Strict-priority solve failed. Retrying once with flexible channel priority..." | tee -a "$RUN_LOG"
      if command -v mamba >/dev/null 2>&1; then
        mamba create -n "${ENV_NAME}" \
          -c conda-forge -c bioconda \
          python=3.11 emu minimap2 "seqkit>=2.6" \
          "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" -y
      else
        conda create -n "${ENV_NAME}" \
          -c conda-forge -c bioconda \
          python=3.11 emu minimap2 "seqkit>=2.6" \
          "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" -y
      fi
      echo ">>> Created '${ENV_NAME}' with flexible priority." | tee -a "$RUN_LOG"
    else
      echo ">>> Environment '${ENV_NAME}' created successfully (strict priority)." | tee -a "$RUN_LOG"
    fi

    # Activate the newly created env
    conda activate "${ENV_NAME}"
    echo ">>> Environment '${ENV_NAME}' is now active." | tee -a "$RUN_LOG"
  fi

  # show which Rscript is in-path for debugging
  log "Rscript path in env: $(which Rscript)"

  # Verify required tools
  command -v emu >/dev/null 2>&1 || die "emu not available in env"
  command -v minimap2 >/dev/null 2>&1 || die "minimap2 not available in env"
  command -v seqkit >/dev/null 2>&1 || die "seqkit not available in env"

  log "Env ready: $(which python); emu=$(which emu)"
}

# ----------------------------------------------------------
# Export the active emu env to envs/emu-env.yml
# ----------------------------------------------------------
export_emu_env() {
  if [ ! -d envs ]; then
      log "Directory 'envs/' not found. Creating it..."
      mkdir -p envs
  fi
  log "Exporting environment to envs/${ENV_NAME}.yml ..."
  conda env export --name "${ENV_NAME}" > "envs/${ENV_NAME}.yml"
  log "Environment exported: envs/${ENV_NAME}.yml"
}

# ----------------------------------------------------------
# Emu 16S DB: auto-download (bacteria + archaea)
# ----------------------------------------------------------
ensure_emu16s_db() {
  if [ -s "${EMU_DB16S_DIR}/species_taxid.fasta" ] && [ -s "${EMU_DB16S_DIR}/taxonomy.tsv" ]; then
    log "Emu 16S DB already present at: ${EMU_DB16S_DIR}"
    return 0
  fi
  log "Fetching Emu default 16S DB to ${EMU_DB16S_DIR} ..."
  mkdir -p "${EMU_DB16S_DIR}"
  pip install --quiet osfclient
  export EMU_DATABASE_DIR="${EMU_DB16S_DIR}"
  pushd "${EMU_DB16S_DIR}" >/dev/null
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
    tar -xvf emu.tar
  popd >/dev/null
  [ -s "${EMU_DB16S_DIR}/species_taxid.fasta" ] || die "Failed to populate Emu 16S DB"
  [ -s "${EMU_DB16S_DIR}/taxonomy.tsv" ] || die "Failed to populate Emu 16S DB taxonomy"
  log "Emu 16S DB ready."
}

# ----------------------------------------------------------
# ITS/LSU DBs (custom): just check presence if user provided
# ----------------------------------------------------------
check_optional_db() {
  local d="$1" ; local label="$2"
  if [ -n "${d}" ] && [ -s "${d}/species_taxid.fasta" ] && [ -s "${d}/taxonomy.tsv" ]; then
    log "${label} DB OK at: ${d}"
    return 0
  fi
  warn "No usable ${label} DB at '${d:-<empty>}' — ${label} classification will be skipped."
  return 1
}

# Compute which markers are enabled based on DB availability
set_marker_flags() {
  ENABLE_16S=1
  ENABLE_ITS=0
  ENABLE_LSU=0
  if [ -n "${EMU_DB_ITS_DIR}" ] && [ -s "${EMU_DB_ITS_DIR}/species_taxid.fasta" ] && [ -s "${EMU_DB_ITS_DIR}/taxonomy.tsv" ]; then
    ENABLE_ITS=1
  fi
  if [ -n "${EMU_DB_LSU_DIR}" ] && [ -s "${EMU_DB_LSU_DIR}/species_taxid.fasta" ] && [ -s "${EMU_DB_LSU_DIR}/taxonomy.tsv" ]; then
    ENABLE_LSU=1
  fi
  log "Markers enabled -> 16S:${ENABLE_16S} ITS:${ENABLE_ITS} LSU:${ENABLE_LSU}"
}

# ----------------------------------------------------------
# FASTQ discovery
# ----------------------------------------------------------
FASTQ_DIR_DEFAULT="${FASTQ_DIR_DEFAULT:-results/filtered}"
FASTQ_MANIFEST="${FASTQ_MANIFEST:-${METADIR}/fastq_meta.tsv}"

discover_fastqs() {
  log "Scanning ${FASTQ_DIR_DEFAULT}/ for filtered FASTQ files..."
  shopt -s nullglob
  FASTQS=( "${FASTQ_DIR_DEFAULT}"/*.fastq.gz "${FASTQ_DIR_DEFAULT}"/*.fq.gz \
           "${FASTQ_DIR_DEFAULT}"/*.fastq     "${FASTQ_DIR_DEFAULT}"/*.fq )
  shopt -u nullglob
  [ ${#FASTQS[@]} -gt 0 ] || die "No FASTQ files found in ${FASTQ_DIR_DEFAULT}/"
  log "Found ${#FASTQS[@]} FASTQ files."
}

# ----------------------------------------------------------
# Build or reuse manifest with parsed metadata
# ----------------------------------------------------------
build_fastq_meta() {
  if [ -s "$FASTQ_MANIFEST" ]; then
    log "Manifest already exists; reusing: $FASTQ_MANIFEST"
    return 0
  fi

  log "Building manifest with sample, replicate, and molecular_feature columns"
  mkdir -p "$(dirname "$FASTQ_MANIFEST")"

  {
    echo -e "file\tsample\treplicate\tmolecular_feature"
    for f in "${FASTQS[@]}"; do
      b=$(basename "$f")
      b=${b%.fastq.gz}; b=${b%.fq.gz}; b=${b%.fastq}; b=${b%.fq}
      IFS='_' read -r sample replicate molecular_feature _ <<< "$b"
      sample=${sample:-NA}
      replicate=${replicate:-NA}
      molecular_feature=${molecular_feature:-NA}
      echo -e "${f}\t${sample}\t${replicate}\t${molecular_feature}"
    done
  } > "$FASTQ_MANIFEST"

  log "Wrote manifest to $FASTQ_MANIFEST with $((${#FASTQS[@]})) entries."
}

# ----------------------------------------------------------
# Split by length bins using seqkit (gz-aware)
# - Auto-merge group labels if nominal ranges overlap
# ----------------------------------------------------------
emit_group_label() {
  local min=$1 max=$2 wanted="$3" # ITS / LSU / 16S
  echo "${wanted}_${min}-${max}"
}

maybe_merge_labels() {
  local labA="$1" minA="$2" maxA="$3"
  local labB="$4" minB="$5" maxB="$6"
  if (( minB < maxA )) && (( minA < maxB )); then
    echo "${labA}_${labB}_merged_${minA}-${maxB}"
  else
    echo ""
  fi
}

split_fastq_bins() {
  local f="$1"; local base="$2"
  local outprefix="${BINS_ROOT}/${base}"
  mkdir -p "${outprefix}/bins"

  # Build labels (no overlaps by our defaults)
  LABEL_ITS=$(emit_group_label "$ITS_MIN" "$ITS_MAX" "ITS")
  LABEL_LSU=$(emit_group_label "$LSU_MIN" "$LSU_MAX" "LSU")
  LABEL_16S=$(emit_group_label "$S16_MIN" "$S16_MAX" "16S")
  LABEL_AMB=$(emit_group_label "$AMB_MIN" "$AMB_MAX" "Ambiguous")

  MERGE_ITS_LSU=$(maybe_merge_labels "ITS" "$ITS_MIN" "$ITS_MAX" "LSU" "$LSU_MIN" "$LSU_MAX")
  if [ -n "$MERGE_ITS_LSU" ]; then LABEL_ITS="$MERGE_ITS_LSU"; LABEL_LSU="$MERGE_ITS_LSU"; fi

  log "Splitting ${f} into bins:"
  log "  - ${LABEL_ITS}"
  log "  - ${LABEL_LSU}"
  log "  - ${LABEL_16S}"
  log "  - ${LABEL_AMB}"

  # Outputs go to ${BINS_ROOT}/${base}/bins/
  local out_its="${outprefix}/bins/${base}.${LABEL_ITS}.fastq.gz"
  local out_lsu="${outprefix}/bins/${base}.${LABEL_LSU}.fastq.gz"
  local out_16s="${outprefix}/bins/${base}.${LABEL_16S}.fastq.gz"
  local out_amb="${outprefix}/bins/${base}.${LABEL_AMB}.fastq.gz"

  seqkit seq -g -m "$ITS_MIN" -M "$ITS_MAX" -o "$out_its" "$f" 2>>"$RUN_LOG" || true
  seqkit seq -g -m "$LSU_MIN" -M "$LSU_MAX" -o "$out_lsu" "$f" 2>>"$RUN_LOG" || true
  seqkit seq -g -m "$S16_MIN" -M "$S16_MAX" -o "$out_16s" "$f" 2>>"$RUN_LOG" || true
  seqkit seq -g -m "$AMB_MIN" -M "$AMB_MAX" -o "$out_amb" "$f" 2>>"$RUN_LOG" || true

  # Return TSV saved alongside the bins under ${BINS_ROOT}/${base}/
  awk -v l1="$LABEL_ITS" -v f1="$out_its" \
      -v l2="$LABEL_LSU" -v f2="$out_lsu" \
      -v l3="$LABEL_16S" -v f3="$out_16s" \
      -v l4="$LABEL_AMB" -v f4="$out_amb" '
      function nonempty(ff){ cmd="bash -lc \"[ -s " ff " ] && echo 1 || echo 0\""; cmd | getline r; close(cmd); return r==1 }
      BEGIN{
        if (nonempty(f1)) print l1 "\t" f1;
        if (nonempty(f2)) print l2 "\t" f2;
        if (nonempty(f3)) print l3 "\t" f3;
        if (nonempty(f4)) print l4 "\t" f4;
      }'
}

# ----------------------------------------------------------
# Run Emu for a given bin (chooses DB based on label)
# ----------------------------------------------------------
run_emu_for_bin() {
  local label="$1"; local binfastq="$2"; local base="$3"
  local outdir="${OUTDIR}/emu_runs/${base}/emu_${label}"
  mkdir -p "$(dirname "$outdir")"

  case "$label" in
    16S* )
      mkdir -p "$outdir"
      if ! emu abundance \
            --threads "$THREADS" \
            --db "$EMU_DB16S_DIR" \
            --output-dir "$outdir" \
            --keep-counts \
            --keep-read-assignments \
            --output-unclassified \
            --type "$EMU_TYPE" \
            "$binfastq" 2>>"$RUN_LOG" | tee -a "$RUN_LOG"; then
        warn "Emu 16S failed for ${base} (${label}); removing empty outdir."
        rmdir "$outdir" >/dev/null 2>&1 || true
      fi
      ;;
    ITS* )
      if [ "$ENABLE_ITS" -eq 1 ]; then
        mkdir -p "$outdir"
        if ! emu abundance \
              --threads "$THREADS" \
              --db "$EMU_DB_ITS_DIR" \
              --output-dir "$outdir" \
              --keep-counts \
              --keep-read-assignments \
              --output-unclassified \
              --type "$EMU_TYPE" \
              "$binfastq" 2>>"$RUN_LOG" | tee -a "$RUN_LOG"; then
          warn "Emu ITS failed for ${base} (${label}); removing empty outdir."
          rmdir "$outdir" >/dev/null 2>&1 || true
        fi
      else
        warn "Skipping ITS bin for ${base} (ITS DB not available)."
      fi
      ;;
    LSU* )
      if [ "$ENABLE_LSU" -eq 1 ]; then
        mkdir -p "$outdir"
        if ! emu abundance \
              --threads "$THREADS" \
              --db "$EMU_DB_LSU_DIR" \
              --output-dir "$outdir" \
              --keep-counts \
              --keep-read-assignments \
              --output-unclassified \
              --type "$EMU_TYPE" \
              "$binfastq" 2>>"$RUN_LOG" | tee -a "$RUN_LOG"; then
          warn "Emu LSU failed for ${base} (${label}); removing empty outdir."
          rmdir "$outdir" >/dev/null 2>&1 || true
        fi
      else
        warn "Skipping LSU bin for ${base} (LSU DB not available)."
      fi
      ;;
    Ambiguous* )
      warn "Skipping Ambiguous_long bin (${binfastq}) from classification (set RUN_AMBIGUOUS=1 to force)."
      # no mkdir here, so no empty emu_Ambiguous_* dirs appear
      ;;
    * )
      warn "Unknown label '${label}' for ${binfastq} — skipping."
      ;;
  esac
}

# ----------------------------------------------------------
# Split all FASTQs into bins (per-sample)
# ----------------------------------------------------------
split_bins_all() {
  for f in "${FASTQS[@]}"; do
    local base
    base=$(basename "$f")
    base=${base%.fastq.gz}; base=${base%.fq.gz}; base=${base%.fastq}; base=${base%.fq}
    mkdir -p "${BINS_ROOT}/${base}"
    local BIN_TSV="${BINS_ROOT}/${base}/bins_index.tsv"
    split_fastq_bins "$f" "$base" > "$BIN_TSV"
  done
}

# ----------------------------------------------------------
# Run Emu across all non-empty bins (per-sample)
# ----------------------------------------------------------
run_emu_all() {
  shopt -s nullglob
  local BIN_INDEXES=( "${BINS_ROOT}"/*/bins_index.tsv )
  shopt -u nullglob

  if [ ${#BIN_INDEXES[@]} -eq 0 ]; then
    warn "[run_emu_all] No bins_index.tsv found under ${BINS_ROOT}/*/ — nothing to run."
    return 0
  fi

  for BIN_TSV in "${BIN_INDEXES[@]}"; do
    local base
    base="$(basename "$(dirname "$BIN_TSV")")"
    [ -s "$BIN_TSV" ] || { warn "[run_emu_all] Empty ${BIN_TSV}; skipping ${base}."; continue; }

    while IFS=$'\t' read -r label binfq; do
      [ -n "$label" ] || continue
      # trim & require known labels
      label="$(echo "$label" | sed -E 's/^[[:space:]]+|[[:space:]]+$//g')"
      [[ "$label" =~ ^(16S|ITS|LSU) ]] || continue
      # skip Ambiguous here, too (defensive)
      [[ "$label" =~ ^Ambiguous ]] && continue
      # require a non-empty bin fastq
      [ -s "$binfq" ] || { warn "[run_emu_all] Missing/empty bin FASTQ for ${base} (${label}): $binfq"; continue; }

      # Respect enabled markers
      if [[ "$label" =~ ^ITS ]] && [ "$ENABLE_ITS" -ne 1 ]; then
        log "Skipping ITS bin '${label}' for ${base} (ITS DB missing)."; continue
      fi
      if [[ "$label" =~ ^LSU ]] && [ "$ENABLE_LSU" -ne 1 ]; then
        log "Skipping LSU bin '${label}' for ${base} (LSU DB missing)."; continue
      fi

      run_emu_for_bin "$label" "$binfq" "$base"
    done < "$BIN_TSV"
  done
}

# ----------------------------------------------------------
# Count reads quickly (gz or plain). Uses seqkit if available, else wc/4.
# ----------------------------------------------------------
count_reads() {
  local f="$1"
  [ -s "$f" ] || { echo 0; return 0; }
  if command -v seqkit >/dev/null 2>&1; then
    seqkit stats -T "$f" 2>/dev/null | awk 'NR==2{print $4+0; exit}'
  else
    local n
    if [[ "$f" == *.gz ]]; then
      n=$(gzip -cd "$f" | wc -l)
    else
      n=$(wc -l < "$f")
    fi
    echo $(( n / 4 ))
  fi
}

# ----------------------------------------------------------
# Build bin counts manifest (metadata/bin_counts.tsv)
# ----------------------------------------------------------
make_bin_counts_manifest() {
  log "Building bin counts manifest ..."
  mkdir -p "$METADIR"

  local out="${METADIR}/bin_counts.tsv"
  echo -e "library\tfiltered_fastq\ttotal_reads\tbin_label\tbin_fastq\tbin_reads\tunbinned_reads\tcheck_sum_ok" > "$out"

  for f in "${FASTQS[@]}"; do
    local base
    base=$(basename "$f")
    base=${base%.fastq.gz}; base=${base%.fq.gz}; base=${base%.fastq}; base=${base%.fq}

    local filtered="$f"
    local total
    total=$(count_reads "$filtered")

    # bins live under ${BINS_ROOT}/${base}/bins/
    local bin_dir="${BINS_ROOT}/${base}/bins"
    local sum_bins=0

    shopt -s nullglob
    local bin_files=( "$bin_dir"/*.fastq.gz "$bin_dir"/*.fq.gz "$bin_dir"/*.fastq "$bin_dir"/*.fq )
    shopt -u nullglob

    if [ ${#bin_files[@]} -gt 0 ]; then
      for bf in "${bin_files[@]}"; do
        local b_reads
        b_reads=$(count_reads "$bf")
        sum_bins=$(( sum_bins + b_reads ))
        local blabel
        blabel=$(basename "$bf")
        blabel=${blabel%.fastq.gz}; blabel=${blabel%.fq.gz}; blabel=${blabel%.fastq}; blabel=${blabel%.fq}
        blabel=${blabel#${base}.}
        echo -e "${base}\t${filtered}\t${total}\t${blabel}\t${bf}\t${b_reads}\tNA\tNA" >> "$out"
      done
    fi

    local unbinned=$(( total - sum_bins ))
    local ok="FALSE"
    if [ "$sum_bins" -eq "$total" ]; then ok="TRUE"; fi
    echo -e "${base}\t${filtered}\t${total}\tUNBINNED\tNA\t0\t${unbinned}\t${ok}" >> "$out"
  done

  log "Wrote bin counts manifest: ${out}"
}

# ----------------------------------------------------------
# Replicates strategy note (I / II / III)
# ----------------------------------------------------------
write_replicate_guidance() {
  cat > "${LOGDIR}/REPLICATE_STRATEGY.txt" <<EOF
Replicates strategy:
- Each FASTQ (I/II/III) was classified independently.
- For ecological inference, keep replicates separate when computing beta/alpha diversity.
- For per-site summaries, aggregate replicate abundances (median of compositional vectors
  or sum of estimated counts followed by renormalization).
EOF
  log "Wrote replicate guidance to ${LOGDIR}/REPLICATE_STRATEGY.txt"
}

# ----------------------------------------------------------
# Final report
# ----------------------------------------------------------
log_run_report() {
  {
    echo "=========================================================="
    echo " Emu Amplicons Run Report - $(date)"
    echo " Env: ${ENV_NAME} | Threads: ${THREADS}"
    echo " Enabled markers: 16S=${ENABLE_16S} ITS=${ENABLE_ITS} LSU=${ENABLE_LSU}"
    echo " DB16S: ${EMU_DB16S_DIR}"
    echo " DB_ITS: ${EMU_DB_ITS_DIR:-<none>}"
    echo " DB_LSU: ${EMU_DB_LSU_DIR:-<none>}"
    echo " Length bins:"
    echo "  ITS: ${ITS_MIN}-${ITS_MAX} | LSU: ${LSU_MIN}-${LSU_MAX} | 16S: ${S16_MIN}-${S16_MAX} | Ambiguous: ${AMB_MIN}-${AMB_MAX}"
    echo " FASTQ manifest: ${METADIR}/fastq_meta.tsv"
    echo " Log file: ${RUN_LOG}"
    echo "=========================================================="
  } | tee -a "$RUN_LOG" > "${LOGDIR}/RUN_REPORT.txt"
}

# ==========================================================
# Calling the functions
# ==========================================================
START_TIME=$(date +%s)

# Reset per-run timing file so reports don't accumulate across runs
[ -d "$LOGDIR" ] || mkdir -p "$LOGDIR"
: > "${LOGDIR}/.timing.tsv"

# Environment setup -------------------------------------------------------
time_function create_env_emu
time_function ensure_emu16s_db
time_function export_emu_env

# Optional DB presence checks (ITS/LSU) ----------------------------------
time_function 'check_optional_db "'"$EMU_DB_ITS_DIR"'" "ITS"' || true
time_function 'check_optional_db "'"$EMU_DB_LSU_DIR"'" "LSU"' || true

#Decide which markers are enabled (16S always; ITS/LSU if DBs valid)
time_function set_marker_flags

# Input discovery ---------------------------------------------------------
time_function discover_fastqs
time_function build_fastq_meta

# Binning by read length --------------------------------------------------
time_function split_bins_all

# Per-library/bin read counts manifest -----------------------------------
time_function make_bin_counts_manifest

# Bin-length boxplots (R) -------------------------------------------------
if command -v Rscript >/dev/null 2>&1 && [ -s "workflow/bin_length_boxplots.R" ]; then
  time_function "Rscript workflow/bin_length_boxplots.R --bins-root=${BINS_ROOT}"
else
  warn "[skip] Rscript or workflow/bin_length_boxplots.R not found; skipping bin length boxplots."
fi

#Emu runs per bin ------------------------------------------------------
time_function run_emu_all

# Replicate guidance & final report --------------------------------------
time_function write_replicate_guidance
time_function log_run_report

log "All done."
