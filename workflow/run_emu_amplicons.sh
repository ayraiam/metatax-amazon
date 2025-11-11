#!/usr/bin/env bash
# ==========================================================
# Script: run_emu_amplicons.sh
# Purpose: Install/check Emu + DBs, run Emu per FASTQ (no length bins),
#          collect abundance + mapping stats, and plot genus-level stacks.
# Notes:
#   - Requires conda/mamba (Miniforge/Conda) + seqkit
#   - Input FASTQ files live in results/filtered/ (default) or FASTQ_DIR_DEFAULT
#   - Outputs to results/, logs to logs/
# ==========================================================

set -euo pipefail

# ---------------- User-tunable defaults -------------------
THREADS="${THREADS:-8}"

# Databases:
EMU_DB16S_DIR="${EMU_DB16S_DIR:-$PWD/refdb/emu_16S}"

EMU_DB_ITS_DIR="${EMU_DB_ITS_DIR:-}"
EMU_DB_LSU_DIR="${EMU_DB_LSU_DIR:-}"

ENV_NAME="${ENV_NAME:-emu-env}"
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

log()  { echo -e ">>> $*" | tee -a "$RUN_LOG" >&2; }
warn() { echo -e "!!! $*" | tee -a "$RUN_LOG" >&2; }
die()  { echo -e "xxx $*" | tee -a "$RUN_LOG" >&2; exit 1; }

# ----------------------------------------------------------
# Function: ensure_channels
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
# ----------------------------------------------------------
create_env_emu() {
  log "Checking conda/mamba..."
  if ! command -v conda >/dev/null 2>&1; then
    die "conda not found. Install Miniforge/Conda first."
  fi
  source "$(conda info --base)/etc/profile.d/conda.sh"

  if conda env list | grep -qE "^${ENV_NAME}\s"; then
    log "Environment '${ENV_NAME}' already exists. Activating it..."
    conda activate "${ENV_NAME}"
    log "Environment '${ENV_NAME}' is now active."
  else
    ensure_channels
    log "Environment '${ENV_NAME}' not found. Creating it now (strict priority)..."
    set +e
    if command -v mamba >/dev/null 2>&1; then
      mamba create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        python=3.11 emu minimap2 "seqkit>=2.6" \
        "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" \
        pandas>=1.5 -y
    else
      conda create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        python=3.11 emu minimap2 "seqkit>=2.6" \
        "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" \
        pandas>=1.5 -y
    fi
    status=$?
    set -e

    if [ $status -ne 0 ]; then
      echo "!!! Strict-priority solve failed. Retrying once with flexible channel priority..." | tee -a "$RUN_LOG"
      if command -v mamba >/dev/null 2>&1; then
        mamba create -n "${ENV_NAME}" \
          -c conda-forge -c bioconda \
          python=3.11 emu minimap2 "seqkit>=2.6" \
          "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" \
          pandas>=1.5 -y
      else
        conda create -n "${ENV_NAME}" \
          -c conda-forge -c bioconda \
          python=3.11 emu minimap2 "seqkit>=2.6" \
          "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" \
          pandas>=1.5 -y
      fi
      echo ">>> Created '${ENV_NAME}' with flexible priority." | tee -a "$RUN_LOG"
    else
      echo ">>> Environment '${ENV_NAME}' created successfully (strict priority)." | tee -a "$RUN_LOG"
    fi

    conda activate "${ENV_NAME}"
    echo ">>> Environment '${ENV_NAME}' is now active." | tee -a "$RUN_LOG"
  fi

  log "Rscript path in env: $(which Rscript)"
  command -v emu >/dev/null 2>&1 || die "emu not available in env"
  command -v minimap2 >/dev/null 2>&1 || die "minimap2 not available in env"
  command -v seqkit >/dev/null 2>&1 || die "seqkit not available in env"
  log "Env ready: $(which python); emu=$(which emu)"
}

# ----------------------------------------------------------
# Export env
# ----------------------------------------------------------
export_emu_env() {
  mkdir -p envs
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
# Optional DB presence checks (not used here)
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
# FASTQ discovery (expects trimmed/filtered fastqs in results/filtered)
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

# limit to an OFFSET + N window of FASTQs
limit_to_n_fastqs() {
  local OFFSET="${OFFSET_FASTQS:-0}"
  local N="${LIMIT_FASTQS:-0}"

  # apply offset first
  if [[ "$OFFSET" -gt 0 && ${#FASTQS[@]} -gt "$OFFSET" ]]; then
    FASTQS=( "${FASTQS[@]:$OFFSET}" )
  elif [[ "$OFFSET" -gt 0 ]]; then
    FASTQS=()
  fi

  # then cap to N
  if [[ "$N" -gt 0 && ${#FASTQS[@]} -gt "$N" ]]; then
    FASTQS=( "${FASTQS[@]:0:$N}" )
    log "Limiting to ${N} FASTQs starting at offset ${OFFSET}."
  else
    log "Using ${#FASTQS[@]} FASTQs starting at offset ${OFFSET}."
  fi

  printf ">>> Using FASTQs:\n" | tee -a "$RUN_LOG"
  printf "    %s\n" "${FASTQS[@]}" | tee -a "$RUN_LOG"
}

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
# QUICK read counter
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
# Run Emu per FASTQ (no binning), keep outputs
# ----------------------------------------------------------
run_emu_per_fastq() {
  for fq in "${FASTQS[@]}"; do
    local base
    base=$(basename "$fq"); base=${base%.fastq.gz}; base=${base%.fq.gz}; base=${base%.fastq}; base=${base%.fq}
    local outdir="${OUTDIR}/emu_runs/${base}"
    mkdir -p "$outdir"

    log "Running Emu (16S) on ${base} ..."
    emu abundance \
      --threads "$THREADS" \
      --db "$EMU_DB16S_DIR" \
      --output-dir "$outdir" \
      --keep-counts \
      --keep-read-assignments \
      --output-unclassified \
      --type "$EMU_TYPE" \
      "$fq" 2>>"$RUN_LOG" | tee -a "$RUN_LOG" || warn "Emu failed for ${base}"

    # write a small per-sample info file with total reads for convenience
    echo -e "file\ttotal_reads" > "${outdir}/input_reads.tsv"
    echo -e "${base}\t$(count_reads "$fq")" >> "${outdir}/input_reads.tsv"

    # --- normalize abundance file to ${outdir}/abundance.tsv ---------------
    # try to find abundance file if Emu named it differently / gzipped / CSV
    if [ ! -s "${outdir}/abundance.tsv" ]; then
      # look up to one level deeper (some tools create nested folders)
      cand=$(find "$outdir" -maxdepth 2 -type f \
               \( -iname "*abundance*.tsv" -o -iname "*abundance*.csv" -o -iname "*abundance*.tsv.gz" -o -iname "*abundance*.csv.gz" \) \
               | head -n1)
      if [ -n "$cand" ]; then
        case "$cand" in
          *.tsv.gz)    gzip -cd "$cand" > "${outdir}/abundance.tsv" ;;
          *.csv.gz)    gzip -cd "$cand" | tr ',' '\t' > "${outdir}/abundance.tsv" ;;
          *.csv)       tr ',' '\t' < "$cand" > "${outdir}/abundance.tsv" ;;
          *.tsv|*.txt) cp -f "$cand" "${outdir}/abundance.tsv" ;;
          *)           : ;;
        esac
      fi
    fi
    # ----------------------------------------------------------------------
  done
}

collate_with_python() {
  local runs="${OUTDIR}/emu_runs"
  local tables="${OUTDIR}/tables"
  local dicts="${METADIR}"
  local minprob="${MIN_ASSIGN_PROB:-0.5}"
  # Save JSON? 1=yes, 0=no (default no)
  local save_json="${SAVE_JSON:-0}"

  if [[ "$save_json" -eq 1 ]]; then
    python workflow/emu_collect.py \
      --runs-dir "$runs" \
      --outdir "$tables" \
      --dictdir "$dicts" \
      --min-prob "$minprob"
  else
    python workflow/emu_collect.py \
      --runs-dir "$runs" \
      --outdir "$tables" \
      --min-prob "$minprob" \
      --no-json
  fi
}

# ----------------------------------------------------------
# Call the R plotting script
# ----------------------------------------------------------
plot_genus_stacks() {
  local abund="${OUTDIR}/tables/abundance_combined.tsv"
  local outd="${OUTDIR}/plots"
  [ -s "$abund" ] || { warn "No abundance table at ${abund}; skipping plots."; return 0; }
  if command -v Rscript >/dev/null 2>&1; then
    Rscript workflow/plot_genus_stacks.R "$abund" "$outd"
    log "Genus stacks saved in ${outd}/emu_genus_stacks_*.png/pdf"
  else
    warn "Rscript not found in env; skipping plots."
  fi
}

# ----------------------------------------------------------
# Replicates strategy note
# ----------------------------------------------------------
write_replicate_guidance() {
  cat > "${LOGDIR}/REPLICATE_STRATEGY.txt" <<EOF
Replicates strategy:
- Each FASTQ (I/II/III) was classified independently against the 16S database.
- For per-site summaries, aggregate replicate abundances (centered log-ratio or
  median across compositional vectors) after converting to relative abundance.
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
    echo " FASTQ manifest: ${METADIR}/fastq_meta.tsv"
    echo " Log file: ${RUN_LOG}"
    echo "=========================================================="
  } | tee -a "$RUN_LOG" > "${LOGDIR}/RUN_REPORT.txt"
}

# ==========================================================
# Calling the functions
# ==========================================================
START_TIME=$(date +%s)
[ -d "$LOGDIR" ] || mkdir -p "$LOGDIR"
: > "${LOGDIR}/.timing.tsv"

# Environment setup -------------------------------------------------------
time_function create_env_emu
time_function ensure_emu16s_db
time_function export_emu_env

# Optional DB presence checks (ITS/LSU) ----------------------------------
time_function 'check_optional_db "'"$EMU_DB_ITS_DIR"'" "ITS"' || true
time_function 'check_optional_db "'"$EMU_DB_LSU_DIR"'" "LSU"' || true
time_function set_marker_flags

# Input discovery ---------------------------------------------------------
time_function discover_fastqs
# limit to 3 for this test
time_function limit_to_n_fastqs
time_function build_fastq_meta

# Emu runs (per FASTQ) ----------------------------------------------------
time_function run_emu_per_fastq

# Collate + plot ----------------------------------------------------------
time_function collate_with_python
time_function plot_genus_stacks

# Guidance & final report -------------------------------------------------
time_function write_replicate_guidance
time_function log_run_report

log "All done."
