#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"

EMU_DB16S_DIR="${EMU_DB16S_DIR:-$PWD/refdb/emu_16S}"
EMU_DB_ITS_DIR="${EMU_DB_ITS_DIR:-}"
EMU_DB_LSU_DIR="${EMU_DB_LSU_DIR:-}"

ENV_NAME="${ENV_NAME:-emu-env}"
EMU_TYPE="${EMU_TYPE:-map-ont}"

LOGDIR="logs"
OUTDIR="results"
METADIR="metadata"
mkdir -p "$OUTDIR" "$METADIR" refdb "$LOGDIR"

RUN_LOG="${LOGDIR}/run_emu_amplicons.$(date +%Y%m%d_%H%M%S).log"
touch "$RUN_LOG"

ENABLE_16S=1
ENABLE_ITS=0
ENABLE_LSU=0

# ---------- BATCH SCOPING ----------
# Always get a deterministic batch tag like b050_n025
: "${OFFSET_FASTQS:=0}"
: "${LIMIT_FASTQS:=0}"
printf -v _OFF "%03d" "$OFFSET_FASTQS"
printf -v _LIM "%03d" "$LIMIT_FASTQS"
BATCH_TAG="${BATCH_TAG:-b${_OFF}_n${_LIM}}"

RUNS_DIR="${RUNS_DIR:-${OUTDIR}/emu_runs_${BATCH_TAG}}"
TABLES_DIR="${TABLES_DIR:-${OUTDIR}/tables_${BATCH_TAG}}"
PLOTS_DIR="${PLOTS_DIR:-${OUTDIR}/plots_${BATCH_TAG}}"
mkdir -p "$RUNS_DIR" "$TABLES_DIR" "$PLOTS_DIR"

SAVE_ASSIGN="${SAVE_ASSIGN:-0}"
SKIP_ASSIGN="${SKIP_ASSIGN:-$((1 - SAVE_ASSIGN))}"

time_function() {
  local fn="$1"; local start=$(date +%s); echo ">>> Running $fn ..." | tee -a "$RUN_LOG"
  local status=0; set +e; $fn; status=$?; set -e
  local end=$(date +%s); local dur=$((end - start))
  echo -e "${fn}\t${dur}" >> "${LOGDIR}/.timing.tsv"
  echo ">>> $fn completed in ${dur}s" | tee -a "$RUN_LOG"
  return $status
}
log()  { echo -e ">>> $*" | tee -a "$RUN_LOG" >&2; }
warn() { echo -e "!!! $*" | tee -a "$RUN_LOG" >&2; }
die()  { echo -e "xxx $*" | tee -a "$RUN_LOG" >&2; exit 1; }

ensure_channels() {
  conda config --remove-key channels >/dev/null 2>&1 || true
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda config --add channels defaults
  conda config --set channel_priority strict
  conda clean --index-cache -y >/dev/null 2>&1 || true
}

create_env_emu() {
  if ! command -v conda >/dev/null 2>&1; then die "conda not found."; fi
  source "$(conda info --base)/etc/profile.d/conda.sh"
  if conda env list | grep -qE "^${ENV_NAME}\s"; then
    conda activate "${ENV_NAME}"
  else
    ensure_channels
    set +e
    if command -v mamba >/dev/null 2>&1; then M=mamba; else M=conda; fi
    $M create -n "${ENV_NAME}" -c conda-forge -c bioconda \
      python=3.11 emu minimap2 "seqkit>=2.6" "r-base>=4.3" \
      "r-ggplot2>=3.4" "r-data.table" pandas>=1.5 -y || \
    $M create -n "${ENV_NAME}" -c conda-forge -c bioconda \
      python=3.11 emu minimap2 "seqkit>=2.6" "r-base>=4.3" \
      "r-ggplot2>=3.4" "r-data.table" pandas>=1.5 -y
    set -e
    conda activate "${ENV_NAME}"
  fi
  command -v emu >/dev/null 2>&1 || die "emu not available"
  command -v minimap2 >/dev/null 2>&1 || die "minimap2 not available"
  command -v seqkit >/dev/null 2>&1 || die "seqkit not available"
}

export_emu_env() {
  mkdir -p envs
  conda env export --name "${ENV_NAME}" > "envs/${ENV_NAME}.yml"
}

ensure_emu16s_db() {
  if [ -s "${EMU_DB16S_DIR}/species_taxid.fasta" ] && [ -s "${EMU_DB16S_DIR}/taxonomy.tsv" ]; then return 0; fi
  mkdir -p "${EMU_DB16S_DIR}"
  pip install --quiet osfclient
  export EMU_DATABASE_DIR="${EMU_DB16S_DIR}"
  pushd "${EMU_DB16S_DIR}" >/dev/null
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
    tar -xvf emu.tar
  popd >/dev/null
  [ -s "${EMU_DB16S_DIR}/species_taxid.fasta" ] || die "Failed 16S DB"
  [ -s "${EMU_DB16S_DIR}/taxonomy.tsv" ] || die "Failed taxonomy"
}

check_optional_db() {
  local d="$1" ; local label="$2"
  if [ -n "${d}" ] && [ -s "${d}/species_taxid.fasta" ] && [ -s "${d}/taxonomy.tsv" ]; then
    log "${label} DB OK: ${d}"; return 0
  fi
  warn "No usable ${label} DB at '${d:-<empty>}' — skipping ${label}."
  return 1
}

set_marker_flags() {
  ENABLE_16S=1; ENABLE_ITS=0; ENABLE_LSU=0
  if [ -n "${EMU_DB_ITS_DIR}" ] && [ -s "${EMU_DB_ITS_DIR}/species_taxid.fasta" ] && [ -s "${EMU_DB_ITS_DIR}/taxonomy.tsv" ]; then ENABLE_ITS=1; fi
  if [ -n "${EMU_DB_LSU_DIR}" ] && [ -s "${EMU_DB_LSU_DIR}/species_taxid.fasta" ] && [ -s "${EMU_DB_LSU_DIR}/taxonomy.tsv" ]; then ENABLE_LSU=1; fi
  log "Markers enabled -> 16S:${ENABLE_16S} ITS:${ENABLE_ITS} LSU:${ENABLE_LSU}"
}

FASTQ_DIR_DEFAULT="${FASTQ_DIR_DEFAULT:-results/filtered}"
FASTQ_MANIFEST="${FASTQ_MANIFEST:-${METADIR}/fastq_meta.${BATCH_TAG}.tsv}"

discover_fastqs() {
  log "Scanning ${FASTQ_DIR_DEFAULT}/ for filtered FASTQ files..."
  shopt -s nullglob
  FASTQS=( "${FASTQ_DIR_DEFAULT}"/*.fastq.gz "${FASTQ_DIR_DEFAULT}"/*.fq.gz \
           "${FASTQ_DIR_DEFAULT}"/*.fastq     "${FASTQ_DIR_DEFAULT}"/*.fq )
  shopt -u nullglob
  [ ${#FASTQS[@]} -gt 0 ] || die "No FASTQ files found in ${FASTQ_DIR_DEFAULT}/"
  log "Found ${#FASTQS[@]} FASTQ files."
}

limit_to_n_fastqs() {
  local OFFSET="${OFFSET_FASTQS:-0}"
  local N="${LIMIT_FASTQS:-0}"
  if [[ "$OFFSET" -gt 0 && ${#FASTQS[@]} -gt "$OFFSET" ]]; then FASTQS=( "${FASTQS[@]:$OFFSET}" ); elif [[ "$OFFSET" -gt 0 ]]; then FASTQS=(); fi
  if [[ "$N" -gt 0 && ${#FASTQS[@]} -gt "$N" ]]; then FASTQS=( "${FASTQS[@]:0:$N}" ); fi
  printf ">>> Using FASTQs:\n" | tee -a "$RUN_LOG"; printf "    %s\n" "${FASTQS[@]}" | tee -a "$RUN_LOG"
}

build_fastq_meta() {
  if [ -s "$FASTQ_MANIFEST" ]; then log "Reusing manifest: $FASTQ_MANIFEST"; return 0; fi
  mkdir -p "$(dirname "$FASTQ_MANIFEST")"
  {
    echo -e "file\tsample\treplicate\tmolecular_feature"
    for f in "${FASTQS[@]}"; do
      b=$(basename "$f"); b=${b%.fastq.gz}; b=${b%.fq.gz}; b=${b%.fastq}; b=${b%.fq}
      IFS='_' read -r sample replicate molecular_feature _ <<< "$b"
      echo -e "${f}\t${sample:-NA}\t${replicate:-NA}\t${molecular_feature:-NA}"
    done
  } > "$FASTQ_MANIFEST"
  log "Wrote manifest: $FASTQ_MANIFEST"
}

count_reads() {
  local f="$1"; [ -s "$f" ] || { echo 0; return 0; }
  if command -v seqkit >/dev/null 2>&1; then
    seqkit stats -T "$f" 2>/dev/null | awk 'NR==2{print $4+0; exit}'
  else
    local n; if [[ "$f" == *.gz ]]; then n=$(gzip -cd "$f" | wc -l); else n=$(wc -l < "$f"); fi
    echo $(( n / 4 ))
  fi
}

run_emu_per_fastq() {
  for fq in "${FASTQS[@]}"; do
    local base; base=$(basename "$fq"); base=${base%.fastq.gz}; base=${base%.fq.gz}; base=${base%.fastq}; base=${base%.fq}
    local outdir="${RUNS_DIR}/${base}"
    mkdir -p "$outdir"
    log "Running Emu (16S) on ${base} ..."
    if [[ "${SAVE_ASSIGN}" -eq 1 ]]; then
      emu abundance \
        --threads "$THREADS" \
        --db "$EMU_DB16S_DIR" \
        --output-dir "$outdir" \
        --keep-counts \
        --keep-read-assignments \
        --output-unclassified \
        --type "$EMU_TYPE" \
        "$fq" 2>>"$RUN_LOG" | tee -a "$RUN_LOG" || warn "Emu failed for ${base}"
    else
      emu abundance \
        --threads "$THREADS" \
        --db "$EMU_DB16S_DIR" \
        --output-dir "$outdir" \
        --keep-counts \
        --output-unclassified \
        --type "$EMU_TYPE" \
        "$fq" 2>>"$RUN_LOG" | tee -a "$RUN_LOG" || warn "Emu failed for ${base}"
    fi
    echo -e "file\ttotal_reads" > "${outdir}/input_reads.tsv"
    echo -e "${base}\t$(count_reads "$fq")" >> "${outdir}/input_reads.tsv"

    if [ ! -s "${outdir}/abundance.tsv" ]; then
      cand=$(find "$outdir" -maxdepth 2 -type f \( -iname "*abundance*.tsv" -o -iname "*abundance*.csv" -o -iname "*abundance*.tsv.gz" -o -iname "*abundance*.csv.gz" \) | head -n1)
      if [ -n "$cand" ]; then
        case "$cand" in
          *.tsv.gz) gzip -cd "$cand" > "${outdir}/abundance.tsv" ;;
          *.csv.gz) gzip -cd "$cand" | tr ',' '\t' > "${outdir}/abundance.tsv" ;;
          *.csv)    tr ',' '\t' < "$cand" > "${outdir}/abundance.tsv" ;;
          *.tsv|*.txt) cp -f "$cand" "${outdir}/abundance.tsv" ;;
        esac
      fi
    fi
  done
}

collate_with_python() {
  local runs="${RUNS_DIR}"
  local tables="${TABLES_DIR}"
  local dicts="${METADIR}"
  local minprob="${MIN_ASSIGN_PROB:-0.5}"
  local save_json="${SAVE_JSON:-0}"
  mkdir -p "$tables"
  local SKIP_FLAG=""
  [[ "${SKIP_ASSIGN}" -eq 1 ]] && SKIP_FLAG="--skip-assign"
  if [[ "$save_json" -eq 1 ]]; then
    python workflow/emu_collect.py --runs-dir "$runs" --outdir "$tables" --dictdir "$dicts" --min-prob "$minprob" $SKIP_FLAG
  else
    python workflow/emu_collect.py --runs-dir "$runs" --outdir "$tables" --min-prob "$minprob" --no-json $SKIP_FLAG
  fi
}

plot_genus_stacks() {
  local abund="${TABLES_DIR}/abundance_combined.tsv"
  local outd="${PLOTS_DIR}"
  [ -s "$abund" ] || { warn "No abundance table at ${abund}; skipping plots."; return 0; }
  if command -v Rscript >/dev/null 2>&1; then
    Rscript workflow/plot_genus_stacks.R "$abund" "$outd"
  else
    warn "Rscript not found; skipping plots."
  fi
}

write_replicate_guidance() {
  cat > "${LOGDIR}/REPLICATE_STRATEGY.txt" <<EOF
Replicates strategy:
- Each FASTQ (I/II/III) classified independently against the 16S DB.
- For per-site summaries, aggregate replicate abundances after converting to relative abundance.
EOF
}

log_run_report() {
  {
    echo "=========================================================="
    echo " Emu Amplicons Run Report - $(date)"
    echo " Batch tag      : ${BATCH_TAG}"
    echo " Runs dir       : ${RUNS_DIR}"
    echo " Tables dir     : ${TABLES_DIR}"
    echo " Plots dir      : ${PLOTS_DIR}"
    echo "=========================================================="
  } > "${LOGDIR}/RUN_REPORT.txt"
}

START_TIME=$(date +%s)
[ -d "$LOGDIR" ] || mkdir -p "$LOGDIR"
: > "${LOGDIR}/.timing.tsv"

time_function create_env_emu
time_function ensure_emu16s_db
time_function export_emu_env
time_function 'check_optional_db "'"$EMU_DB_ITS_DIR"'" "ITS"' || true
time_function 'check_optional_db "'"$EMU_DB_LSU_DIR"'" "LSU"' || true
time_function set_marker_flags
time_function discover_fastqs
time_function limit_to_n_fastqs
time_function build_fastq_meta
time_function run_emu_per_fastq
time_function collate_with_python
time_function plot_genus_stacks
time_function write_replicate_guidance
time_function log_run_report

log "All done."
