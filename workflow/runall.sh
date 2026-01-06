#!/usr/bin/env bash
# ==========================================================
# Script: workflow/runall.sh
# Purpose: Submit libsQC via Slurm and optionally run Emu
#          and downstream diversity analysis.
# ==========================================================
set -euo pipefail

# --- Default values ---
PARTITION="short"
TIME="04:00:00"
CPUS="4"
MEM="16G"
WDIR="$PWD"
PRIMER_FWD=""
PRIMER_REV=""
SEQ_SUMMARY=""

RUN_EMU=1
RUN_LIBSQC=1
RUN_DOWNSTREAM=1
BUILD_MARKER_DBS=1
ONLY_BUILD_MARKER_DBS=0

EMU_PARTITION=""
EMU_TIME=""
EMU_CPUS=""
EMU_MEM=""
EMU_DB_ITS_DIR=""
EMU_DB_LSU_DIR=""

DOWNSTREAM_PARTITION=""
DOWNSTREAM_TIME=""
DOWNSTREAM_CPUS=""
DOWNSTREAM_MEM=""
DOWNSTREAM_ENV_NAME="emu-downstream"

### downstream input file becomes MODE-dependent (16S vs ITS)
MODE="16S"          ### 16S or ITS
USE_COUNTS_0_4="0"  # Parts 0-4: abundance
USE_COUNTS_5="1"    # Part 5: estimated_counts

# Dedicated flag for ANCOM-BC2 input (default: estimated_counts)
USE_COUNTS_ANCOM="1"  # Steps 7-8 (ANCOM-BC2): 1=estimated_counts, 0=abundance

DOWNSTREAM_INFILE_16S="/home/t.sousa/metataxonomy_rds/metatax-amazon/results/tables/abundance_combined.tsv"
DOWNSTREAM_INFILE_ITS="/home/t.sousa/metataxonomy_rds/metatax-amazon/results/tables_ITS/abundance_combined.tsv"
DOWNSTREAM_INFILE="$DOWNSTREAM_INFILE_16S"

DOWNSTREAM_OUTDIR="/home/t.sousa/metataxonomy_rds/metatax-amazon/results/plots"
DOWNSTREAM_BASENAME="downstream"

usage() {
  echo "Usage: bash workflow/runall.sh [options]"
  echo
  echo "General:"
  echo "  --partition STR       Partition/queue (default: short)"
  echo "  --time HH:MM:SS       Walltime (default: 04:00:00)"
  echo "  --cpus INT            CPUs per task (default: 4)"
  echo "  --mem STR             Memory (default: 16G)"
  echo "  --wd PATH             Working directory (default: current dir)"
  echo "  --primer-fwd SEQ      Forward primer (optional)"
  echo "  --primer-rev SEQ      Reverse primer (optional)"
  echo "  --seq-summary PATH    sequencing_summary.txt (optional)"
  echo
  echo "Stage control:"
  echo "  --no-qc               Skip libsQC"
  echo "  --no-emu              Skip Emu"
  echo "  --no-downstream       Skip downstream analysis"
  echo "  --no-build-marker-dbs Do not attempt to build ITS/LSU marker databases"
  echo "  --only-build-marker-dbs  Build ITS/LSU DBs ONLY (no QC, no Emu, no downstream)"
  echo
  echo "Emu options:"
  echo "  --emu-partition STR   Emu partition (default: inherit libsQC)"
  echo "  --emu-time HH:MM:SS   Emu walltime (default: inherit libsQC)"
  echo "  --emu-cpus INT        Emu CPUs (default: inherit libsQC)"
  echo "  --emu-mem STR         Emu memory (default: inherit libsQC)"
  echo "  --emu-db-its PATH     Emu ITS DB directory (optional)"
  echo "  --emu-db-lsu PATH     Emu LSU DB directory (optional)"
  echo
  echo "Downstream options:"
  echo "  --downstream-partition STR  Partition (default: inherit libsQC)"
  echo "  --downstream-time HH:MM:SS  Walltime (default: inherit libsQC)"
  echo "  --downstream-cpus INT       CPUs (default: inherit libsQC)"
  echo "  --downstream-mem STR        Memory (default: inherit libsQC)"
  echo "  --downstream-infile PATH    Input abundance table (optional override)"
  echo "  --downstream-outdir PATH    Output plots directory"
  echo "  --downstream-env STR        Conda env name (default: emu-downstream)"
  echo "  --downstream-basename STR   Basename prefix for outputs"
  echo
  echo "Mode options (kingdom runs separately by directory):"
  echo "  --mode STR             16S or ITS (default: 16S)"
  echo "  --use-counts-0-4 INT   Parts 0–4: 1=estimated_counts, 0=abundance (default: 0)"
  echo "  --use-counts-5   INT   Part 5:   1=estimated_counts, 0=abundance (default: 1)"
  echo "  --use-counts-ancom INT Steps 7–8: 1=estimated_counts, 0=abundance (default: 1)"  
  echo
  echo "  -h, --help            Show this help message"
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --partition) PARTITION="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --wd) WDIR="$2"; shift 2 ;;
    --primer-fwd) PRIMER_FWD="$2"; shift 2 ;;
    --primer-rev) PRIMER_REV="$2"; shift 2 ;;
    --seq-summary) SEQ_SUMMARY="$2"; shift 2 ;;
    --no-emu) RUN_EMU=0; shift 1 ;;
    --no-qc|--skip-libsQC) RUN_LIBSQC=0; shift 1 ;;
    --no-downstream) RUN_DOWNSTREAM=0; shift 1 ;;
    --no-build-marker-dbs) BUILD_MARKER_DBS=0; shift 1 ;;
    --only-build-marker-dbs)
      RUN_LIBSQC=0
      RUN_EMU=0
      RUN_DOWNSTREAM=0
      BUILD_MARKER_DBS=1
      ONLY_BUILD_MARKER_DBS=1
      shift 1 ;;
    --emu-partition) EMU_PARTITION="$2"; shift 2 ;;
    --emu-time) EMU_TIME="$2"; shift 2 ;;
    --emu-cpus) EMU_CPUS="$2"; shift 2 ;;
    --emu-mem) EMU_MEM="$2"; shift 2 ;;
    --emu-db-its) EMU_DB_ITS_DIR="$2"; shift 2 ;;
    --emu-db-lsu) EMU_DB_LSU_DIR="$2"; shift 2 ;;
    --downstream-partition) DOWNSTREAM_PARTITION="$2"; shift 2 ;;
    --downstream-time) DOWNSTREAM_TIME="$2"; shift 2 ;;
    --downstream-cpus) DOWNSTREAM_CPUS="$2"; shift 2 ;;
    --downstream-mem) DOWNSTREAM_MEM="$2"; shift 2 ;;
    --downstream-infile) DOWNSTREAM_INFILE="$2"; shift 2 ;;  # user override still allowed
    --downstream-outdir) DOWNSTREAM_OUTDIR="$2"; shift 2 ;;
    --downstream-env) DOWNSTREAM_ENV_NAME="$2"; shift 2 ;;
    --downstream-basename) DOWNSTREAM_BASENAME="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --use-counts-0-4) USE_COUNTS_0_4="$2"; shift 2 ;;
    --use-counts-5)   USE_COUNTS_5="$2"; shift 2 ;;
    --use-counts-ancom) USE_COUNTS_ANCOM="$2"; shift 2 ;;  

    -h|--help) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

EMU_PARTITION="${EMU_PARTITION:-$PARTITION}"
EMU_TIME="${EMU_TIME:-$TIME}"
EMU_CPUS="${EMU_CPUS:-$CPUS}"
EMU_MEM="${EMU_MEM:-$MEM}"

DOWNSTREAM_PARTITION="${DOWNSTREAM_PARTITION:-$PARTITION}"
DOWNSTREAM_TIME="${DOWNSTREAM_TIME:-$TIME}"
DOWNSTREAM_CPUS="${DOWNSTREAM_CPUS:-$CPUS}"
DOWNSTREAM_MEM="${DOWNSTREAM_MEM:-$MEM}"

### if user did NOT override --downstream-infile explicitly,
### switch default infile based on MODE (16S vs ITS directory)
if [[ "${DOWNSTREAM_INFILE}" == "${DOWNSTREAM_INFILE_16S}" || "${DOWNSTREAM_INFILE}" == "${DOWNSTREAM_INFILE_ITS}" ]]; then
  if [[ "$(echo "$MODE" | tr '[:lower:]' '[:upper:]')" == "ITS" ]]; then
    DOWNSTREAM_INFILE="$DOWNSTREAM_INFILE_ITS"
  else
    DOWNSTREAM_INFILE="$DOWNSTREAM_INFILE_16S"
  fi
fi

mkdir -p logs metadata

# Primer lists ----------------------------
PRIMERS_FWD_FILE="metadata/primers_fwd.list"
PRIMERS_REV_FILE="metadata/primers_rev.list"
cat > "$PRIMERS_FWD_FILE" <<'EOF'
TTCCGGTTGATCCTGCCGGA
TCCGTAGGTGAACCTGCGG
AGAGTTTGATCMTGGCTCAG
TACTACCACCAAGATCT
EOF
cat > "$PRIMERS_REV_FILE" <<'EOF'
TACGGWTACCTTGTTACGACTT
TCCTCCGCTTATTGATATGC
GGTTACCTTGTTACGACTT
ACCCGCTGAACTTAAGC
EOF
if [[ -n "$PRIMER_FWD" ]]; then echo "$PRIMER_FWD" | tr '[:lower:]' '[:upper:]' >> "$PRIMERS_FWD_FILE"; fi
if [[ -n "$PRIMER_REV" ]]; then echo "$PRIMER_REV" | tr '[:lower:]' '[:upper:]' >> "$PRIMERS_REV_FILE"; fi

TS=$(date +%Y%m%d_%H%M%S)
OUT_LOG="logs/run_${TS}.out"
ERR_LOG="logs/run_${TS}.err"
EMU_OUT_LOG="logs/emu_${TS}.out"
EMU_ERR_LOG="logs/emu_${TS}.err"
DOWN_OUT_LOG="logs/downstream_${TS}.out"
DOWN_ERR_LOG="logs/downstream_${TS}.err"

echo "============================================"
echo "libsQC:       ${RUN_LIBSQC}"
echo "Emu:          ${RUN_EMU}"
echo "Downstream:   ${RUN_DOWNSTREAM}"
echo "--------------------------------------------"
echo "Partition : $PARTITION"
echo "Time      : $TIME"
echo "CPUs      : $CPUS"
echo "Memory    : $MEM"
echo "Work dir  : $WDIR"
echo "MODE      : $MODE"
echo "USE_COUNTS_0_4: $USE_COUNTS_0_4"
echo "USE_COUNTS_5  : $USE_COUNTS_5"
echo "USE_COUNTS_ANCOM: $USE_COUNTS_ANCOM"   
echo "DOWN_IN   : $DOWNSTREAM_INFILE"
echo "============================================"
echo

export OMP_NUM_THREADS="$CPUS"
export MKL_NUM_THREADS="$CPUS"
export NUMEXPR_NUM_THREADS="$CPUS"

# --- libsQC step ---
if [[ "$RUN_LIBSQC" -eq 1 ]]; then
  echo ">>> Launching libsQC..."
  srun \
    --partition="$PARTITION" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --time="$TIME" \
    --chdir="$WDIR" \
    --export=ALL,THREADS="$CPUS",PRIMER_FWD="$PRIMER_FWD",PRIMER_REV="$PRIMER_REV",SEQ_SUMMARY="$SEQ_SUMMARY",PRIMERS_FWD_FILE="$PRIMERS_FWD_FILE",PRIMERS_REV_FILE="$PRIMERS_REV_FILE" \
    /bin/bash workflow/run_libsQC.sh \
    1>"$OUT_LOG" \
    2>"$ERR_LOG"
else
  echo ">>> Skipping libsQC (--no-qc)"
fi

### ITS/LSU DB build step (UNITE ITS + SILVA LSU) --------------------
if [[ "$BUILD_MARKER_DBS" -eq 1 ]]; then
  echo ">>> Ensuring ITS/LSU Emu databases (UNITE ITS + SILVA LSU)..."
  srun \
    --partition="$PARTITION" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --time="$TIME" \
    --chdir="$WDIR" \
    --export=ALL,ENV_NAME="emu-env",EMU_DB_ITS_DIR="$EMU_DB_ITS_DIR",EMU_DB_LSU_DIR="$EMU_DB_LSU_DIR",ITS_FASTA="${ITS_FASTA:-}",LSU_FASTA="${LSU_FASTA:-}" \
    /bin/bash workflow/run_build_ITS_LSU_dbs.sh
else
  echo ">>> Skipping ITS/LSU DB build (--no-build-marker-dbs)"
fi
### ----------------------------------------------------------------------

# --- Emu step ---
if [[ "$RUN_EMU" -eq 1 ]]; then
  echo ">>> Launching Emu Amplicons..."
  export OMP_NUM_THREADS="$EMU_CPUS"
  export MKL_NUM_THREADS="$EMU_CPUS"
  export NUMEXPR_NUM_THREADS="$EMU_CPUS"

  # Build export list for Emu step (avoid clobbering DB defaults; pass marker flags if set)
  EMU_EXPORT="ALL,THREADS=$EMU_CPUS,FASTQ_DIR_DEFAULT=${FASTQ_DIR_DEFAULT:-results/filtered},LIMIT_FASTQS=${LIMIT_FASTQS:-1}"

  # Only pass custom DB paths if the user explicitly set them
  if [[ -n "$EMU_DB_ITS_DIR" ]]; then
    EMU_EXPORT+=",EMU_DB_ITS_DIR=$EMU_DB_ITS_DIR"
  fi
  if [[ -n "$EMU_DB_LSU_DIR" ]]; then
    EMU_EXPORT+=",EMU_DB_LSU_DIR=$EMU_DB_LSU_DIR"
  fi

  # Pass marker-selection flags if they exist in the environment
  if [[ -n "${ENABLE_16S:-}" ]]; then EMU_EXPORT+=",ENABLE_16S=$ENABLE_16S"; fi
  if [[ -n "${ENABLE_ITS:-}" ]]; then EMU_EXPORT+=",ENABLE_ITS=$ENABLE_ITS"; fi
  if [[ -n "${ENABLE_LSU:-}" ]]; then EMU_EXPORT+=",ENABLE_LSU=$ENABLE_LSU"; fi

  srun \
    --partition="$EMU_PARTITION" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$EMU_CPUS" \
    --mem="$EMU_MEM" \
    --time="$EMU_TIME" \
    --chdir="$WDIR" \
    --export="$EMU_EXPORT" \
    /bin/bash workflow/run_emu_amplicons.sh \
    1>"$EMU_OUT_LOG" \
    2>"$EMU_ERR_LOG"
else
  echo ">>> Skipping Emu (--no-emu)"
fi

# --- Downstream step ---
if [[ "$RUN_DOWNSTREAM" -eq 1 ]]; then
  echo ">>> Launching downstream analysis..."
  srun \
    --partition="$DOWNSTREAM_PARTITION" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$DOWNSTREAM_CPUS" \
    --mem="$DOWNSTREAM_MEM" \
    --time="$DOWNSTREAM_TIME" \
    --chdir="$WDIR" \
    --export=ALL,ENV_NAME="$DOWNSTREAM_ENV_NAME",INFILE="$DOWNSTREAM_INFILE",OUTDIR="$DOWNSTREAM_OUTDIR",BASENAME="$DOWNSTREAM_BASENAME",MODE="$MODE",USE_COUNTS_0_4="$USE_COUNTS_0_4",USE_COUNTS_5="$USE_COUNTS_5",USE_COUNTS_ANCOM="$USE_COUNTS_ANCOM" \
    /bin/bash workflow/run_downstream_analysis.sh \
    1>"$DOWN_OUT_LOG" \
    2>"$DOWN_ERR_LOG"
else
  echo ">>> Skipping downstream (--no-downstream)"
fi

echo
echo ">>> Jobs finished. Check logs:"
[ "$RUN_LIBSQC" -eq 1 ]     && echo "  $OUT_LOG / $ERR_LOG"
[ "$RUN_EMU" -eq 1 ]        && echo "  $EMU_OUT_LOG / $EMU_ERR_LOG"
[ "$RUN_DOWNSTREAM" -eq 1 ] && echo "  $DOWN_OUT_LOG / $DOWN_ERR_LOG"
