#!/usr/bin/env bash
# ==========================================================
# Script: workflow/runall.sh
# Purpose: Submit libsQC run via Slurm (srun) and optionally
#          run Emu Amplicons afterward.
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

EMU_PARTITION=""
EMU_TIME=""
EMU_CPUS=""
EMU_MEM=""
EMU_DB_ITS_DIR=""
EMU_DB_LSU_DIR=""

# Batching controls
OFFSET_FASTQS="${OFFSET_FASTQS:-0}"
# (LIMIT_FASTQS is read in run_emu_amplicons.sh via env)

# derive a stable batch tag from OFFSET/LIMIT for per-batch outputs
LIMIT_FASTQS_SHOW="${LIMIT_FASTQS:-0}"                                # for echo only
LIMIT_SAFE="${LIMIT_FASTQS:-0}"
printf -v OFFSET_PAD "%03d" "${OFFSET_FASTQS}"
printf -v LIMIT_PAD  "%03d" "${LIMIT_SAFE}"
BATCH_TAG="${BATCH_TAG:-b${OFFSET_PAD}_n${LIMIT_PAD}}"                # e.g., b050_n026

# per-batch output roots (used by Emu step)
RUNS_DIR_DEFAULT="results/emu_runs_${BATCH_TAG}"
TABLES_DIR_DEFAULT="results/tables_${BATCH_TAG}"
PLOTS_DIR_DEFAULT="results/plots_${BATCH_TAG}"

# --- Help message ---
usage() {
  echo "Usage: bash workflow/runall.sh [options]"
  echo
  echo "General:"
  echo "  --partition STR       Partition/queue name (default: short)"
  echo "  --time HH:MM:SS       Walltime (default: 04:00:00)"
  echo "  --cpus INT            CPUs per task (default: 4)"
  echo "  --mem STR             Memory (default: 16G)"
  echo "  --wd PATH             Working directory (default: current directory)"
  echo "  --primer-fwd SEQ      Forward primer sequence (optional)"
  echo "  --primer-rev SEQ      Reverse primer sequence (optional)"
  echo "  --seq-summary PATH    Path to sequencing_summary.txt (optional)"
  echo
  echo "Batching:"
  echo "  --offset-fastqs INT   Skip this many FASTQs from the start before selecting (default: 0)"
  echo "                        (Set LIMIT_FASTQS=N via env to cap how many to run.)"
  echo
  echo "Stage control:"
  echo "  --no-qc               Skip libsQC step (run Emu only)"
  echo "  --no-emu              Skip Emu step (run libsQC only)"
  echo
  echo "Emu options:"
  echo "  --emu-partition STR   Emu partition (default: inherit libsQC)"
  echo "  --emu-time HH:MM:SS   Emu walltime (default: inherit libsQC)"
  echo "  --emu-cpus INT        Emu CPUs (default: inherit libsQC)"
  echo "  --emu-mem STR         Emu memory (default: inherit libsQC)"
  echo "  --emu-db-its PATH     Emu ITS DB directory (optional)"
  echo "  --emu-db-lsu PATH     Emu LSU DB directory (optional)"
  echo
  echo "  -h, --help            Show this help message"
  exit 0
}

# --- Parse CLI arguments ---
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
    --offset-fastqs) OFFSET_FASTQS="$2"; shift 2 ;;
    --no-emu) RUN_EMU=0; shift 1 ;;
    --no-qc|--skip-libsQC) RUN_LIBSQC=0; shift 1 ;;
    --emu-partition) EMU_PARTITION="$2"; shift 2 ;;
    --emu-time) EMU_TIME="$2"; shift 2 ;;
    --emu-cpus) EMU_CPUS="$2"; shift 2 ;;
    --emu-mem) EMU_MEM="$2"; shift 2 ;;
    --emu-db-its) EMU_DB_ITS_DIR="$2"; shift 2 ;;
    --emu-db-lsu) EMU_DB_LSU_DIR="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

EMU_PARTITION="${EMU_PARTITION:-$PARTITION}"
EMU_TIME="${EMU_TIME:-$TIME}"
EMU_CPUS="${EMU_CPUS:-$CPUS}"
EMU_MEM="${EMU_MEM:-$MEM}"

mkdir -p logs

# write default multi-primer lists and expose as files to QC >>>
mkdir -p metadata
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
# If single primers provided via CLI, append them as well (uppercased).
if [[ -n "$PRIMER_FWD" ]]; then echo "$PRIMER_FWD" | tr '[:lower:]' '[:upper:]' >> "$PRIMERS_FWD_FILE"; fi
if [[ -n "$PRIMER_REV" ]]; then echo "$PRIMER_REV" | tr '[:lower:]' '[:upper:]' >> "$PRIMERS_REV_FILE"; fi

TS=$(date +%Y%m%d_%H%M%S)
OUT_LOG="logs/run_${TS}.out"
ERR_LOG="logs/run_${TS}.err"
EMU_OUT_LOG="logs/emu_${TS}.out"
EMU_ERR_LOG="logs/emu_${TS}.err"

echo "============================================"
echo "libsQC:  ${RUN_LIBSQC}"
echo "Emu:     ${RUN_EMU}"
echo "--------------------------------------------"
echo "Partition : $PARTITION"
echo "Time      : $TIME"
echo "CPUs      : $CPUS"
echo "Memory    : $MEM"
echo "Work dir  : $WDIR"
echo "Offset    : $OFFSET_FASTQS"
echo "Limit     : ${LIMIT_FASTQS_SHOW} (env)"
echo "Batch tag : ${BATCH_TAG}"
echo "Runs dir  : ${RUNS_DIR_DEFAULT}"
echo "Tables dir: ${TABLES_DIR_DEFAULT}"
echo "Plots dir : ${PLOTS_DIR_DEFAULT}"
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

# --- Emu step ---
if [[ "$RUN_EMU" -eq 1 ]]; then
  echo ">>> Launching Emu Amplicons..."
  export OMP_NUM_THREADS="$EMU_CPUS"
  export MKL_NUM_THREADS="$EMU_CPUS"
  export NUMEXPR_NUM_THREADS="$EMU_CPUS"

  # default behavior—do NOT keep massive read-assignment matrices
  SAVE_ASSIGN="${SAVE_ASSIGN:-0}"            # 1 to keep Emu read-assign matrices
  SKIP_ASSIGN="${SKIP_ASSIGN:-$((1-SAVE_ASSIGN))}"  # collector will skip loading them

	srun \
    --partition="$EMU_PARTITION" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$EMU_CPUS" \
    --mem="$EMU_MEM" \
    --time="$EMU_TIME" \
    --chdir="$WDIR" \
    --export=ALL,THREADS="$EMU_CPUS",EMU_DB_ITS_DIR="$EMU_DB_ITS_DIR",EMU_DB_LSU_DIR="$EMU_DB_LSU_DIR",FASTQ_DIR_DEFAULT="${FASTQ_DIR_DEFAULT:-results/filtered}",LIMIT_FASTQS="${LIMIT_FASTQS:-0}",OFFSET_FASTQS="${OFFSET_FASTQS:-0}",BATCH_TAG="${BATCH_TAG}",RUNS_DIR="${RUNS_DIR_DEFAULT}",TABLES_DIR="${TABLES_DIR_DEFAULT}",PLOTS_DIR="${PLOTS_DIR_DEFAULT}",SAVE_JSON="${SAVE_JSON:-0}",SAVE_ASSIGN="${SAVE_ASSIGN}",SKIP_ASSIGN="${SKIP_ASSIGN}" \
    /bin/bash workflow/run_emu_amplicons.sh \
    1>"$EMU_OUT_LOG" \
    2>"$EMU_ERR_LOG"
else
  echo ">>> Skipping Emu (--no-emu)"
fi

echo
echo ">>> Jobs finished. Check logs:"
[ "$RUN_LIBSQC" -eq 1 ] && echo "  $OUT_LOG / $ERR_LOG"
[ "$RUN_EMU" -eq 1 ] && echo "  $EMU_OUT_LOG / $EMU_ERR_LOG"
