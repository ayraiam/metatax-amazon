#!/usr/bin/env bash
set -euo pipefail

PARTITION="short"; TIME="04:00:00"; CPUS="4"; MEM="16G"; WDIR="$PWD"
PRIMER_FWD=""; PRIMER_REV=""; SEQ_SUMMARY=""
RUN_EMU=1; RUN_LIBSQC=1
EMU_PARTITION=""; EMU_TIME=""; EMU_CPUS=""; EMU_MEM=""
EMU_DB_ITS_DIR=""; EMU_DB_LSU_DIR=""

OFFSET_FASTQS="${OFFSET_FASTQS:-0}"

# ---------- BATCH TAG ----------
: "${LIMIT_FASTQS:=0}"
printf -v OFFSET_PAD "%03d" "${OFFSET_FASTQS}"
printf -v LIMIT_PAD  "%03d" "${LIMIT_FASTQS}"
BATCH_TAG="${BATCH_TAG:-b${OFFSET_PAD}_n${LIMIT_PAD}}"

RUNS_DIR_DEFAULT="results/emu_runs_${BATCH_TAG}"
TABLES_DIR_DEFAULT="results/tables_${BATCH_TAG}"
PLOTS_DIR_DEFAULT="results/plots_${BATCH_TAG}"

usage(){ echo "Usage: bash workflow/runall.sh [options]"; exit 0; }

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
    --no-emu) RUN_EMU=0; shift ;;
    --no-qc|--skip-libsQC) RUN_LIBSQC=0; shift ;;
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

mkdir -p logs metadata

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
[[ -n "$PRIMER_FWD" ]] && echo "$PRIMER_FWD" | tr '[:lower:]' '[:upper:]' >> "$PRIMERS_FWD_FILE"
[[ -n "$PRIMER_REV" ]] && echo "$PRIMER_REV" | tr '[:lower:]' '[:upper:]' >> "$PRIMERS_REV_FILE"

TS=$(date +%Y%m%d_%H%M%S)
OUT_LOG="logs/run_${TS}.out"; ERR_LOG="logs/run_${TS}.err"
EMU_OUT_LOG="logs/emu_${TS}.out"; EMU_ERR_LOG="logs/emu_${TS}.err"

echo "============================================"
echo "libsQC:  ${RUN_LIBSQC}"
echo "Emu:     ${RUN_EMU}"
echo "Partition : $PARTITION"
echo "Time      : $TIME"
echo "CPUs      : $CPUS"
echo "Memory    : $MEM"
echo "Offset    : $OFFSET_FASTQS"
echo "Limit     : ${LIMIT_FASTQS}"
echo "Batch tag : ${BATCH_TAG}"
echo "Runs dir  : ${RUNS_DIR_DEFAULT}"
echo "Tables dir: ${TABLES_DIR_DEFAULT}"
echo "Plots dir : ${PLOTS_DIR_DEFAULT}"
echo "============================================"

export OMP_NUM_THREADS="$CPUS" MKL_NUM_THREADS="$CPUS" NUMEXPR_NUM_THREADS="$CPUS"

if [[ "$RUN_LIBSQC" -eq 1 ]]; then
  srun --partition="$PARTITION" --nodes=1 --ntasks=1 --cpus-per-task="$CPUS" \
       --mem="$MEM" --time="$TIME" --chdir="$WDIR" \
       --export=ALL,THREADS="$CPUS",PRIMER_FWD="$PRIMER_FWD",PRIMER_REV="$PRIMER_REV",SEQ_SUMMARY="$SEQ_SUMMARY",PRIMERS_FWD_FILE="$PRIMERS_FWD_FILE",PRIMERS_REV_FILE="$PRIMERS_REV_FILE" \
       /bin/bash workflow/run_libsQC.sh 1>"$OUT_LOG" 2>"$ERR_LOG"
fi

if [[ "$RUN_EMU" -eq 1 ]]; then
  export OMP_NUM_THREADS="$EMU_CPUS" MKL_NUM_THREADS="$EMU_CPUS" NUMEXPR_NUM_THREADS="$EMU_CPUS"
  SAVE_ASSIGN="${SAVE_ASSIGN:-0}"
  SKIP_ASSIGN="${SKIP_ASSIGN:-$((1 - SAVE_ASSIGN))}"

  srun --partition="$EMU_PARTITION" --nodes=1 --ntasks=1 --cpus-per-task="$EMU_CPUS" \
       --mem="$EMU_MEM" --time="$EMU_TIME" --chdir="$WDIR" \
       --export=ALL,THREADS="$EMU_CPUS",EMU_DB_ITS_DIR="$EMU_DB_ITS_DIR",EMU_DB_LSU_DIR="$EMU_DB_LSU_DIR",FASTQ_DIR_DEFAULT="${FASTQ_DIR_DEFAULT:-results/filtered}",LIMIT_FASTQS="${LIMIT_FASTQS}",OFFSET_FASTQS="${OFFSET_FASTQS}",BATCH_TAG="${BATCH_TAG}",RUNS_DIR="${RUNS_DIR_DEFAULT}",TABLES_DIR="${TABLES_DIR_DEFAULT}",PLOTS_DIR="${PLOTS_DIR_DEFAULT}",SAVE_JSON="${SAVE_JSON:-0}",SAVE_ASSIGN="${SAVE_ASSIGN}",SKIP_ASSIGN="${SKIP_ASSIGN}" \
       /bin/bash workflow/run_emu_amplicons.sh 1>"$EMU_OUT_LOG" 2>"$EMU_ERR_LOG"
fi

echo ">>> Done. Logs:"
[[ "$RUN_LIBSQC" -eq 1 ]] && echo "  $OUT_LOG / $ERR_LOG"
[[ "$RUN_EMU"  -eq 1 ]] && echo "  $EMU_OUT_LOG / $EMU_ERR_LOG"
