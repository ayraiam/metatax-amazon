#!/usr/bin/env bash
# ==========================================================
# Script: workflow/runall.sh
# Purpose: Submit libsQC run via Slurm (srun) and log output
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

# --- Help message ---
usage() {
  echo "Usage: bash workflow/runall.sh [options]"
  echo
  echo "Options:"
  echo "  --partition STR       Partition/queue name (default: short)"
  echo "  --time HH:MM:SS       Walltime (default: 04:00:00)"
  echo "  --cpus INT            CPUs per task (default: 4)"
  echo "  --mem STR             Memory (default: 16G)"
  echo "  --wd PATH             Working directory (default: current directory)"
  echo "  --primer-fwd SEQ      Forward primer sequence (optional)"
  echo "  --primer-rev SEQ      Reverse primer sequence (optional)"
  echo "  --seq-summary PATH    Path to sequencing_summary.txt (optional)"
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
    -h|--help) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

# --- Ensure logs directory exists ---
mkdir -p logs

# --- Create timestamped log files ---
TS=$(date +%Y%m%d_%H%M%S)
OUT_LOG="logs/run_${TS}.out"
ERR_LOG="logs/run_${TS}.err"

# --- Print run configuration ---
echo "============================================"
echo "Running libsQC with the following parameters:"
echo "  Partition : $PARTITION"
echo "  Time      : $TIME"
echo "  CPUs      : $CPUS"
echo "  Memory    : $MEM"
echo "  Work dir  : $WDIR"
[ -n "$PRIMER_FWD" ] && echo "  Primer FWD: $PRIMER_FWD"
[ -n "$PRIMER_REV" ] && echo "  Primer REV: $PRIMER_REV"
[ -n "$SEQ_SUMMARY" ] && echo "  Seq summary: $SEQ_SUMMARY"
echo "============================================"
echo "Logs will be written to:"
echo "  STDOUT -> $OUT_LOG"
echo "  STDERR -> $ERR_LOG"
echo

# --- Export thread variables for good measure ---
export OMP_NUM_THREADS="$CPUS"
export MKL_NUM_THREADS="$CPUS"
export NUMEXPR_NUM_THREADS="$CPUS"

# --- Submit the job ---
srun \
  --partition="$PARTITION" \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task="$CPUS" \
  --mem="$MEM" \
  --time="$TIME" \
  --chdir="$WDIR" \
  --export=ALL,THREADS="$CPUS",PRIMER_FWD="$PRIMER_FWD",PRIMER_REV="$PRIMER_REV",SEQ_SUMMARY="$SEQ_SUMMARY" \
  /bin/bash workflow/run_libsQC.sh \
  1>"$OUT_LOG" \
  2>"$ERR_LOG"

echo
echo ">>> Job finished. Check logs:"
echo "    $OUT_LOG"
echo "    $ERR_LOG"
