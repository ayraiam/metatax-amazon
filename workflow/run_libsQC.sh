#!/usr/bin/env bash
# ==========================================================
# Script: run_libsQC.sh
# Purpose: Create, activate, and snapshot "libsQC" environment for QC + plotting
# Notes:
#   - Run with bash: bash run_libsQC.sh
#   - Requires conda/mamba (Miniforge recommended)
# ==========================================================

set -euo pipefail

# Initialize conda in non-interactive shells -----------------------
if ! command -v conda >/dev/null 2>&1; then  #
  for CAND in "$HOME/mambaforge" "$HOME/miniforge3" "$HOME/miniconda3" "/opt/conda"; do
    if [ -f "$CAND/etc/profile.d/conda.sh" ]; then
      # shellcheck disable=SC1091
      source "$CAND/etc/profile.d/conda.sh"
      break
    fi
  done
fi
# ---------------------------------------------------------------------------

#User-tunable vars for threads and optional primer/summary inputs
THREADS="${THREADS:-4}"
PRIMER_FWD="${PRIMER_FWD:-}"
PRIMER_REV="${PRIMER_REV:-}"
SEQ_SUMMARY="${SEQ_SUMMARY:-}"

# ----------------------------------------------------------
# Function: ensure_channels
#   Reset channels to a clean, known-good order and strict priority.
#   Also clear index cache to avoid stale solver metadata.
# ----------------------------------------------------------
ensure_channels() {
    echo ">>> Ensuring Conda channels (conda-forge, bioconda, defaults) with strict priority..."
    conda config --remove-key channels >/dev/null 2>&1 || true
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --add channels defaults
    conda config --set channel_priority strict
    echo ">>> Clearing index cache..."
    conda clean --index-cache -y >/dev/null 2>&1 || true

    echo ">>> Current channels:"
    conda config --show channels
}

# ----------------------------------------------------------
# Function: create_env_libsQC
#   1) If 'libsQC' exists, activate it.
#   2) Else, create with preferred channel order and pins.
#   3) If solver fails, retry once with flexible channel priority (fallback).
# ----------------------------------------------------------
create_env_libsQC() {
    echo ">>> Checking if conda/mamba environment 'libsQC' exists..."
    if conda env list | grep -qE '^libsQC\s'; then
        echo ">>> Environment 'libsQC' already exists. Activating it..."
        # shellcheck disable=SC1091
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate libsQC
        echo ">>> Environment 'libsQC' is now active."
        return 0
    fi

    # Only ensure channels & clean index if we need to CREATE the env
    ensure_channels

    # Choose solver with fallback to conda ------------------------
    local SOLVER="mamba"
    command -v mamba >/dev/null 2>&1 || SOLVER="conda"
    echo ">>> Using solver: ${SOLVER}"
    # -----------------------------------------------------------------------

    echo ">>> Environment 'libsQC' not found. Creating it now (strict priority)..."
    set +e
    ${SOLVER} create -n libsQC \
      -c conda-forge -c bioconda \
      python=3.11 "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" \
      "seqkit>=2.6" fastqc=0.12.1 multiqc=1.21 \
      nanostat nanoplot \
      "cutadapt>=4.5" nanofilt pycoqc -y
    status=$?
    set -e

    if [ $status -ne 0 ]; then
        echo "!!! Strict-priority solve failed. Retrying once with flexible channel priority..."
        ${SOLVER} create -n libsQC \
          -c conda-forge -c bioconda \
          python=3.11 "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" \
          "seqkit>=2.6" fastqc=0.12.1 multiqc=1.21 \
          nanostat nanoplot \
          "cutadapt>=4.5" nanofilt pycoqc -y
        echo ">>> Created 'libsQC' with flexible priority."
    else
        echo ">>> Environment 'libsQC' created successfully (strict priority)."
    fi

    # Activate the newly created env
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate libsQC
    echo ">>> Environment 'libsQC' is now active."
}

# ----------------------------------------------------------
# Function: check_versions
#   Print all packages/versions in the active environment.
# ----------------------------------------------------------
check_versions() {
    echo ">>> Listing all package versions in the current environment..."
    conda list
    echo ">>> Done."
}

# ----------------------------------------------------------
# Function: export_env
#   Save an exact snapshot of the active environment to envs/libsQC.yml
#   Creates envs/ only if it doesn't exist.
# ----------------------------------------------------------
export_env() {
    if [ ! -d envs ]; then
        echo ">>> Directory 'envs/' not found. Creating it..."
        mkdir -p envs
    fi
    echo ">>> Exporting environment to envs/libsQC.yml ..."
    conda env export --name libsQC > envs/libsQC.yml
    echo ">>> Environment exported: envs/libsQC.yml"
}

# ----------------------------------------------------------
# Function: gather_fastq_files
#   - Collects FASTQ paths once into a global array FASTQ_FILES
#   - Used by downstream functions to avoid recomputing the list
# ----------------------------------------------------------
gather_fastq_files() {
    echo ">>> Gathering FASTQ files from data/ ..."
    shopt -s nullglob
    FASTQ_FILES=( data/*.fastq.gz data/*.fq.gz data/*.fastq data/*.fq )
    shopt -u nullglob
    if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
        echo "!!! No FASTQ files found in data/ (expected *.fastq.gz, *.fq.gz, *.fastq, *.fq)"
        return 1
    fi
    echo ">>> Found ${#FASTQ_FILES[@]} FASTQ files."
}

# ----------------------------------------------------------
# Function: build_fastq_meta
#   - Scans global FASTQ_FILES
#   - Parses names with patterns: "_I_", "_II_", or "_III_"
#   - Outputs TSV: metadata/fastq_meta.tsv
# ----------------------------------------------------------
build_fastq_meta() {
    echo ">>> Building FASTQ manifest..."

    # Create metadata dir only if missing
    if [ ! -d metadata ]; then
        echo ">>> Directory 'metadata/' not found. Creating it..."
        mkdir -p metadata
    fi

    out="metadata/fastq_meta.tsv"
    echo -e "sample\treplicate\tmolecular_feature\tfile_name" > "$out"

    if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
        echo "!!! No FASTQ files available in FASTQ_FILES"
        return 1
    fi

    for f in "${FASTQ_FILES[@]}"; do
        bn="${f##*/}"                         # file_name with extension
        # strip known FASTQ extensions for parsing
        base="${bn%.fastq.gz}"
        base="${base%.fq.gz}"
        base="${base%.fastq}"
        base="${base%.fq}"

        # match: <sample>_(I|II|III)_<rest>
        if [[ "$base" =~ ^(.+)_((I|II|III))_(.+)$ ]]; then
            sample="${BASH_REMATCH[1]}"
            replicate="${BASH_REMATCH[2]}"    # I / II / III
            molecular_feature="${BASH_REMATCH[4]}"
            # normalize accidental leading/trailing underscores
            molecular_feature="${molecular_feature##_}"
            molecular_feature="${molecular_feature%%_}"
            printf "%s\t%s\t%s\t%s\n" "$sample" "$replicate" "$molecular_feature" "$bn" >> "$out"
        else
            echo "!!! Skipping (pattern not found): $bn" >&2
        fi
    done

    echo ">>> Manifest written: $out"
}

# ----------------------------------------------------------
# Function: run_fastqc_all
#   - Runs FastQC on all FASTQ files in data/
#   - Outputs to results/qc_raw
#   - THREADS env var controls parallelism (default 4)
# ----------------------------------------------------------
run_fastqc_all() {
    echo ">>> Running FastQC on data/ ..."
    mkdir -p results/qc_raw

    THREADS="${THREADS:-4}"

    if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
        echo "!!! No FASTQ files available in FASTQ_FILES"
        return 1
    fi

    fastqc -t "$THREADS" -o results/qc_raw "${FASTQ_FILES[@]}"
    echo ">>> FastQC done. Reports in results/qc_raw"
}

# ----------------------------------------------------------
# Function: run_multiqc
#   - Aggregates FastQC outputs (and any other supported tool logs) with MultiQC
#   - Scans results/qc_raw and writes to results/multiqc
# ----------------------------------------------------------
run_multiqc() {
    echo ">>> Running MultiQC over results/qc_raw ..."
    mkdir -p results/multiqc
    multiqc -o results/multiqc results/qc_raw
    echo ">>> MultiQC report: results/multiqc/multiqc_report.html"
}

# ----------------------------------------------------------
# Function: quick_len_qual_overview
#   - Combined NanoPlot HTML for all FASTQs
#   - Per-file NanoPlot raw TSVs so we can aggregate per library
# ----------------------------------------------------------
quick_len_qual_overview() {
    local LABEL="${1:-pre}"   # "pre" (default) or "post"
    echo ">>> Quick length & quality overview with NanoStat/NanoPlot (${LABEL}) ..."
    mkdir -p results/nanoplot results/nanostat results/nanoplot/per_file

    # Debug header to confirm inputs
    echo ">>> FASTQ_FILES count: ${#FASTQ_FILES[@]}"
    printf '    %s\n' "${FASTQ_FILES[@]:0:3}"
    [ ${#FASTQ_FILES[@]} -gt 3 ] && echo "    ..."

    # Per-file NanoStat TSVs (non-fatal)
    for f in "${FASTQ_FILES[@]}"; do
        base="$(basename "$f")"
        base="${base%.fastq.gz}"
        base="${base%.fq.gz}"
        base="${base%.fastq}"
        base="${base%.fq}"
        base="${base%.}"  # remove any stray trailing dot

        set +e
        NanoStat --fastq "$f" --threads "${THREADS}" \
          --outdir "results/nanostat/${base}.stat" \
          --name "${base}" >/dev/null
        ns_status=$?
        set -e
        if [ "$ns_status" -ne 0 ]; then
            echo "!!! NanoStat FAILED on $base (exit $ns_status) — continuing."
        else
            echo "    NanoStat -> results/nanostat/${base}.stat"
        fi
    done

    # Phase-specific combined NanoPlot (HTML + per-read TSV)
    local ALL_OUTDIR="results/nanoplot/all_${LABEL}"
    mkdir -p "${ALL_OUTDIR}"
    set +e
    NanoPlot --threads "${THREADS}" --fastq "${FASTQ_FILES[@]}" \
      --N50 --loglength --plots hex dot kde --tsv_stats --raw \
      -o "${ALL_OUTDIR}" \
      1>"${ALL_OUTDIR}/NanoPlot.stdout.log" \
      2>"${ALL_OUTDIR}/NanoPlot.stderr.log"
    np_status=$?
    set -e
    if [ "$np_status" -eq 0 ]; then
      echo ">>> Combined NanoPlot (${LABEL}): HTML/PNGs (+ per-read TSV) in ${ALL_OUTDIR}"
    else
      echo "!!! Combined NanoPlot (${LABEL}) FAILED (exit $np_status). Skipping combined overview."
      echo "    See ${ALL_OUTDIR}/NanoPlot.stderr.log for details."
    fi

    # Per-file NanoPlot raw TSVs (normalize base naming; non-fatal)
    echo ">>> Generating per-file NanoPlot raw TSVs (${LABEL}) ..."
    for f in "${FASTQ_FILES[@]}"; do
        b="$(basename "$f")"
        b="${b%.fastq.gz}"
        b="${b%.fq.gz}"
        b="${b%.fastq}"
        b="${b%.fq}"
        b="${b%.}"  # remove any stray trailing dot

        outdir="results/nanoplot/per_file/${b}"
        mkdir -p "$outdir"

        set +e
        NanoPlot --threads "${THREADS}" --fastq "$f" \
          --N50 --loglength --tsv_stats --raw \
          -o "$outdir" >/devnull 2>&1
        p_status=$?
        set -e
        if [ "$p_status" -ne 0 ]; then
            echo "!!! NanoPlot raw FAILED for $b (exit $p_status) — continuing."
        else
            echo "    NanoPlot raw -> $outdir/NanoPlot-data.tsv(.gz)"
        fi
    done
    echo ">>> Per-file NanoPlot raw TSVs ready under results/nanoplot/per_file/"
}

# ----------------------------------------------------------
# Function: run_pycoqc_optional
# ----------------------------------------------------------
run_pycoqc_optional() {
    if [ -n "$SEQ_SUMMARY" ] && [ -f "$SEQ_SUMMARY" ] && command -v pycoQC >/dev/null 2>&1; then
        echo ">>> Running pycoQC on $SEQ_SUMMARY ..."
        mkdir -p results/pycoqc
        pycoQC -f "$SEQ_SUMMARY" -o results/pycoqc/pycoqc_report.html
        echo ">>> pycoQC report: results/pycoqc/pycoqc_report.html"
    else
        echo ">>> pycoQC skipped (set SEQ_SUMMARY and ensure 'pycoqc' is installed)."
    fi
}

# ----------------------------------------------------------
# Function: primer_spotcheck
# ----------------------------------------------------------
primer_spotcheck() {
    if [ -z "$PRIMER_FWD" ] && [ -z "$PRIMER_REV" ]; then
        echo ">>> Primer spot-check skipped (set PRIMER_FWD/PRIMER_REV to enable)."
        return 0
    fi
    echo ">>> Primer/adapter residuals spot-check with cutadapt (no trimming)..."
    mkdir -p results/primer_checks
    for f in "${FASTQ_FILES[@]}"; do
        base="${f##*/}"; base="${base%.gz}"
        args=()
        [ -n "$PRIMER_FWD" ] && args+=( -g "^${PRIMER_FWD}" )
        [ -n "$PRIMER_REV" ] && args+=( -a "${PRIMER_REV}\$" )
        cutadapt "${args[@]}" --report=minimal --no-trim -j "${THREADS}" \
          -o /dev/null "$f" > "results/primer_checks/${base%.fastq}.cutadapt_report.txt"
        echo "    cutadapt -> results/primer_checks/${base%.fastq}.cutadapt_report.txt"
    done
    echo ">>> Review hit rates; high matches mean more trimming may be needed."
}

# ----------------------------------------------------------
# Function: filter_amplicons  (global mode)
# ----------------------------------------------------------
filter_amplicons() {
  echo ">>> Global length + Q filtering (multi-marker libraries: ITS + 16S + LSU) ..."
  mkdir -p results/filtered

  local MEANQ="${MEANQ:-10}"          # minimum mean quality
  local LEN_MIN="${LEN_MIN:-200}"     # min read length to keep
  local LEN_MAX="${LEN_MAX:-3300}"    # max read length to keep (extended for LSU)

  for f in "${FASTQ_FILES[@]}"; do
    base="${f##*/}"; base="${base%.gz}"; base="${base%.fastq}"; base="${base%.fq}"
    out="results/filtered/${base}.filtered.fastq.gz"

    echo "    Filtering $base  (keep ${LEN_MIN}-${LEN_MAX} bp, Q ≥ ${MEANQ})"

    # optional primer trim if sequences are known
    tmp="results/filtered/${base}.tmp.fastq.gz"
    if [ -n "$PRIMER_FWD" ] || [ -n "$PRIMER_REV" ]; then
      cutadapt -j "${THREADS}" \
        $( [ -n "$PRIMER_FWD" ] && echo -g "^${PRIMER_FWD}" ) \
        $( [ -n "$PRIMER_REV" ] && echo -a "${PRIMER_REV}\$" ) \
        -o "$tmp" "$f" >/dev/null
    else
      cp "$f" "$tmp"
    fi

    # main length + Q filter
    zcat "$tmp" \
      | NanoFilt -q "$MEANQ" -l "$LEN_MIN" --maxlength "$LEN_MAX" \
      | gzip > "$out"
    rm -f "$tmp"

    echo "    -> $out"
  done

  echo ">>> Filtered FASTQs written to results/filtered/"
}

# ----------------------------------------------------------
# Function: re_qc_filtered
#   - Re-run QC on results/filtered/*.fastq.gz
#   - Produces parallel reports so you can compare before/after
# ----------------------------------------------------------
re_qc_filtered() {
  echo ">>> Re-QC on filtered FASTQs ..."
  shopt -s nullglob
  FASTQ_FILES=( results/filtered/*.fastq.gz )
  shopt -u nullglob
  if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "!!! No filtered FASTQs found in results/filtered/"
    return 1
  fi

  # FastQC/MultiQC on filtered
  mkdir -p results/qc_filtered
  fastqc -t "${THREADS}" -o results/qc_filtered "${FASTQ_FILES[@]}"
  mkdir -p results/multiqc_filtered
  multiqc -o results/multiqc_filtered results/qc_filtered

  # NanoPlot/NanoStat + summaries on filtered (post)
  quick_len_qual_overview post
  make_fastq_summary
  qc_flags_from_nanoplot
  plot_fastq_length_boxplots post

  echo ">>> Re-QC complete. See results/multiqc_filtered and results/nanoplot/*"
}

# ----------------------------------------------------------
# Function: make_fastq_summary
#   - Build compact per-FASTQ summary with SeqKit (incl. N50)
#   - Ensure NanoPlot per-read TSV exists for downstream analysis
#   - Relaxed: do NOT force-create combined 'all' anymore; accept pre/post dirs
# ----------------------------------------------------------
make_fastq_summary() {
    echo ">>> Summarizing FASTQs with SeqKit and NanoPlot TSV ..."
    mkdir -p results/summary

    seqkit stats -a -T "${FASTQ_FILES[@]}" > results/summary/seqkit_stats.tsv
    echo ">>> Wrote results/summary/seqkit_stats.tsv"

    # If at least one combined NanoPlot exists (pre or post), we're good; else just skip creating.
    if [ -f results/nanoplot/all_pre/NanoPlot-data.tsv ] || [ -f results/nanoplot/all_post/NanoPlot-data.tsv ]; then
        echo ">>> Found combined NanoPlot per-read TSV (pre/post)."
    else
        echo ">>> No combined NanoPlot per-read TSV found (pre/post). Skipping generation."
    fi
}

# ----------------------------------------------------------
# Function: qc_flags_from_nanoplot
#   - Aggregate per-file NanoPlot raw TSVs (quals/lengths) into summary per library
#   - Columns: file_base, total_reads, pct_len_200_700, pct_len_ge_1000, mean_read_q
#   - Output: results/summary/qc_flags.tsv
# ----------------------------------------------------------
qc_flags_from_nanoplot() {
    echo ">>> Deriving QC flags from per-file NanoPlot raw TSVs ..."
    local STATS_TSV="results/summary/seqkit_stats.tsv"
    [ -f "$STATS_TSV" ] || { echo "!!! Missing $STATS_TSV (run make_fastq_summary first)"; return 1; }

    echo ">>> Aggregating per-file NanoPlot raw TSVs ..."
    mkdir -p results/summary

    # Write header first
    echo -e "file_base\ttotal_reads\tpct_len_200_700\tpct_len_ge_1000\tmean_read_q" > results/summary/qc_flags.tsv

    shopt -s nullglob
    local found=false
    for tsv in results/nanoplot/per_file/*/NanoPlot-data.tsv results/nanoplot/per_file/*/NanoPlot-data.tsv.gz; do
        [ -e "$tsv" ] || continue
        found=true
        base="$(basename "$(dirname "$tsv")")"

        if [[ "$tsv" == *.gz ]]; then
            reader="gunzip -c \"$tsv\""
        else
            reader="cat \"$tsv\""
        fi

        eval $reader | awk -F'\t' -v base="$base" '
          BEGIN{IGNORECASE=1}
          NR==1{
            for(i=1;i<=NF;i++){
              if($i~/^lengths?$|length/) L=i
              else if($i~/^quals$|mean_q|qmean|average_q|quality/) Q=i
            }
            next
          }
          {
            len=$L+0; q=$Q+0
            total++; sumq+=q
            if(len>=200 && len<=700) its++
            if(len>=1000) full16s++
          }
          END{
            if(total>0){
              pct_its  = (its/total)*100
              pct_16s  = (full16s/total)*100
              meanq    = (sumq/total)
              printf "%s\t%d\t%.2f\t%.2f\t%.3f\n", base, total, pct_its, pct_16s, meanq
            }
          }' >> results/summary/qc_flags.tsv
    done
    shopt -u nullglob

    if [ "$found" = false ]; then
        echo "!!! No per-file NanoPlot raw TSVs found. Did quick_len_qual_overview run?"
        rm -f results/summary/qc_flags.tsv
        return 1
    fi

    echo ">>> Wrote simplified QC summary: results/summary/qc_flags.tsv"
}

# ----------------------------------------------------------
# Function: plot_fastq_length_boxplots
#   - Builds results/lengths/all_lengths.tsv for R plotting (consumed by the R script)
#   - Then renames outputs to *_<LABEL>.*
#   - Also writes a copy of the input table to all_lengths_<LABEL>.tsv for bookkeeping
# ----------------------------------------------------------
plot_fastq_length_boxplots() {
    local LABEL="${1:-pre}"   # "pre" (default) or "post"
    echo ">>> Preparing read-length boxplots for FASTQs in data/ (${LABEL}) ..."
    [ -d results/lengths ] || mkdir -p results/lengths

    if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
        echo "!!! No FASTQ files available in FASTQ_FILES"
        return 1
    fi

    SAMPLE_N="${SAMPLE_N:-0}"

    echo ">>> Extracting read lengths (this may take a bit for large files)..."
    for f in "${FASTQ_FILES[@]}"; do
        base="${f##*/}"
        base="${base%.gz}"
        base="${base%.fastq}"
        base="${base%.fq}"

        if [ "$SAMPLE_N" -gt 0 ]; then
          seqkit sample -n "$SAMPLE_N" "$f" \
            | seqkit fx2tab -n -l - \
            | awk -F'\t' -v OFS='\t' -v s="$base" '{print s,$NF}' \
            > "results/lengths/${base}.len.tsv"
        else
          seqkit fx2tab -n -l "$f" \
            | awk -F'\t' -v OFS='\t' -v s="$base" '{print s,$NF}' \
            > "results/lengths/${base}.len.tsv"
        fi
        echo "    wrote results/lengths/${base}.len.tsv"
    done

    {
      echo -e "sample\tlength"
      cat results/lengths/*.len.tsv
    } > results/lengths/all_lengths.tsv

    echo ">>> Combined table: results/lengths/all_lengths.tsv"
    rm -f results/lengths/*.len.tsv
    echo ">>> Removed temporary per-file .len.tsv files; only all_lengths.tsv retained."

    if [ ! -f workflow/plot_fastq_lengths.R ]; then
        echo "!!! Missing workflow/plot_fastq_lengths.R — please create it."
        return 1
    fi

    echo ">>> Running R to generate the boxplot figure..."
    Rscript workflow/plot_fastq_lengths.R

    # Keep phase-specific copies
    cp results/lengths/all_lengths.tsv "results/lengths/all_lengths_${LABEL}.tsv" || true
    [ -f results/lengths/read_length_boxplots.png ] && mv results/lengths/read_length_boxplots.png "results/lengths/read_length_boxplots_${LABEL}.png"
    [ -f results/lengths/read_length_boxplots.pdf ] && mv results/lengths/read_length_boxplots.pdf "results/lengths/read_length_boxplots_${LABEL}.pdf"

    echo ">>> Plots saved to results/lengths/read_length_boxplots_${LABEL}.png and .pdf"
}

# ----------------------------------------------------------
# Function: log_run_report
#   - Records runtime, threads, system info, and per-function timing
#   - Saves a summary text file to logs/run_report_<timestamp>.txt
# ----------------------------------------------------------
log_run_report() {
    # Ensure logs directory exists
    if [ ! -d logs ]; then
        echo ">>> Creating logs/ directory..."
        mkdir -p logs
    fi

    local logfile="logs/run_report_$(date +%Y%m%d_%H%M%S).txt"
    local end_time=$(date +%s)
    local runtime=$((end_time - START_TIME))
    local minutes=$((runtime / 60))
    local seconds=$((runtime % 60))

    # System and environment info
    local host=$(hostname)
    local env_hash=$(conda env export --name libsQC 2>/dev/null | sha256sum | cut -c1-12)

    echo "=================================================="    > "$logfile"
    echo " QC RUN REPORT — $(date)"                              >> "$logfile"
    echo "=================================================="   >> "$logfile"
    echo "Environment : libsQC"                                  >> "$logfile"
    echo "System       : ${host}"                                >> "$logfile"
    echo "Conda env fingerprint: ${env_hash}"                    >> "$logfile"
    echo "Threads used: ${THREADS}"                              >> "$logfile"
    echo                                                        >> "$logfile"

    # Aggregated timing summary (Total, Count, Avg)
    echo "Function timing summary (aggregated):"                 >> "$logfile"
    printf "%-30s %10s %8s %10s\n" "Function" "Total(s)" "n" "Avg(s)" >> "$logfile"
    printf "%-30s %10s %8s %10s\n" "--------" "--------" "--" "------" >> "$logfile"
    awk -F'\t' '
      {tot[$1]+=$2; cnt[$1]++}
      END{
        for (k in tot) {
          printf "%-30s %10d %8d %10.1f\n", k, tot[k], cnt[k], (cnt[k]? tot[k]/cnt[k] : 0)
        }
      }' logs/.timing.tsv >> "$logfile" 2>/dev/null || true

    echo                                                        >> "$logfile"
    echo "Report saved to: $logfile"                             >> "$logfile"
    echo ">>> Run report complete."
    echo ">>> Total runtime: ${minutes}m ${seconds}s — report saved to ${logfile}"
}

# ----------------------------------------------------------
# Helper: time_function
#   - Wraps a function call and logs its duration
#   - Usage: time_function run_fastqc_all
# ----------------------------------------------------------
time_function() {
    local fn="$1"
    local start=$(date +%s)
    echo ">>> Running $fn ..."
    local status=0

    set +e
    $fn
    status=$?
    set -e

    local end=$(date +%s)
    local dur=$((end - start))
    mkdir -p logs
    echo -e "${fn}\t${dur}" >> logs/.timing.tsv
    echo ">>> $fn completed in ${dur}s"
    return $status
}

# ==========================================================
# Calling the functions (QC pipeline main flow)
# ==========================================================

START_TIME=$(date +%s)

# Reset per-run timing file so reports don't accumulate across runs
[ -d logs ] || mkdir -p logs
: > logs/.timing.tsv

#Environment setup --------------------------------------------------------
time_function create_env_libsQC
time_function check_versions
time_function export_env

#Input discovery ---------------------------------------------------------
time_function gather_fastq_files
time_function build_fastq_meta

#Initial QC on raw data --------------------------------------------------
time_function run_fastqc_all
time_function run_multiqc

# Nanopore-specific QC (pre-filter) -----------------------------------
time_function 'quick_len_qual_overview pre'

# Summaries and QC flags --------------------------------------------------
time_function make_fastq_summary
time_function qc_flags_from_nanoplot
time_function 'plot_fastq_length_boxplots pre'

# Filtering step ----------------------------------------------------------
time_function filter_amplicons
time_function re_qc_filtered

#Final report ------------------------------------------------------------
log_run_report
