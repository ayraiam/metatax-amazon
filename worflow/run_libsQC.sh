#!/usr/bin/env bash
# ==========================================================
# Script: run_libsQC.sh
# Purpose: Create, activate, and snapshot "libsQC" environment for QC + plotting
# Notes:
#   - Run with bash: bash run_libsQC.sh
#   - Requires conda/mamba (Miniforge recommended)
# ==========================================================

set -euo pipefail

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

    echo ">>> Environment 'libsQC' not found. Creating it now (strict priority)..."
    set +e
    mamba create -n libsQC \
      -c conda-forge -c bioconda \
      python=3.11 "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" \
      "seqkit>=2.6" fastqc=0.12.1 multiqc=1.21 \
      nanostat nanoplot -y
    status=$?
    set -e

    if [ $status -ne 0 ]; then
        echo "!!! Strict-priority solve failed. Retrying once with flexible channel priority..."
        mamba create -n libsQC \
          -c conda-forge -c bioconda --channel-priority flexible \
          python=3.11 "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" \
          "seqkit>=2.6" fastqc=0.12.1 multiqc=1.21 \
          nanostat nanoplot -y
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
#   - Quick length & quality overview for ONT using NanoStat & NanoPlot
# ----------------------------------------------------------
quick_len_qual_overview() {
    echo ">>> Quick length & quality overview with NanoStat/NanoPlot ..."
    mkdir -p results/nanoplot results/nanostat

    # Per-file NanoStat TSVs
    for f in "${FASTQ_FILES[@]}"; do
        base="${f##*/}"; base="${base%.gz}"
        NanoStat --fastq "$f" --threads "${THREADS}" \
          --outdir "results/nanostat/${base}.stat" \
          --name "${base%.fastq}" >/dev/null
        echo "    NanoStat -> results/nanostat/${base}.stat"
    done

    # Combined NanoPlot: including --tsv_stats to also get a per-read table
    NanoPlot --threads "${THREADS}" --fastq "${FASTQ_FILES[@]}" \
      --N50 --loglength --plots hex dot kde --tsv_stats --raw \
      -o results/nanoplot/all
    echo ">>> NanoPlot HTML/PNGs + TSV in results/nanoplot/all"
}

# ----------------------------------------------------------
# Function: run_pycoqc_optional
# ----------------------------------------------------------
run_pycoqc_optional() {
    if [ -n "$SEQ_SUMMARY" ] && [ -f "$SEQ_SUMMARY" ]; then
        echo ">>> Running pycoQC on $SEQ_SUMMARY ..."
        mkdir -p results/pycoqc
        pycoQC -f "$SEQ_SUMMARY" -o results/pycoqc/pycoqc_report.html
        echo ">>> pycoQC report: results/pycoqc/pycoqc_report.html"
    else
        echo ">>> pycoQC skipped (set SEQ_SUMMARY to an existing sequencing_summary.txt)"
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
# Function: make_fastq_summary
#   - Build compact per-FASTQ summary with SeqKit (incl. N50)
#   - Ensure NanoPlot per-read TSV exists for downstream analysis
# ----------------------------------------------------------
make_fastq_summary() {
    echo ">>> Summarizing FASTQs with SeqKit and NanoPlot TSV ..."
    mkdir -p results/summary

    # Per-file stats (tab-separated) with SeqKit
    # Columns (seqkit -a -T): file, format, type, num_seqs, sum_len, min_len, avg_len, max_len, Q1, Q2, Q3, sum_gap, N50
    seqkit stats -a -T "${FASTQ_FILES[@]}" > results/summary/seqkit_stats.tsv
    echo ">>> Wrote results/summary/seqkit_stats.tsv"

    # Ensure NanoPlot per-read TSV was generated (quick_len_qual_overview should have produced this)
    if [ ! -f results/nanoplot/all/NanoPlot-data.tsv ]; then
        echo ">>> NanoPlot per-read TSV missing; generating now ..."
        NanoPlot --threads "${THREADS}" --fastq "${FASTQ_FILES[@]}" \
          --N50 --loglength --tsv_stats --raw \
          -o results/nanoplot/all >/dev/null
    fi
    echo ">>> NanoPlot per-read TSV: results/nanoplot/all/NanoPlot-data.tsv"
}

# ----------------------------------------------------------
# Function: qc_flags_from_nanoplot
#   - Derive simple flags per library from NanoPlot per-read TSV (no primers needed)
#   - Outputs: results/summary/qc_flags.tsv with:
#       file_base, total_reads, pct_len_200_700, pct_len_ge_1000, mean_read_q, N50
#   - N50 is joined from results/summary/seqkit_stats.tsv
# ----------------------------------------------------------
qc_flags_from_nanoplot() {
    echo ">>> Deriving QC flags from NanoPlot per-read TSV ..."
    local NP_TSV="results/nanoplot/all/NanoPlot-data.tsv"
    local STATS_TSV="results/summary/seqkit_stats.tsv"
    [ -f "$NP_TSV" ] || { echo "!!! Missing $NP_TSV (run quick_len_qual_overview first)"; return 1; }
    [ -f "$STATS_TSV" ] || { echo "!!! Missing $STATS_TSV (run make_fastq_summary first)"; return 1; }

    # Build a map: file_base -> N50 from seqkit_stats.tsv
    # seqkit_stats.tsv header contains "file" and "N50" columns.
    awk -F'\t' 'NR==1{
        for(i=1;i<=NF;i++){h[$i]=i}
        fn=h["file"]; n50=h["N50"];
        next
    }
    {
        # Extract basename only to match NanoPlot filename field
        f=$fn; gsub(/^.*\//,"",f);
        print f"\t"$n50
    }' "$STATS_TSV" > results/summary/.n50_map.tsv

    # Parse NanoPlot-data.tsv: detect column indices by header (robust across versions)
    # Expect columns containing (case-insensitive) names: "length", "mean_q", "filename" (or "source")
    awk -F'\t' '
    BEGIN{IGNORECASE=1}
    NR==1{
        for(i=1;i<=NF;i++){
            if($i ~ /length/) L=i
            else if($i ~ /mean_q|qmean|average_q/) Q=i
            else if($i ~ /filename|source|file/) F=i
        }
        if(!L || !Q || !F){
            print "!!! Could not find length/mean_q/filename columns in NanoPlot-data.tsv" > "/dev/stderr"; exit 1
        }
        next
    }
    {
        len=$L+0; q=$Q+0; file=$F
        gsub(/^.*\//,"",file)         # keep basename only
        key=file
        total[key]++
        sumq[key]+=q
        if(len>=200 && len<=700) its[key]++
        if(len>=1000) full16s[key]++
    }
    END{
        print "file_base\ttotal_reads\tpct_len_200_700\tpct_len_ge_1000\tmean_read_q\tN50"
        while((getline < "results/summary/.n50_map.tsv")>0){
            split($0,a,"\t"); n50[a[1]]=a[2]
        }
        for(k in total){
            pct_its = (its[k]>0? (its[k]/total[k])*100 : 0)
            pct_16s = (full16s[k]>0? (full16s[k]/total[k])*100 : 0)
            meanq   = (sumq[k]/total[k])
            printf "%s\t%d\t%.2f\t%.2f\t%.3f\t%s\n", k, total[k], pct_its, pct_16s, meanq, (k in n50 ? n50[k] : "NA")
        }
    }' "$NP_TSV" | sort -k1,1 > results/summary/qc_flags.tsv

    rm -f results/summary/.n50_map.tsv
    echo ">>> Wrote results/summary/qc_flags.tsv"
}

# ----------------------------------------------------------
# Function: plot_fastq_length_boxplots
#   - Builds results/lengths/all_lengths.tsv from FASTQ_FILES
#   - Then calls workflow/plot_fastq_lengths.R to make one figure with all boxplots
#   - Env var SAMPLE_N (optional): subsample N reads per file (0 = all reads)
# ----------------------------------------------------------
plot_fastq_length_boxplots() {
    echo ">>> Preparing read-length boxplots for FASTQs in data/ ..."
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

		# combine to a tidy table with header
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
        echo ">>> Plots saved to results/lengths/read_length_boxplots.png and .pdf"
}

# ==========================================================
# Calling the functions:
# ==========================================================
create_env_libsQC
check_versions
export_env
gather_fastq_files
build_fastq_meta
run_fastqc_all
run_multiqc
quick_len_qual_overview
run_pycoqc_optional
primer_spotcheck
make_fastq_summary
qc_flags_from_nanoplot
plot_fastq_length_boxplots
