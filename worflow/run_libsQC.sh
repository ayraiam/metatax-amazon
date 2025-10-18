#!/usr/bin/env bash
# ==========================================================
# Script: run_libsQC.sh
# Purpose: Create, activate, and snapshot "libsQC" environment for QC + plotting
# Notes:
#   - Run with bash: bash run_libsQC.sh
#   - Requires conda/mamba (Miniforge recommended)
# ==========================================================

set -euo pipefail

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
      "seqkit>=2.6" fastqc=0.12.1 multiqc=1.21 -y
    status=$?
    set -e

    if [ $status -ne 0 ]; then
        echo "!!! Strict-priority solve failed. Retrying once with flexible channel priority..."
        mamba create -n libsQC \
          -c conda-forge -c bioconda --channel-priority flexible \
          python=3.11 "r-base>=4.3" "r-ggplot2>=3.4" "r-data.table" \
          "seqkit>=2.6" fastqc=0.12.1 multiqc=1.21 -y
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
plot_fastq_length_boxplots
