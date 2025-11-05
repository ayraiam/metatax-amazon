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
if ! command -v conda >/dev/null 2>&1; then
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

# primer trimming tolerance, list support, and report dirs
PRIMER_ERR="${PRIMER_ERR:-0.10}"        # max mismatch rate for cutadapt matches
PRIMER_FWD_LIST="${PRIMER_FWD_LIST:-$PRIMER_FWD}"   # CSV or single value
PRIMER_REV_LIST="${PRIMER_REV_LIST:-$PRIMER_REV}"   # CSV or single value

RESULTS="${RESULTS:-results}"

PRIMER_CHECK_DIR="${RESULTS}/primer_checks"
PRIMER_TRIM_DIR="${RESULTS}/primer_trimming"

# If runall.sh exported primer list files, read them, uppercase, and
# convert to CSV so parse_primers() below can consume them.
if [[ -n "${PRIMERS_FWD_FILE:-}" && -s "${PRIMERS_FWD_FILE}" ]]; then
  PRIMER_FWD_LIST="$(awk 'NF{gsub(/[ \t\r]/,""); print toupper($0)}' "${PRIMERS_FWD_FILE}" | paste -sd, -)"
fi
if [[ -n "${PRIMERS_REV_FILE:-}" && -s "${PRIMERS_REV_FILE}" ]]; then
  PRIMER_REV_LIST="$(awk 'NF{gsub(/[ \t\r]/,""); print toupper($0)}' "${PRIMERS_REV_FILE}" | paste -sd, -)"
fi

# internal arrays built from lists (parsed later)
FORWARD_PRIMERS=()
REVERSE_PRIMERS=()

# parse CSV env vars into arrays (trim spaces, drop empties)
parse_primers() {
  IFS=',' read -ra _tmpF <<< "$PRIMER_FWD_LIST"
  IFS=',' read -ra _tmpR <<< "$PRIMER_REV_LIST"
  FORWARD_PRIMERS=()
  REVERSE_PRIMERS=()
  for p in "${_tmpF[@]:-}"; do
    p="${p//[[:space:]]/}"; [[ -n "$p" ]] && FORWARD_PRIMERS+=("$p")
  done
  for p in "${_tmpR[@]:-}"; do
    p="${p//[[:space:]]/}"; [[ -n "$p" ]] && REVERSE_PRIMERS+=("$p")
  done
}

# emit cutadapt -g/-a flags for all primers (anchored)
cutadapt_primer_flags() {
  local flags=()
  for p in "${FORWARD_PRIMERS[@]}";  do flags+=( -g "^${p}" ); done
  for p in "${REVERSE_PRIMERS[@]}";  do flags+=( -a "${p}\$" ); done
  printf '%q ' "${flags[@]}"
}

# --- Convenience: per-group primer pairs (FOR CLASSIFICATION) ------------
# canonical group primer definitions for demultiplexing
PRIM_FWD_ARCHAEA="TTCCGGTTGATCCTGCCGGA"
PRIM_REV_ARCHAEA="TACGGWTACCTTGTTACGACTT"

PRIM_FWD_ASCOMIC="TCCGTAGGTGAACCTGCGG"
PRIM_REV_ASCOMIC="TCCTCCGCTTATTGATATGC"

PRIM_FWD_BAC="AGAGTTTGATCMTGGCTCAG"
PRIM_REV_BAC="GGTTACCTTGTTACGACTT"

PRIM_FWD_BASID="TACTACCACCAAGATCT"
PRIM_REV_BASID="ACCCGCTGAACTTAAGC"

# ----------------------------------------------------------
# Function: ensure_channels
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
# ----------------------------------------------------------
create_env_libsQC() {
    echo ">>> Checking if conda/mamba environment 'libsQC' exists..."
    if conda env list | grep -qE '^libsQC\s'; then
        echo ">>> Environment 'libsQC' already exists. Activating it..."
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate libsQC
        echo ">>> Environment 'libsQC' is now active."
        return 0
    fi

    ensure_channels

    local SOLVER="mamba"
    command -v mamba >/dev/null 2>&1 || SOLVER="conda"
    echo ">>> Using solver: ${SOLVER}"

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

    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate libsQC
    echo ">>> Environment 'libsQC' is now active."
}

check_versions() { echo ">>> Listing all package versions in the current environment..."; conda list; echo ">>> Done."; }

export_env() {
    [ -d envs ] || mkdir -p envs
    echo ">>> Exporting environment to envs/libsQC.yml ..."
    conda env export --name libsQC > envs/libsQC.yml
    echo ">>> Environment exported: envs/libsQC.yml"
}

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

### limit to the first 3 fastqs as requested
limit_to_three_fastqs() {
    if [ ${#FASTQ_FILES[@]} -gt 3 ]; then
        echo ">>> Limiting to first 3 FASTQs as requested."
        FASTQ_FILES=( "${FASTQ_FILES[@]:0:3}" )
    fi
    printf ">>> Using FASTQs:\n"; printf "    %s\n" "${FASTQ_FILES[@]}"
}

build_fastq_meta() {
    echo ">>> Building FASTQ manifest..."
    [ -d metadata ] || mkdir -p metadata
    out="metadata/fastq_meta.tsv"
    echo -e "sample\treplicate\tmolecular_feature\tfile_name" > "$out"

    for f in "${FASTQ_FILES[@]}"; do
        bn="${f##*/}"
        base="${bn%.fastq.gz}"; base="${base%.fq.gz}"; base="${base%.fastq}"; base="${base%.fq}"
        if [[ "$base" =~ ^(.+)_((I|II|III))_(.+)$ ]]; then
            sample="${BASH_REMATCH[1]}"; replicate="${BASH_REMATCH[2]}"; molecular_feature="${BASH_REMATCH[4]}"
            molecular_feature="${molecular_feature##_}"; molecular_feature="${molecular_feature%%_}"
            printf "%s\t%s\t%s\t%s\n" "$sample" "$replicate" "$molecular_feature" "$bn" >> "$out"
        else
            echo "!!! Skipping (pattern not found): $bn" >&2
        fi
    done
    echo ">>> Manifest written: $out"
}

# ----------------------------------------------------------
# Function: demux_by_primers (classify into 5 groups)
#   Strategy: iterative partition with cutadapt
#   - For each input FASTQ:
#       remaining := file
#       for group in Archaea, Ascomic, Bac, Basid:
#           grab reads with -g F -a R (unanchored; tolerant with -e)
#           write to results/groups/<GROUP>/data/<file>.fastq.gz
#           write non-matching to a new remaining; continue
#       leftover -> Unknown
#   (We do not trim here; only classify.)
# ----------------------------------------------------------
### function for grouping
demux_by_primers() {
  echo ">>> Demultiplexing reads into groups by primer pairs ..."
  local groups_root="${RESULTS}/groups"
  mkdir -p "${groups_root}"

  for inF in "${FASTQ_FILES[@]}"; do
    local bn="${inF##*/}"; local b="${bn%.gz}"; b="${b%.fastq}"; b="${b%.fq}"
    local rem="${groups_root}/_tmp_${b}.remaining.fastq.gz"
    # ensure gz input
    if [[ "$inF" == *.gz ]]; then
      cp "$inF" "$rem"
    else
      gzip -c "$inF" > "$rem"
    fi

    # helper to run one group extraction
    _grab_group() {
      local grp="$1"; local FWD="$2"; local REV="$3"
      local outdir="${groups_root}/${grp}/data"
      mkdir -p "$outdir"
      local outF="${outdir}/${b}.${grp}.fastq.gz"
      local next="${groups_root}/_tmp_${b}.${grp}_remaining.fastq.gz"
      # capture matches to OUT; non-matches to NEXT
      cutadapt -j "${THREADS}" -e "${PRIMER_ERR}" \
        -g "${FWD}" -a "${REV}" \
        -o "$outF" \
        --untrimmed-output "$next" \
        "$rem" \
        > "${outdir}/${b}.${grp}.demux_report.txt"
      mv -f "$next" "$rem"
    }

    _grab_group "Archaea" "${PRIM_FWD_ARCHAEA}" "${PRIM_REV_ARCHAEA}"
    _grab_group "Ascomic" "${PRIM_FWD_ASCOMIC}" "${PRIM_REV_ASCOMIC}"
    _grab_group "Bac"     "${PRIM_FWD_BAC}"     "${PRIM_REV_BAC}"
    _grab_group "Basid"   "${PRIM_FWD_BASID}"   "${PRIM_REV_BASID}"

    # leftover = Unknown
    local outdirU="${groups_root}/Unknown/data"
    mkdir -p "$outdirU"
    mv -f "$rem" "${outdirU}/${b}.Unknown.fastq.gz"
  done

  echo ">>> Demultiplexing complete. See ${groups_root}/*/data/"
}

run_fastqc_all() {
    echo ">>> Running FastQC ..."
    mkdir -p "${RESULTS}/qc_raw"
    if [ ${#FASTQ_FILES[@]} -eq 0 ]; then echo "!!! No FASTQs"; return 1; fi
    fastqc -t "$THREADS" -o "${RESULTS}/qc_raw" "${FASTQ_FILES[@]}"
    echo ">>> FastQC done -> ${RESULTS}/qc_raw"
}

run_multiqc() {
    echo ">>> Running MultiQC ..."
    mkdir -p "${RESULTS}/multiqc"
    multiqc -o "${RESULTS}/multiqc" "${RESULTS}/qc_raw"
    echo ">>> MultiQC -> ${RESULTS}/multiqc/multiqc_report.html"
}

quick_len_qual_overview() {
    local LABEL="${1:-pre}"
    echo ">>> NanoStat/NanoPlot (${LABEL}) ..."
    mkdir -p "${RESULTS}/nanoplot" "${RESULTS}/nanostat" "${RESULTS}/nanoplot/per_file"

    # per-file NanoStat
    for f in "${FASTQ_FILES[@]}"; do
        base="$(basename "$f")"; base="${base%.fastq.gz}"; base="${base%.fq.gz}"; base="${base%.fastq}"; base="${base%.fq}"; base="${base%.}"
        set +e
        NanoStat --fastq "$f" --threads "${THREADS}" \
          --outdir "${RESULTS}/nanostat/${base}.stat" \
          --name "${base}" >/dev/null
        set -e
    done

    local ALL_OUTDIR="${RESULTS}/nanoplot/all_${LABEL}"
    mkdir -p "${ALL_OUTDIR}"
    set +e
    NanoPlot --threads "${THREADS}" --fastq "${FASTQ_FILES[@]}" \
      --N50 --loglength --plots hex dot kde --tsv_stats --raw \
      -o "${ALL_OUTDIR}" \
      1>"${ALL_OUTDIR}/NanoPlot.stdout.log" \
      2>"${ALL_OUTDIR}/NanoPlot.stderr.log"
    set -e

    # Per-file NanoPlot raw TSV
    for f in "${FASTQ_FILES[@]}"; do
        b="$(basename "$f")"; b="${b%.fastq.gz}"; b="${b%.fq.gz}"; b="${b%.fastq}"; b="${b%.fq}"; b="${b%.}"
        outdir="${RESULTS}/nanoplot/per_file/${b}"
        mkdir -p "$outdir"
        set +e
        NanoPlot --threads "${THREADS}" --fastq "$f" \
          --N50 --loglength --tsv_stats --raw \
          -o "$outdir" \
          1>"$outdir/NanoPlot.stdout.log" \
          2>"$outdir/NanoPlot.stderr.log"
        set -e
    done
    echo ">>> Per-file TSVs -> ${RESULTS}/nanoplot/per_file/"
}

run_pycoqc_optional() {
    if [ -n "$SEQ_SUMMARY" ] && [ -f "$SEQ_SUMMARY" ] && command -v pycoQC >/dev/null 2>&1; then
        echo ">>> Running pycoQC on $SEQ_SUMMARY ..."
        mkdir -p "${RESULTS}/pycoqc"
        pycoQC -f "$SEQ_SUMMARY" -o "${RESULTS}/pycoqc/pycoqc_report.html"
        echo ">>> pycoQC report: ${RESULTS}/pycoqc/pycoqc_report.html"
    else
        echo ">>> pycoQC skipped (set SEQ_SUMMARY and ensure 'pycoqc' is installed)."
    fi
}

primer_spotcheck() {
    parse_primers
    if [ ${#FORWARD_PRIMERS[@]} -eq 0 ] && [ ${#REVERSE_PRIMERS[@]} -eq 0 ]; then
        echo ">>> Primer spot-check skipped (no primers provided)."
        return 0
    fi
    echo ">>> Primer/adapter residuals spot-check (no trimming) ..."
    mkdir -p "${PRIMER_CHECK_DIR}"
    for f in "${FASTQ_FILES[@]}"; do
        base="${f##*/}"; base="${base%.gz}"
        flags=$(cutadapt_primer_flags)
        eval cutadapt $flags -e "$PRIMER_ERR" --no-trim -j "${THREADS}" \
          -o /dev/null "$f" > "${PRIMER_CHECK_DIR}/${base%.fastq}.cutadapt_pre_check.txt"
    done
    echo ">>> Pre-check reports -> ${PRIMER_CHECK_DIR}"
}

filter_amplicons() {
  echo ">>> Global length + Q filtering ..."
  mkdir -p "${RESULTS}/filtered" "${PRIMER_TRIM_DIR}"

  local MEANQ="${MEANQ:-10}"
  local LEN_MIN="${LEN_MIN:-200}"
  local LEN_MAX="${LEN_MAX:-3300}"

  parse_primers

  for f in "${FASTQ_FILES[@]}"; do
    base="${f##*/}"; base="${base%.gz}"; base="${base%.fastq}"; base="${base%.fq}"
    out="${RESULTS}/filtered/${base}.filtered.fastq.gz"
    tmp="${RESULTS}/filtered/${base}.tmp.fastq.gz"

    if [ ${#FORWARD_PRIMERS[@]} -gt 0 ] || [ ${#REVERSE_PRIMERS[@]} -gt 0 ]; then
      flags=$(cutadapt_primer_flags)
      eval cutadapt $flags -j "${THREADS}" -e "$PRIMER_ERR" -o "$tmp" "$f" \
        > "${PRIMER_TRIM_DIR}/${base}.cutadapt_trim_report.txt"
    else
      if [[ "$f" == *.gz ]]; then cp "$f" "$tmp"; else gzip -c "$f" > "$tmp"; fi
    fi

    zcat "$tmp" | NanoFilt -q "$MEANQ" -l "$LEN_MIN" --maxlength "$LEN_MAX" | gzip > "$out"
    rm -f "$tmp"
  done

  echo ">>> Filtered -> ${RESULTS}/filtered/"
}

re_qc_filtered() {
  echo ">>> Re-QC on filtered FASTQs ..."
  shopt -s nullglob
  FASTQ_FILES=( "${RESULTS}/filtered/"*.fastq.gz )
  shopt -u nullglob
  [ ${#FASTQ_FILES[@]} -eq 0 ] && { echo "!!! No filtered FASTQs"; return 1; }

  mkdir -p "${RESULTS}/qc_filtered"
  fastqc -t "${THREADS}" -o "${RESULTS}/qc_filtered" "${FASTQ_FILES[@]}"
  mkdir -p "${RESULTS}/multiqc_filtered"
  multiqc -o "${RESULTS}/multiqc_filtered" "${RESULTS}/qc_filtered"

  quick_len_qual_overview post
  make_fastq_summary
  qc_flags_from_nanoplot
  plot_fastq_length_boxplots post "${RESULTS}/lengths"
}

make_fastq_summary() {
    echo ">>> Summarizing FASTQs with SeqKit ..."
    mkdir -p "${RESULTS}/summary"
    seqkit stats -a -T "${FASTQ_FILES[@]}" > "${RESULTS}/summary/seqkit_stats.tsv"
    if [ -f "${RESULTS}/nanoplot/all_pre/NanoPlot-data.tsv" ] || [ -f "${RESULTS}/nanoplot/all_post/NanoPlot-data.tsv" ] \
      || [ -f "${RESULTS}/nanoplot/all_pre/NanoPlot-data.txt" ] || [ -f "${RESULTS}/nanoplot/all_post/NanoPlot-data.txt" ]; then
        echo ">>> Found combined NanoPlot per-read TSV (pre/post)."
    else
        echo ">>> No combined NanoPlot per-read TSV found (pre/post). Skipping generation."
    fi
}

qc_flags_from_nanoplot() {
    echo ">>> Deriving QC flags from per-file NanoPlot raw TSVs ..."
    local STATS_TSV="${RESULTS}/summary/seqkit_stats.tsv"
    [ -f "$STATS_TSV" ] || { echo "!!! Missing $STATS_TSV"; return 1; }

    mkdir -p "${RESULTS}/summary"
    echo -e "file_base\ttotal_reads\tpct_len_200_700\tpct_len_ge_1000\tmean_read_q" > "${RESULTS}/summary/qc_flags.tsv"

    shopt -s nullglob
    local found=false
    for tsv in \
      "${RESULTS}/nanoplot/per_file/"*/NanoPlot-data.tsv \
      "${RESULTS}/nanoplot/per_file/"*/NanoPlot-data.tsv.gz \
      "${RESULTS}/nanoplot/per_file/"*/NanoPlot-data.txt \
      "${RESULTS}/nanoplot/per_file/"*/NanoPlot-data.txt.gz; do

        [ -e "$tsv" ] || continue
        found=true
        base="$(basename "$(dirname "$tsv")")"

        if [[ "$tsv" == *.gz ]]; then reader="gunzip -c \"$tsv\""; else reader="cat \"$tsv\""; fi

        eval $reader | awk -F'\t' -v base="$base" '
          BEGIN{IGNORECASE=1}
          NR==1{for(i=1;i<=NF;i++){if($i~/^lengths?$|length/)L=i;else if($i~/^quals$|mean_q|qmean|average_q|quality/)Q=i} next}
          {len=$L+0; q=$Q+0; total++; sumq+=q; if(len>=200 && len<=700) its++; if(len>=1000) full16s++}
          END{if(total>0){printf "%s\t%d\t%.2f\t%.2f\t%.3f\n", base, total, (its/total)*100, (full16s/total)*100, (sumq/total)}}
        ' >> "${RESULTS}/summary/qc_flags.tsv"
    done
    shopt -u nullglob

    if [ "$found" = false ]; then
        echo "!!! No per-file NanoPlot raw TSVs found."
        rm -f "${RESULTS}/summary/qc_flags.tsv"
        return 0
    fi
    echo ">>> Wrote ${RESULTS}/summary/qc_flags.tsv"
}

plot_fastq_length_boxplots() {
    local LABEL="${1:-pre}"
    local OUT_DIR="${2:-results/lengths}"
    echo ">>> Preparing read-length boxplots for FASTQs (${LABEL}) -> ${OUT_DIR} ..."
    [ -d "${OUT_DIR}" ] || mkdir -p "${OUT_DIR}"

    if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
        echo "!!! No FASTQ files available in FASTQ_FILES — skipping ${LABEL} plot for OUT_DIR=${OUT_DIR}"
        return 0
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
            > "${OUT_DIR}/${base}.len.tsv"
        else
          seqkit fx2tab -n -l "$f" \
            | awk -F'\t' -v OFS='\t' -v s="$base" '{print s,$NF}' \
            > "${OUT_DIR}/${base}.len.tsv"
        fi
        echo "    wrote ${OUT_DIR}/${base}.len.tsv"
    done

    {
      echo -e "sample\tlength"
      cat "${OUT_DIR}"/*.len.tsv
    } > "${OUT_DIR}/all_lengths.tsv"

    echo ">>> Combined table: ${OUT_DIR}/all_lengths.tsv"
    rm -f "${OUT_DIR}"/*.len.tsv
    echo ">>> Removed temporary per-file .len.tsv files; only all_lengths.tsv retained."

    if [ ! -f workflow/plot_fastq_lengths.R ]; then
        echo "!!! Missing workflow/plot_fastq_lengths.R — skipping plot."
        return 0
    fi

    echo ">>> Running R to generate the boxplot figure for ${LABEL} in ${OUT_DIR} ..."
    set +e
    Rscript workflow/plot_fastq_lengths.R "${OUT_DIR}/all_lengths.tsv" "${OUT_DIR}" "${LABEL}"
    rstat=$?
    set -e
    if [ $rstat -ne 0 ]; then
      echo "!!! R failed for OUT_DIR=${OUT_DIR}, LABEL=${LABEL} (exit ${rstat}). Listing directory for clues:"
      ls -lh "${OUT_DIR}" || true
      return 0
    fi

    # suffix copy (for downstream aggregation if you want it)
    if [ -f "${OUT_DIR}/all_lengths.tsv" ]; then
      cp "${OUT_DIR}/all_lengths.tsv" "${OUT_DIR}/all_lengths_${LABEL}.tsv"
      echo ">>> Wrote ${OUT_DIR}/all_lengths_${LABEL}.tsv"
    fi

    # verify files are really there, then print exact paths
    if [ -f "${OUT_DIR}/read_length_boxplots_${LABEL}.png" ]; then
      echo ">>> Plots saved:"
      echo "    ${OUT_DIR}/read_length_boxplots_${LABEL}.png"
      echo "    ${OUT_DIR}/read_length_boxplots_${LABEL}.pdf"
    else
      echo "!!! Plot files not found in ${OUT_DIR} after R ran. Contents:"
      ls -lh "${OUT_DIR}" || true
    fi
}

# verify primers after filtering and summarize to TSV
verify_primer_removal() {
  parse_primers
  if [ ${#FORWARD_PRIMERS[@]} -eq 0 ] && [ ${#REVERSE_PRIMERS[@]} -eq 0 ]; then
    echo ">>> verify_primer_removal: skipped (no primers for this group)"
    return 0
  fi
  echo ">>> Verifying primer removal ..."
  mkdir -p "${PRIMER_CHECK_DIR}"
  local out="${PRIMER_CHECK_DIR}/post_filter_summary.tsv"
  if [ ! -f "$out" ]; then echo -e "file\treads_processed\treads_with_adapters\tpct_with_adapters" > "$out"; fi

  shopt -s nullglob
  local F=( "${RESULTS}/filtered/"*.fastq.gz )
  shopt -u nullglob
  [ ${#F[@]} -eq 0 ] && { echo "!!! No filtered FASTQs"; return 1; }

  for ff in "${F[@]}"; do
    b="${ff##*/}"; b="${b%.gz}"; b="${b%.fastq}"; b="${b%.fq}"
    flags=$(cutadapt_primer_flags)
    rpt="${PRIMER_CHECK_DIR}/${b}.cutadapt_post_check.txt"
    eval cutadapt $flags -e "$PRIMER_ERR" --no-trim -j "${THREADS}" -o /dev/null "$ff" > "$rpt"

    reads_proc=$(awk -F': *' '/Total reads processed:/ {gsub(/,/, "", $2); print $2}' "$rpt")
    reads_adp=$(awk -F': *' '/Reads with adapters:/ {gsub(/,/, "", $2); print $2}' "$rpt")
    pct_adp=$(awk -F'[(%)]' '/Reads with adapters:/ {gsub(/ /, "", $2); print $2}' "$rpt")
    [ -z "$reads_proc" ] && reads_proc=NA
    [ -z "$reads_adp" ] && reads_adp=NA
    [ -z "$pct_adp" ] && pct_adp=NA
    echo -e "${b}\t${reads_proc}\t${reads_adp}\t${pct_adp}" >> "$out"
  done
  echo ">>> Post-filter primer check summary -> ${out}"
}

# helper to aggregate per-group length tables to the legacy top-level path
aggregate_group_lengths() {
  echo ">>> Aggregating per-group length tables into legacy combined files ..."
  mkdir -p results/lengths

  #build BOTH pre and post combined TSVs from results/groups/<GROUP>/lengths
  local out_pre="results/lengths/all_lengths.tsv"
  local out_post="results/lengths/all_lengths_post.tsv"
  echo -e "sample\tlength" > "$out_pre"
  echo -e "sample\tlength" > "$out_post"

  local found_pre=0
  local found_post=0

  for g in Archaea Ascomic Bac Basid Unknown; do
    local d="results/groups/${g}/lengths"

    # PRE
    if [ -s "${d}/all_lengths.tsv" ]; then
      awk -v grp="$g" 'NR>1{print grp "/" $1 "\t" $2}' "${d}/all_lengths.tsv" >> "$out_pre"
      echo "    PRE: included ${d}/all_lengths.tsv"
      found_pre=1
    fi

    # POST
    if [ -s "${d}/all_lengths_post.tsv" ]; then
      awk -v grp="$g" 'NR>1{print grp "/" $1 "\t" $2}' "${d}/all_lengths_post.tsv" >> "$out_post"
      echo "    POST: included ${d}/all_lengths_post.tsv"
      found_post=1
    fi
  done

  if [ "$found_pre" -eq 1 ]; then
    echo ">>> Wrote results/lengths/all_lengths.tsv"
  else
    echo "!!! No group PRE length tables found — results/lengths/all_lengths.tsv is empty."
  fi

  if [ "$found_post" -eq 1 ]; then
    echo ">>> Wrote results/lengths/all_lengths_post.tsv"
  else
    echo "!!! No group POST length tables found — results/lengths/all_lengths_post.tsv is empty."
  fi
}

#helper to render combined PRE and POST plots after aggregation
render_combined_plots() {
  if [ ! -f workflow/plot_fastq_lengths.R ]; then
    echo "!!! Missing workflow/plot_fastq_lengths.R — cannot render combined plots."
    return 1
  fi
  mkdir -p results/lengths

  # PRE
  if [ -s results/lengths/all_lengths.tsv ]; then
    echo ">>> Rendering combined PRE plot ..."
    Rscript workflow/plot_fastq_lengths.R
    mv -f results/lengths/read_length_boxplots.png results/lengths/read_length_boxplots_pre.png || true
    mv -f results/lengths/read_length_boxplots.pdf results/lengths/read_length_boxplots_pre.pdf || true
  else
    echo ">>> Skipping PRE plot (no results/lengths/all_lengths.tsv)"
  fi

  # POST
  if [ -s results/lengths/all_lengths_post.tsv ]; then
    echo ">>> Rendering combined POST plot ..."
    cp -f results/lengths/all_lengths_post.tsv results/lengths/all_lengths.tsv
    Rscript workflow/plot_fastq_lengths.R
    mv -f results/lengths/read_length_boxplots.png results/lengths/read_length_boxplots_post.png || true
    mv -f results/lengths/read_length_boxplots.pdf results/lengths/read_length_boxplots_post.pdf || true
    # restore PRE table if it existed (harmless if absent)
    if [ -s results/lengths/all_lengths_pre.tsv ]; then
      cp -f results/lengths/all_lengths_pre.tsv results/lengths/all_lengths.tsv || true
    else
      # or restore the aggregated PRE file if it's still there
      cp -f results/lengths/all_lengths_post.tsv results/lengths/all_lengths.tsv >/dev/null 2>&1 || true
    fi
  else
    echo ">>> Skipping POST plot (no results/lengths/all_lengths_post.tsv)"
  fi
}

# ----------------------------------------------------------
# Function: log_run_report
# ----------------------------------------------------------
log_run_report() {
    [ -d logs ] || mkdir -p logs
    local logfile="logs/run_report_$(date +%Y%m%d_%H%M%S).txt"
    local end_time=$(date +%s); local runtime=$((end_time - START_TIME))
    local minutes=$((runtime / 60)); local seconds=$((runtime % 60))
    local host=$(hostname)
    local env_hash=$(conda env export --name libsQC 2>/dev/null | sha256sum | cut -c1-12)

    {
      echo "=================================================="
      echo " QC RUN REPORT — $(date)"
      echo "=================================================="
      echo "Environment : libsQC"
      echo "System       : ${host}"
      echo "Conda env fingerprint: ${env_hash}"
      echo "Threads used: ${THREADS}"
      echo
      echo "Function timing summary (aggregated):"
      printf "%-30s %10s %8s %10s\n" "Function" "Total(s)" "n" "Avg(s)"
      printf "%-30s %10s %8s %10s\n" "--------" "--------" "--" "------"
      awk -F'\t' '{tot[$1]+=$2; cnt[$1]++} END{for(k in tot){printf "%-30s %10d %8d %10.1f\n", k, tot[k], cnt[k], (cnt[k]? tot[k]/cnt[k] : 0)}}' logs/.timing.tsv 2>/dev/null || true
      echo
      echo "Report saved to: $logfile"
    } > "$logfile"

    echo ">>> Run report complete."
    echo ">>> Total runtime: ${minutes}m ${seconds}s — report saved to ${logfile}"
}

# ----------------------------------------------------------
# Helper: time_function
# ----------------------------------------------------------
time_function() {
    local fn="$1"
    local start=$(date +%s)
    echo ">>> Running $fn ..."
    local status=0
    set +e; $fn; status=$?; set -e
    local end=$(date +%s); local dur=$((end - start))
    [ -d logs ] || mkdir -p logs
    echo -e "${fn}\t${dur}" >> logs/.timing.tsv
    echo ">>> $fn completed in ${dur}s"
    return $status
}

#single-function runner (optional mode)
single_function_runner() {
  local fn="$1"                # name of the bash function to call
  local label="${LABEL:-pre}"  # e.g., pre or post for plotting
  echo ">>> SINGLE-FUNCTION MODE: ${fn}  (LABEL=${label})"

  # we assume demultiplexing already done; operate on existing groups
  for GROUP in Archaea Ascomic Bac Basid Unknown; do
    echo
    echo "================= GROUP: ${GROUP} (ONLY_FN=${fn}) ================="
    RESULTS="results/groups/${GROUP}"
    PRIMER_CHECK_DIR="${RESULTS}/primer_checks"
    PRIMER_TRIM_DIR="${RESULTS}/primer_trimming"

    # pick inputs per function/label
    shopt -s nullglob
    if [[ "${fn}" == "plot_fastq_length_boxplots" ]]; then
      # For plotting: PRE uses raw group data; POST uses filtered
      if [[ "${label}" == "post" ]]; then
        FASTQ_FILES=( "${RESULTS}/filtered/"*.fastq.gz )
      else
        FASTQ_FILES=( "${RESULTS}/data/"*.fastq.gz )
      fi
    else
      # generic: operate on raw group data
      FASTQ_FILES=( "${RESULTS}/data/"*.fastq.gz )
    fi
    shopt -u nullglob

    if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
      echo ">>> No FASTQs for ${GROUP} matching this mode; skipping."
      continue
    fi

    # group-scoped out dir for plotting
    if [[ "${fn}" == "plot_fastq_length_boxplots" ]]; then
      time_function "${fn} ${label} results/groups/${GROUP}/lengths"
    else
      # call the function without extra args
      time_function "${fn}"
    fi
  done

}

# ==========================================================
# Calling the functions (QC pipeline main flow)
# ==========================================================

START_TIME=$(date +%s)
[ -d logs ] || mkdir -p logs
: > logs/.timing.tsv

#if ONLY_FN is set, run just that across all groups and exit
if [[ -n "${ONLY_FN:-}" ]]; then
  # If your function needs the conda env, you can still initialize it:
  time_function create_env_libsQC
  time_function check_versions
  time_function export_env
  # We do NOT demux again; we assume groups already exist.
  single_function_runner "${ONLY_FN}"
  log_run_report
  exit 0
fi

# Environment setup
time_function create_env_libsQC
time_function check_versions
time_function export_env

# Input discovery
time_function gather_fastq_files
# only process first 3 files
#time_function limit_to_three_fastqs

# 3) CLASSIFY reads into groups BEFORE any trimming/QC
# demultiplex by primer pairs
time_function demux_by_primers

# ===== Per-group processing loop ==========================================
# run all subsequent steps per group
for GROUP in Archaea Ascomic Bac Basid Unknown; do
  echo
  echo "================= GROUP: ${GROUP} ================="

  # Set results root for this group
  RESULTS="results/groups/${GROUP}"

  # Pick the group's FASTQs as input for this round
  shopt -s nullglob
  if [[ "$GROUP" == "Unknown" ]]; then
    FASTQ_FILES=( "results/groups/${GROUP}/data/"*.fastq.gz )
  else
    FASTQ_FILES=( "results/groups/${GROUP}/data/"*.fastq.gz )
  fi
  shopt -u nullglob

  if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo ">>> No reads for group ${GROUP}; skipping."
    continue
  fi

  # Set group-specific primer lists for checks & trimming
  case "$GROUP" in
    Archaea) PRIMER_FWD_LIST="${PRIM_FWD_ARCHAEA}"; PRIMER_REV_LIST="${PRIM_REV_ARCHAEA}" ;;
    Ascomic) PRIMER_FWD_LIST="${PRIM_FWD_ASCOMIC}"; PRIMER_REV_LIST="${PRIM_REV_ASCOMIC}" ;;
    Bac)     PRIMER_FWD_LIST="${PRIM_FWD_BAC}";     PRIMER_REV_LIST="${PRIM_REV_BAC}"     ;;
    Basid)   PRIMER_FWD_LIST="${PRIM_FWD_BASID}";   PRIMER_REV_LIST="${PRIM_REV_BASID}"   ;;
    Unknown) PRIMER_FWD_LIST=""; PRIMER_REV_LIST="";;
  esac

  # Update primer dirs to group-scoped locations
  PRIMER_CHECK_DIR="${RESULTS}/primer_checks"
  PRIMER_TRIM_DIR="${RESULTS}/primer_trimming"

  # Initial QC on raw data (per group)
  time_function run_fastqc_all
  time_function run_multiqc

  # Nanopore-specific QC (pre-filter)
  time_function 'quick_len_qual_overview pre'
  time_function primer_spotcheck

  # Summaries and QC flags
  time_function make_fastq_summary
  time_function qc_flags_from_nanoplot
  time_function "plot_fastq_length_boxplots pre results/groups/${GROUP}/lengths"

  # Filtering step (trim primers + NanoFilt) and verify
  time_function filter_amplicons
  time_function verify_primer_removal
  time_function re_qc_filtered
done

#Final report
log_run_report
