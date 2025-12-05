<!-- LOGO -->
<p align="center">
  <img src="https://github.com/ayrdelta/.github/blob/main/profile/logo.png" alt="AY:RΔ Logo" width="120"/>
</p>

<h1 align="center" style="font-weight: normal;">metatax-amazon</h1>

<p align="center">
  <code>data and discovery in flow</code><br/>
  <a href="mailto:ayrabioinf@gmail.com">ayrabioinf@gmail.com</a> · 
  <a href="https://www.linkedin.com/company/aryaiam">linkedin here</a>
</p>

---

<pre>
ABOUT
-----
metatax-amazon is a metataxonomy analysis pipeline developed by AY:RΔ.
It integrates 16S/ITS amplicon data processing, taxonomic classification,
and downstream diversity metrics for environmental microbiome studies.

This project was designed for large-scale Amazonian soil and water datasets,
enabling reproducible exploration of microbial community structure and composition.
</pre>

<pre>
STRUCTURE
---------
 /workflow/     - main pipeline scripts
 /envs/         - Conda/Mamba environment YAMLs
 /metadata/     - sample and sequencing metadata tables
 /config/       - configuration files for parameters and paths
 /docs/         - usage and documentation
 LICENSE        - project license (MIT)
 CITATION.cff   - citation metadata for referencing this work
 bootstrap.sh   - environment and directory setup script
 README.md      - this file
</pre>

<pre>
DEPENDENCIES
------------
 - Python ≥ 3.10
 - R ≥ 4.3
 - Conda or Mamba
 - Slurm (or compatible scheduler)
 - NanoPlot · FastQC · MultiQC · SeqKit · Cutadapt · NanoFilt · Emu
</pre>

<pre>
USAGE
-----
After cloning the repository and running "bash bootstrap.sh", you can execute
the pipeline via the main launcher:

  bash workflow/runall.sh [options]

MAIN OPTIONS
-------------
  --partition STR       Partition/queue name (default: short)
  --time HH:MM:SS       Walltime (default: 04:00:00)
  --cpus INT            CPUs per task (default: 4)
  --mem STR             Memory (default: 16G)
  --wd PATH             Working directory (default: current directory)
  --primer-fwd SEQ      Forward primer sequence (optional)
  --primer-rev SEQ      Reverse primer sequence (optional)
  --seq-summary PATH    Path to sequencing_summary.txt (optional)
  --no-qc               Skip libsQC (run Emu only)
  --no-emu              Skip Emu (run libsQC only)

EMU OPTIONS
------------
  --emu-partition STR   Partition for Emu (default: inherit libsQC)
  --emu-time HH:MM:SS   Walltime for Emu (default: inherit libsQC)
  --emu-cpus INT        CPUs for Emu (default: inherit libsQC)
  --emu-mem STR         Memory for Emu (default: inherit libsQC)
  --emu-db-its PATH     Path to Emu ITS database (optional)
  --emu-db-lsu PATH     Path to Emu LSU database (optional)

If ITS/LSU DBs are not provided, only the 16S step will be performed.
The 16S Emu database (bacteria + archaea) is auto-downloaded at first use.

<pre>
EXAMPLES:
---------

# 1) Run ALL stages (QC → Emu → Downstream) with default settings
bash workflow/runall.sh

# 2) Run ALL stages with custom resources
bash workflow/runall.sh --time 08:00:00 --cpus 8 --mem 32G

# 3) Build ONLY the ITS/LSU marker databases (UNITE ITS + SILVA LSU)
#    Uses default locations unless ITS_FASTA / LSU_FASTA are overridden:
#      ITS_FASTA=/path/to/sh_general_release_dynamic_19.02.2025.fasta
#      LSU_FASTA=/path/to/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
#
#    This step does NOT run QC, Emu, or downstream analysis.
bash workflow/runall.sh --only-build-marker-dbs

# 4) Run ONLY the QC stage
bash workflow/runall.sh --no-emu --no-downstream \
  --time 08:00:00 --cpus 16 --mem 32G


# 5) Run ONLY Emu (QC already done)
#    By default, processes only the first 3 FASTQs for a quick test
bash workflow/runall.sh --no-qc --no-downstream \
  --emu-time 12:00:00 --emu-cpus 20 --emu-mem 32G


# 6) Run Emu on ALL FASTQs (disable the default 3-file test limit)
LIMIT_FASTQS=0 bash workflow/runall.sh --no-qc --no-downstream \
  --emu-time 12:00:00 --emu-cpus 20 --emu-mem 32G


# 7) Run QC + Emu, but give Emu more resources
bash workflow/runall.sh \
  --time 06:00:00 --cpus 8 --mem 32G \
  --emu-time 12:00:00 --emu-cpus 16 --emu-mem 64G


# 8) Run Emu with a custom ITS database
bash workflow/runall.sh --no-qc --no-downstream \
  --emu-db-its /path/to/emu_unite_its


# 9) Run Emu with a custom LSU database
bash workflow/runall.sh --no-qc --no-downstream \
  --emu-db-lsu /path/to/emu_silva_lsu


# 10) Run Emu on a custom FASTQ directory (skip QC)
FASTQ_DIR_DEFAULT=/path/to/custom_filtered \
bash workflow/runall.sh --no-qc --no-downstream \
  --emu-time 10:00:00 --emu-cpus 12 --emu-mem 32G


# 11) Run Emu in batches (recommended for very large datasets)
#     LIMIT_FASTQS = number of FASTQs to process per batch
#     --offset-fastqs = number of FASTQs to skip from the start
#
#     Each batch creates its own directories:
#        results/emu_runs_bXXX_nYYY/
#        results/tables_bXXX_nYYY/
#        results/plots_bXXX_nYYY/

# ---- First 25 files (0–24)
LIMIT_FASTQS=25 bash workflow/runall.sh --no-qc --no-downstream \
  --offset-fastqs 0 \
  --emu-time 05:00:00 --emu-cpus 20 --emu-mem 32G

# ---- Next 25 files (25–49)
LIMIT_FASTQS=25 bash workflow/runall.sh --no-qc --no-downstream \
  --offset-fastqs 25 \
  --emu-time 05:00:00 --emu-cpus 20 --emu-mem 32G

# ---- Last 26 files (50–75)
LIMIT_FASTQS=26 bash workflow/runall.sh --no-qc --no-downstream \
  --offset-fastqs 50 \
  --emu-time 05:00:00 --emu-cpus 20 --emu-mem 32G

# (Optional) Save JSON dictionaries from Emu
SAVE_JSON=1 LIMIT_FASTQS=25 bash workflow/runall.sh --no-qc --no-downstream

# (Optional) Save read-assignment matrices (large files)
SAVE_ASSIGN=1 LIMIT_FASTQS=25 bash workflow/runall.sh --no-qc --no-downstream

# ---------------------------------------------------------
# 12) Run ONLY the downstream diversity analysis
# ---------------------------------------------------------
# NOTE:
#   If the file you provide (via --downstream-infile) exists,
#   the R script will use it directly.
#
#   If it does NOT exist, but you have batch folders such as:
#       results/tables_b000_n025/abundance_combined.tsv
#       results/tables_b025_n025/abundance_combined.tsv
#       ...
#   then downstream_analysis.R will automatically:
#       - detect all batch abundance files
#       - merge them into one master table
#       - write it to: results/tables/abundance_combined.tsv
#       - then continue normally.

bash workflow/runall.sh --no-qc --no-emu \
  --downstream-infile results/tables/abundance_combined.tsv \
  --downstream-outdir results/plots \
  --downstream-basename downstream

# 13) Run QC + Emu, but skip downstream analysis
bash workflow/runall.sh --no-downstream
</pre>
  
<pre>
OUTPUTS [TO BE UPDATED]
-------
 logs/                - timestamped job logs (.out / .err)
 results/qc_raw/      - FastQC reports
 results/multiqc/     - MultiQC aggregated reports
 results/filtered/    - post-filtered FASTQs
 results/emu_runs/    - per-sample Emu abundance outputs
 metadata/bin_counts.tsv  - per-bin read counts manifest
</pre>

<pre>
CITATION
--------
If you use this pipeline, please cite:

Lobo, I. (2025).
Metatax-Amazon: A reproducible long-read metataxonomy pipeline
for Amazonian soil microbiome analysis.
AY:RΔ - data and discovery in flow.
https://github.com/ayraiam/metatax-amazon

See CITATION.cff for full citation metadata.
</pre>

<pre>
CONTACT
-------
ayrabioinf@gmail.com
https://www.linkedin.com/company/ayraiam
</pre>

email: ayrabioinf@gmail.com
linkedin: https://www.linkedin.com/company/ayraiam
</pre>

<p align="center"><sub>© 2025 AY:RΔ — data and discovery in flow</sub></p>
