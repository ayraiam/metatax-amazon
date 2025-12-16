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
After cloning the repository and running `bash bootstrap.sh`, you can execute
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

  --no-qc               Skip libsQC
  --no-emu              Skip Emu
  --no-downstream       Skip downstream analysis

  --no-build-marker-dbs   Do not attempt to build ITS/LSU databases
  --only-build-marker-dbs Build ITS/LSU databases **only**
                          (no QC, no Emu, no downstream)

  ### Downstream mode and numeric source
  --mode STR            Downstream marker mode: 16S or ITS (default: 16S)
                        Selects the default downstream input file:
                          16S -> results/tables/abundance_combined.tsv
                          ITS -> results/tables_ITS/abundance_combined.tsv

  ### Numeric source selection (split behavior)
  ### Parts 0–4 control stacked bars + alpha/beta diversity
  ### Part 5 controls CLR concordance scatter plots
  --use-counts-0-4 INT  Numeric column for Parts 0–4 (default: 0)
                          0 = abundance
                          1 = estimated_counts

  --use-counts-5 INT    Numeric column for Part 5 (default: 1)
                          1 = estimated_counts (recommended for CLR)
                          0 = abundance
</pre>

<pre>
MARKER SELECTION (ENV VARS)
---------------------------
The pipeline can run **any combination** of 16S / ITS / LSU via
environment variables that are passed through to `run_emu_amplicons.sh`:

  ENABLE_16S=0|1    Enable / disable 16S analysis   (default: 1)
  ENABLE_ITS=0|1    Enable / disable ITS analysis   (default: 0, but auto-enabled if ITS DB exists and you set it to 1)
  ENABLE_LSU=0|1    Enable / disable LSU analysis   (default: 0, but auto-enabled if LSU DB exists and you set it to 1)

### NOTE
These variables control the **Emu stage only**.
The downstream analysis is selected independently using `--mode 16S|ITS`.

Examples:

  # ITS + LSU only
  ENABLE_16S=0 ENABLE_ITS=1 ENABLE_LSU=1 bash workflow/runall.sh ...

  # 16S only
  ENABLE_16S=1 ENABLE_ITS=0 ENABLE_LSU=0 bash workflow/runall.sh ...

  # All markers (assuming ITS / LSU DBs exist)
  ENABLE_16S=1 ENABLE_ITS=1 ENABLE_LSU=1 bash workflow/runall.sh ...
</pre>

<pre>
EMU OPTIONS
------------
  --emu-partition STR   Partition for Emu (default: inherit libsQC)
  --emu-time HH:MM:SS   Walltime for Emu
  --emu-cpus INT        CPUs for Emu
  --emu-mem STR         Memory for Emu

  --emu-db-its PATH     Path to ITS Emu DB directory
  --emu-db-lsu PATH     Path to LSU Emu DB directory

If ITS/LSU DBs are not provided or are missing, those markers are skipped.
The 16S Emu database (bacteria + archaea) is auto-downloaded at first use.
</pre>

<pre>
EXAMPLES
---------

# 1) Run ALL stages (QC → Emu → Downstream) for 16S only (default)
bash workflow/runall.sh --mode 16S
(Defaults: Parts 0–4 use abundance; Part 5 uses estimated_counts)

# 2) Run ALL stages with custom resources
bash workflow/runall.sh --time 08:00:00 --cpus 8 --mem 32G

# 3) Build ONLY the ITS/LSU marker databases
bash workflow/runall.sh --only-build-marker-dbs

# 4) Run ONLY the QC stage
bash workflow/runall.sh --no-emu --no-downstream

# 5) Run ONLY Emu (QC already done) on 16S
bash workflow/runall.sh --no-qc --no-downstream

# 6) Run Emu on ALL FASTQs (disable 3-file test limit) for 16S
LIMIT_FASTQS=0 bash workflow/runall.sh --no-qc --no-downstream

# 7) Run QC + Emu, but give Emu more resources
bash workflow/runall.sh \
  --time 06:00:00 --cpus 8 --mem 32G \
  --emu-time 12:00:00 --emu-cpus 16 --emu-mem 64G

# 8) Run Emu with a custom ITS DB
bash workflow/runall.sh --no-qc --no-downstream \
  --emu-db-its /path/to/its_db

# 9) Run Emu with a custom LSU DB
bash workflow/runall.sh --no-qc --no-downstream \
  --emu-db-lsu /path/to/lsu_db

# 10) Run Emu on a custom FASTQ directory (skip QC)
FASTQ_DIR_DEFAULT=/path/to/filtered \
bash workflow/runall.sh --no-qc --no-downstream

# 11) Run ONLY ITS + LSU (skip 16S) on ALL FASTQs
ENABLE_16S=0 ENABLE_ITS=1 ENABLE_LSU=1 \
bash workflow/runall.sh --no-qc --no-downstream

# 12) Run ONLY 16S (explicit)
ENABLE_16S=1 ENABLE_ITS=0 ENABLE_LSU=0 \
bash workflow/runall.sh --no-qc --no-downstream

# 13) Run ITS ONLY in batches (recommended for large datasets)
#     Example: first 25 FASTQs (0–24)
ENABLE_16S=0 ENABLE_ITS=1 ENABLE_LSU=0 \
BATCH_TAG=bITS_b000_n025 \
LIMIT_FASTQS=25 OFFSET_FASTQS=0 \
FASTQ_DIR_DEFAULT=results/filtered \
bash workflow/runall.sh --no-qc --no-downstream \
  --emu-time 05:00:00 --emu-cpus 20 --emu-mem 32G

#     Next 25 FASTQs (25–49)
ENABLE_16S=0 ENABLE_ITS=1 ENABLE_LSU=0 \
BATCH_TAG=bITS_b025_n025 \
LIMIT_FASTQS=25 OFFSET_FASTQS=25 \
FASTQ_DIR_DEFAULT=results/filtered \
bash workflow/runall.sh --no-qc --no-downstream \
  --emu-time 05:00:00 --emu-cpus 20 --emu-mem 32G

# 14) Run ONLY the downstream analysis (16S)
bash workflow/runall.sh --no-qc --no-emu --mode 16S

#     Run ONLY the downstream analysis (ITS)
bash workflow/runall.sh --no-qc --no-emu --mode ITS

# 15) Run QC + Emu but skip downstream analysis
bash workflow/runall.sh --no-downstream

# 16) Override numeric sources explicitly (rare / advanced):
#     - Parts 0–4: stacked bars + alpha/beta
#     - Part 5: CLR concordance
bash workflow/runall.sh --mode 16S \
  --use-counts-0-4 1 \
  --use-counts-5 0
</pre>

<pre>
OUTPUTS
-------
 logs/                             - timestamped job logs
 results/filtered/                 - post-filtered FASTQs
 results/emu_runs_bXXX/            - per-marker Emu output (16S / ITS / LSU separated)
 results/tables_bXXX/              - merged abundance & mapping tables
 results/plots_bXXX/               - genus-level stacked barplots
 results/plots/*_code_concordance*/- per-code genus CLR concordance scatter plots   ### NEW
 metadata/                         - FASTQ manifest, JSON dicts, primer lists
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

<p align="center"><sub>© 2025 AY:RΔ — data and discovery in flow</sub></p>

