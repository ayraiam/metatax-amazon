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

  --no-qc               Skip libsQC
  --no-emu              Skip Emu
  --no-downstream       Skip downstream analysis

  --no-build-marker-dbs   Do not attempt to build ITS/LSU databases
  --only-build-marker-dbs Build ITS/LSU databases **only**
                          (no QC, no Emu, no downstream)


MARKER SELECTION (NEW)
-----------------------
The pipeline can now run **any combination** of 16S / ITS / LSU:

  --enable-16s [0|1]    Enable or disable 16S analysis   (default: 1)
  --enable-its [0|1]    Enable or disable ITS analysis   (default: auto: 1 if ITS DB exists)
  --enable-lsu [0|1]    Enable or disable LSU analysis   (default: auto: 1 if LSU DB exists)

Examples:
  --enable-16s 0 --enable-its 1 --enable-lsu 1     ⟶ ITS + LSU only
  --enable-16s 1 --enable-its 0 --enable-lsu 0     ⟶ 16S only
  --enable-16s 1 --enable-its 1 --enable-lsu 1     ⟶ run all markers


EMU OPTIONS
------------
  --emu-partition STR   Partition for Emu (default: inherit libsQC)
  --emu-time HH:MM:SS   Walltime for Emu
  --emu-cpus INT        CPUs for Emu
  --emu-mem STR         Memory for Emu

  --emu-db-its PATH     Path to ITS Emu DB directory
  --emu-db-lsu PATH     Path to LSU Emu DB directory

If ITS/LSU DBs are not provided or missing, their respective markers are skipped.
The 16S Emu database (bacteria + archaea) is auto-downloaded at first use.


EXAMPLES
---------

# 1) Run ALL stages (QC → Emu → Downstream) for ALL markers
bash workflow/runall.sh


# 2) Run ALL stages with custom resources
bash workflow/runall.sh --time 08:00:00 --cpus 8 --mem 32G


# 3) Build ONLY the ITS/LSU marker databases
bash workflow/runall.sh --only-build-marker-dbs


# 4) Run ONLY the QC stage
bash workflow/runall.sh --no-emu --no-downstream


# 5) Run ONLY Emu (QC already done)
bash workflow/runall.sh --no-qc --no-downstream


# 6) Run Emu on ALL FASTQs (disable 3-file test limit)
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


# 11) Run ONLY ITS + LSU (skip 16S)
bash workflow/runall.sh \
  --enable-16s 0 \
  --enable-its 1 \
  --enable-lsu 1


# 12) Run ONLY 16S
bash workflow/runall.sh \
  --enable-16s 1 \
  --enable-its 0 \
  --enable-lsu 0


# 13) Run Emu in batches (recommended for large datasets)
LIMIT_FASTQS=25 bash workflow/runall.sh --no-qc --no-downstream --offset-fastqs 0
LIMIT_FASTQS=25 bash workflow/runall.sh --no-qc --no-downstream --offset-fastqs 25
LIMIT_FASTQS=26 bash workflow/runall.sh --no-qc --no-downstream --offset-fastqs 50


# 14) Run ONLY the downstream analysis
bash workflow/runall.sh --no-qc --no-emu \
  --downstream-infile results/tables/abundance_combined.tsv


# 15) Run QC + Emu but skip downstream analysis
bash workflow/runall.sh --no-downstream
</pre>

<pre>
OUTPUTS
-------
 logs/                    - timestamped job logs
 results/filtered/        - post-filtered FASTQs
 results/emu_runs_bXXX/   - per-marker Emu output (16S / ITS / LSU separated)
 results/tables_bXXX/     - merged abundance & mapping tables
 results/plots_bXXX/      - genus-level stacked barplots
 metadata/                - FASTQ manifest, JSON dicts, primer lists

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
