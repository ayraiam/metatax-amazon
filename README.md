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
It integrates 16S/ITS/LSU amplicon data processing, taxonomic classification,
and downstream diversity metrics for environmental microbiome studies.

This project was designed for large-scale Amazonian soil and water datasets,
enabling reproducible exploration of microbial community structure and composition.
</pre>

<pre>
MARKER-AWARE QC
---------------
The QC stage performs marker-aware primer trimming based on FASTQ filenames.

Libraries may contain different combinations of:
  - 16SA  (archaeal 16S)
  - 16SB  (bacterial 16S)
  - ITS
  - LSU

Example:
  nanopore_amplicon_SAMPLE_16SA-LSU.fastq.gz
  → archaeal 16S primers + LSU primers are used

  nanopore_amplicon_SAMPLE_16SB-ITS.fastq.gz
  → bacterial 16S primers + ITS primers are used

Primer trimming is performed dynamically per library using Cutadapt.

For each FASTQ, the pipeline records:
  - detected marker groups
  - forward/reverse primers used
  - exact Cutadapt command
  - trimming reports

Outputs:
  results/<batch>/trim_reports/
  results/<batch>/summary/primer_trimming_by_library.tsv
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
 - R ≥ 4.5                    
 - Conda or Mamba
 - Slurm (or compatible scheduler)
 - NanoPlot · FastQC · MultiQC · SeqKit · Cutadapt · NanoFilt · Emu
</pre>

<pre>
BATCH PROCESSING
----------------
Large sequencing runs can be processed in batches using:

  LIMIT_FASTQS
  OFFSET_FASTQS
  BATCH_TAG

This prevents overwriting outputs and enables scalable execution
across large ONT amplicon datasets.

Example:
  BATCH_TAG=b000_n050
  LIMIT_FASTQS=50
  OFFSET_FASTQS=0

Outputs become batch-specific:

  results/b000_n050/
  logs/b000_n050_*

Recommended for:
  - large sequencing projects
  - Slurm cluster execution
  - incremental QC / Emu processing
  - recovery from interrupted runs
</pre>

<pre>
MAIN OPTIONS
-------------
  --partition STR       Partition/queue name (default: short)
  --time HH:MM:SS       Walltime (default: 04:00:00)
  --cpus INT            CPUs per task (default: 4)
  --mem STR             Memory (default: 16G)
  --wd PATH             Working directory (default: current directory)

  --batch-tag STR      Batch identifier used for results and logs
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
  ### Steps 7–8 (ANCOM-BC2) are controlled independently (see below)  

  --use-counts-0-4 INT  Numeric column for Parts 0–4 (default: 0)
                          0 = abundance
                          1 = estimated_counts

  --use-counts-5 INT    Numeric column for Part 5 (default: 1)
                          1 = estimated_counts (recommended for CLR)
                          0 = abundance

  --use-counts-ancom INT Numeric column for Steps 7–8 (ANCOM-BC2) (default: 1)  
                          1 = estimated_counts (recommended; ANCOM-BC2 expects counts)  
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
</pre>

<pre>
DOWNSTREAM NOTES (IMPORTANT)                                   
-------------------------
- Parts 0–4 default to `abundance` (compositional-style plots + diversity metrics).
- Part 5 defaults to `estimated_counts` (recommended for CLR concordance).
- Steps 7–8 (ANCOM-BC2) default to `estimated_counts` via `--use-counts-ancom 1`,
  because ANCOM-BC2 is a count-based model and expects counts rather than relative abundance.  
- Step 9 generates a pheatmap heatmap for ANCOM-BC2 significant genera (CLR-transformed).  
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
(Defaults: Parts 0–4 use abundance; Part 5 uses estimated_counts; ANCOM-BC2 uses estimated_counts)  

# 2) Run ALL stages with custom resources
bash workflow/runall.sh --time 08:00:00 --cpus 8 --mem 32G

# 3) Build ONLY the ITS/LSU marker databases
bash workflow/runall.sh --only-build-marker-dbs

# 4) Run ONLY the QC stage
bash workflow/runall.sh --no-emu --no-downstream

# QC in batches (recommended for large projects)

# 5) First 50 libraries
BATCH_TAG=b000_n050 \
LIMIT_FASTQS=50 OFFSET_FASTQS=0 \
bash workflow/runall.sh \
  --no-emu --no-downstream

# 6) Next 50 libraries
BATCH_TAG=b050_n050 \
LIMIT_FASTQS=50 OFFSET_FASTQS=50 \
bash workflow/runall.sh \
  --no-emu --no-downstream

# 7) Run ONLY Emu (QC already done) on 16S
bash workflow/runall.sh --no-qc --no-downstream

# 8) Run Emu on ALL FASTQs (disable 3-file test limit) for 16S
LIMIT_FASTQS=0 bash workflow/runall.sh --no-qc --no-downstream

# 9) Run QC + Emu, but give Emu more resources
bash workflow/runall.sh \
  --time 06:00:00 --cpus 8 --mem 32G \
  --emu-time 12:00:00 --emu-cpus 16 --emu-mem 64G

# 10) Run Emu with a custom ITS DB
bash workflow/runall.sh --no-qc --no-downstream \
  --emu-db-its /path/to/its_db

# 11) Run Emu with a custom LSU DB
bash workflow/runall.sh --no-qc --no-downstream \
  --emu-db-lsu /path/to/lsu_db

# 12) Run Emu on a custom FASTQ directory (skip QC)
FASTQ_DIR_DEFAULT=/path/to/filtered \
bash workflow/runall.sh --no-qc --no-downstream

# 13) Run ONLY ITS + LSU (skip 16S) on ALL FASTQs
ENABLE_16S=0 ENABLE_ITS=1 ENABLE_LSU=1 \
bash workflow/runall.sh --no-qc --no-downstream

# 14) Run ONLY 16S (explicit)
ENABLE_16S=1 ENABLE_ITS=0 ENABLE_LSU=0 \
bash workflow/runall.sh --no-qc --no-downstream

# 15) Run ITS ONLY in batches (recommended for large datasets)
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

# 16) Run ONLY the downstream analysis (16S)
bash workflow/runall.sh --no-qc --no-emu --mode 16S

#     Run ONLY the downstream analysis (ITS)
bash workflow/runall.sh --no-qc --no-emu --mode ITS

# 17) Run QC + Emu but skip downstream analysis
bash workflow/runall.sh --no-downstream

# 18) Override numeric sources explicitly (rare / advanced):
#     - Parts 0–4: stacked bars + alpha/beta
#     - Part 5: CLR concordance
#     - Steps 7–8: ANCOM-BC2 differential abundance
bash workflow/runall.sh --mode 16S \
  --use-counts-0-4 1 \
  --use-counts-5 0 \
  --use-counts-ancom 1                                  
</pre>

<pre>
OUTPUTS
-------
 logs/<batch>_*                         - batch-specific job logs
 results/<batch>/trimmed/              - primer-trimmed FASTQs
 results/<batch>/untrimmed/            - reads without detected primers
 results/<batch>/filtered/             - post-filtered FASTQs
 results/<batch>/trim_reports/         - Cutadapt reports per library
 results/<batch>/qc_raw/               - FastQC reports
 results/<batch>/multiqc/              - MultiQC reports
 results/<batch>/nanoplot/             - NanoPlot outputs
 results/<batch>/summary/              - QC summaries and flags
 results/<batch>/lengths/              - read-length distributions
 results/emu_runs_*                    - Emu outputs
 results/tables_*                      - merged abundance tables
 results/plots_*                       - downstream figures
 metadata/                             - FASTQ manifests and primer lists
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


