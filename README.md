<!-- LOGO -->
<p align="center">
  <img src="assets/logo-horizontal.png" alt="AY:RΔ Logo" width="160"/>
</p>

<h1 align="center" style="font-weight: normal;">metatax-amazon</h1>

<p align="center">
  <img src="https://img.shields.io/badge/pipeline-HPC--native-green"/>
  <img src="https://img.shields.io/badge/data-Nanopore-orange"/>
  <img src="https://img.shields.io/badge/workflow-metataxonomy-blueviolet"/>
  <img src="https://img.shields.io/badge/markers-16S%20%7C%20ITS%20%7C%20LSU-teal"/>
  <img src="https://img.shields.io/badge/QC-marker--aware-darkgreen"/>
  <img src="https://img.shields.io/badge/classifier-Emu-blue"/>
  <img src="https://img.shields.io/badge/license-MIT-blue"/>
</p>

<p align="center">
  <code>data and discovery in flow</code><br/>
  <a href="mailto:ayrabioinf@gmail.com">ayrabioinf@gmail.com</a> · 
  <a href="https://www.linkedin.com/company/aryaiam">linkedin here</a>
</p>

---
---

<pre>
ABOUT
-----
metatax-amazon is a metataxonomy analysis pipeline developed by AY:RΔ.
It integrates 16S(archaea and bacteria)/ITS/LSU amplicon data processing, taxonomic classification,
and downstream diversity metrics for environmental microbiome studies.

This project was designed for large-scale Amazonian soil and water datasets,
enabling reproducible exploration of microbial community structure and composition.
</pre>

<pre>
WORKFLOW
--------
1. FASTQ discovery from data/
2. Marker-aware primer trimming with Cutadapt
3. Pre-filter QC on primer-trimmed reads
4. Read-length diagnostics and QC flag generation
5. Optional marker-aware NanoFilt filtering
6. Post-filter QC on filtered reads
7. Emu taxonomic classification
8. Abundance table generation
9. Downstream diversity and differential abundance analysis
</pre>

<pre>
QC DIAGNOSTIC MODE
------------------
The QC diagnostic mode performs primer trimming first, then evaluates read quality
and read-length distributions using the primer-trimmed FASTQs.

This mode is useful before applying NanoFilt because it allows empirical inspection
of marker-specific read-length profiles.

Example:
  bash workflow/runall.sh \
    --qc-length-diagnostic-only

This runs:
  - marker-aware primer trimming
  - FastQC on trimmed reads
  - MultiQC on trimmed reads
  - NanoStat / NanoPlot on trimmed reads
  - SeqKit summary statistics
  - read-length boxplots
  - QC flag generation

And stops before:
  - NanoFilt filtering
  - post-filter QC
  - Emu
  - downstream analysis

Main outputs:
  results/<batch>/trimmed/
  results/<batch>/untrimmed/
  results/<batch>/trim_reports/
  results/<batch>/qc_raw/
  results/<batch>/multiqc/
  results/<batch>/nanoplot/
  results/<batch>/lengths/
  results/<batch>/summary/
</pre>

<pre>
MARKER-AWARE QC
---------------
QC is marker-aware and based on marker patterns detected in FASTQ file names.
Supported marker labels include:

  - 16SA  archaeal 16S
  - 16SB  bacterial 16S
  - ITS
  - LSU

Examples:
  nanopore_amplicon_SAMPLE_16SA-LSU.fastq.gz
  → archaeal 16S + LSU primers are used

  nanopore_amplicon_SAMPLE_16SB-ITS.fastq.gz
  → bacterial 16S + ITS primers are used

Primer trimming is performed with Cutadapt using the primer sets selected from
the FASTQ name. The pipeline records, for each library:

  - detected markers
  - forward primers used
  - reverse primers used
  - Cutadapt report path
  - exact Cutadapt command inside the report

Reads without detected primer matches are not discarded. They are written to:

  results/<batch>/untrimmed/

Primer-trimmed reads are written to:

  results/<batch>/trimmed/

The pre-filter QC plots and read-length boxplots are generated from these
primer-trimmed reads, not from the original raw FASTQs.
</pre>

<pre>
OPTIONAL NANOFILT FILTERING
---------------------------
NanoFilt filtering is not applied automatically after diagnostic QC.
To continue from trimmed-read diagnostics into quality/length filtering, use:

  bash workflow/runall.sh \
    --qc-run-filtering \
    --no-emu \
    --no-downstream

NanoFilt uses the primer-trimmed FASTQs as input and writes filtered FASTQs to:

  results/<batch>/filtered/

Filtering is marker-aware. Current default cutoffs are:

  ITS-only libraries:
    mean Q >= 10
    length 150–1000 bp

  16S, LSU, or mixed-marker libraries:
    mean Q >= 10
    length 150–1800 bp

For each library, the exact NanoFilt parameters are recorded in:

  results/<batch>/summary/nanofilt_cutoffs_by_library.tsv

After filtering, the pipeline runs post-filter QC:

  - FastQC on filtered reads
  - MultiQC on filtered reads
  - NanoStat / NanoPlot on filtered reads
  - SeqKit summary statistics
  - post-filter read-length boxplots
  - QC flag generation
</pre>
    
<pre>
STRUCTURE (IMPROVE!!)
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
DEPENDENCIES (IMPROVE!!)
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

# 1) Run QC diagnostics only (trim primers, inspect read lengths, stop)
#     Example: first 25 FASTQs (0–24)
BATCH_TAG=b000_n025 \
LIMIT_FASTQS=25 \
OFFSET_FASTQS=0 \
bash workflow/runall.sh \
  --qc-length-diagnostic-only

# 2) Run QC filtering after inspecting diagnostics
#     Example: same 25 FASTQs (0–24)
BATCH_TAG=b000_n025_filter \
LIMIT_FASTQS=25 \
OFFSET_FASTQS=0 \
bash workflow/runall.sh \
  --qc-run-filtering \
  --no-emu \
  --no-downstream

# 7) Run ONLY Emu (QC already done) on 16S (FIX FROM HERE!!)
bash workflow/runall.sh --no-qc --no-downstream

# 8) Run Emu on ALL filtered FASTQs
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
 results/<batch>/qc_raw/               - FastQC reports for primer-trimmed reads
 results/<batch>/qc_filtered/          - FastQC reports for NanoFilt-filtered reads
 results/<batch>/multiqc/              - MultiQC reports
 results/<batch>/nanoplot/             - NanoPlot outputs
 results/<batch>/summary/              - QC summaries and flags
   results/<batch>/summary/primer_trimming_by_library.tsv
   results/<batch>/summary/nanofilt_cutoffs_by_library.tsv
   results/<batch>/summary/qc_flags.tsv
   results/<batch>/summary/seqkit_stats.tsv
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


