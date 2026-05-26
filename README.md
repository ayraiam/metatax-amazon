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
7. Marker-aware Emu classification (16S / ITS / LSU independently)
8. Marker-specific abundance collation and mapping statistics
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
  results/[batch]/trimmed/
  results/[batch]/untrimmed/
  results/[batch]/trim_reports/
  results/[batch]/qc_raw/
  results/[batch]/multiqc/
  results/[batch]/nanoplot/
  results/[batch]/lengths/
  results/[batch]/summary/
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

  results/[batch]/untrimmed/

Primer-trimmed reads are written to:

  results/[batch]/trimmed/

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

  results/[batch]/filtered/

Filtering is marker-aware. Current default cutoffs are:

  ITS-only libraries:
    mean Q >= 10
    length 150–1000 bp

  16S, LSU, or mixed-marker libraries:
    mean Q >= 10
    length 150–1800 bp

For each library, the exact NanoFilt parameters are recorded in:

  results/[batch]/summary/nanofilt_cutoffs_by_library.tsv

After filtering, the pipeline runs post-filter QC:

  - FastQC on filtered reads
  - MultiQC on filtered reads
  - NanoStat / NanoPlot on filtered reads
  - SeqKit summary statistics
  - post-filter read-length boxplots
  - QC flag generation
</pre>
    
<pre>
STRUCTURE
---------
 /workflow/     - core pipeline scripts
 /envs/         - Conda/Mamba environment YAMLs
 /metadata/     - manifests, primers, and run metadata
 /results/      - QC, Emu, and downstream outputs
 /logs/         - batch-specific execution logs
 LICENSE        - project license (MIT)
 CITATION.cff   - citation metadata
 bootstrap.sh   - environment setup helper
 README.md      - documentation and usage
</pre>

<pre>
DEPENDENCIES
------------
Runtime:
  - Python >= 3.10
  - R >= 4.3
  - Conda or Mamba (recommended)

Execution backend:
  - Slurm (recommended for HPC execution)
    or a compatible scheduler/environment capable of running shell jobs

Core QC tools:
  - Cutadapt
  - NanoFilt
  - NanoPlot / NanoStat
  - FastQC
  - MultiQC
  - SeqKit

Taxonomic classification:
  - Emu
  - minimap2

Python packages:
  - pandas
  - numpy

R packages:
  - ggplot2
  - data.table

Environment management:
  - The pipeline automatically creates and exports Conda environments
    when needed (see /envs/*.yml)

Notes:
  - Emu dependencies are installed automatically through
    workflow/run_emu_amplicons.sh
  - ITS/LSU reference databases can be built automatically through
    workflow/run_build_ITS_LSU_dbs.sh
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
MARKER SELECTION (EMU STAGE)
----------------------------
The Emu step can run any combination of:

  - 16S
  - ITS
  - LSU

Marker selection is controlled through environment variables passed to
`run_emu_amplicons.sh`:

  ENABLE_16S=0|1
  ENABLE_ITS=0|1
  ENABLE_LSU=0|1

Defaults:

  ENABLE_16S=1
  ENABLE_ITS=0
  ENABLE_LSU=0

Examples:

  # 16S only
  ENABLE_16S=1 ENABLE_ITS=0 ENABLE_LSU=0

  # ITS only
  ENABLE_16S=0 ENABLE_ITS=1 ENABLE_LSU=0

  # ITS + LSU
  ENABLE_16S=0 ENABLE_ITS=1 ENABLE_LSU=1

The pipeline automatically routes FASTQs to the appropriate marker
based on filename patterns:

  16SA / 16SB → 16S
  ITS          → ITS
  LSU          → LSU

Libraries without a matching marker pattern are skipped for that marker.

If a marker is enabled but its Emu database is missing, the marker is
disabled automatically with a warning message.

NOTE:
These variables affect the Emu stage only.
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
-----------
Emu runs independently for each enabled marker and processes only
FASTQs whose filenames contain the corresponding marker label.

Markers:
  16SA / 16SB → 16S
  ITS          → ITS
  LSU          → LSU

Resources:
  --emu-partition STR   Partition for Emu (default: inherit libsQC)
  --emu-time HH:MM:SS   Walltime for Emu
  --emu-cpus INT        CPUs for Emu
  --emu-mem STR         Memory for Emu

Databases:
  --emu-db-its PATH     Custom ITS Emu database directory
  --emu-db-lsu PATH     Custom LSU Emu database directory

Database behavior:
  - 16S DB is automatically downloaded at first use
    (bacteria + archaea reference database)

  - ITS and LSU databases are prepared through:

      workflow/run_build_ITS_LSU_dbs.sh

  - Missing ITS/LSU databases result in marker skipping
    with warning messages.

Input FASTQs:
By default, Emu uses:

  results/[batch]/filtered/

This can be overridden using:

  FASTQ_DIR_DEFAULT=/path/to/filtered

Outputs are marker-specific and batch-specific:

  results/emu_runs_16S_[batch]/
  results/emu_runs_ITS_[batch]/
  results/emu_runs_LSU_[batch]/

  results/tables_16S_[batch]/
  results/tables_ITS_[batch]/
  results/tables_LSU_[batch]/

  results/plots_16S_[batch]/
  results/plots_ITS_[batch]/
  results/plots_LSU_[batch]/
</pre>

<pre>
EMU OUTPUTS AND COLLATION
-------------------------
For each sample, Emu outputs are stored in:

  results/emu_runs_<marker>_[batch]/<sample>/

including:

  abundance.tsv
  estimated counts
  read assignments
  input_reads.tsv

After classification, `workflow/emu_collect.py` collates results into
marker-level summary tables.

Generated files:

  abundance_combined.tsv
  mapping_stats.tsv

abundance_combined.tsv
----------------------
A merged table preserving the original Emu abundance columns
(e.g., tax_id, abundance, estimated_counts, species, genus, etc.)
with an additional column:

  file

indicating the originating library.

mapping_stats.tsv
-----------------
Mapping statistics are computed using fractional assignment
probabilities from Emu read-assignment distributions.

For each sample, the table reports:

  total_reads
  assigned_reads
  assigned_frac
  unassigned_reads
  unassigned_frac

Assignment is computed fractionally from posterior probabilities rather
than hard thresholding.

NOTE:
`--min-prob` is retained for backward compatibility but is ignored in
the current implementation.
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

# 3) Run ONLY Emu (QC already done)
#    Run all enabled markers on filtered FASTQs
#    Example: batch 250–274 (25 FASTQs)

BATCH_TAG=b250_n025_filter \
ENABLE_16S=1 ENABLE_ITS=1 ENABLE_LSU=1 \
bash workflow/runall.sh \
  --no-qc \
  --no-downstream \
  --no-build-marker-dbs \
  --emu-cpus 20 \
  --emu-mem 64G \
  --emu-time 12:00:00

# 4) Run ONLY ITS + LSU (explicitly skipping 16S)

BATCH_TAG=b250_n025_filter \
ENABLE_16S=0 ENABLE_ITS=1 ENABLE_LSU=1 \
bash workflow/runall.sh \
  --no-qc \
  --no-downstream \
  --no-build-marker-dbs \
  --emu-cpus 20 \
  --emu-mem 64G \
  --emu-time 12:00:00

# 16) Run ONLY the downstream analysis (16S) (FIX FROM HERE!!)
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
 logs/[batch]_*                         - batch-specific job logs
 results/[batch]/trimmed/              - primer-trimmed FASTQs
 results/[batch]/untrimmed/            - reads without detected primers
 results/[batch]/filtered/             - post-filtered FASTQs
 results/[batch]/trim_reports/         - Cutadapt reports per library
 results/[batch]/qc_raw/               - FastQC reports for primer-trimmed reads
 results/[batch]/qc_filtered/          - FastQC reports for NanoFilt-filtered reads
 results/[batch]/multiqc/              - MultiQC reports
 results/[batch]/nanoplot/             - NanoPlot outputs
 results/[batch]/summary/              - QC summaries and flags
   results/[batch]/summary/primer_trimming_by_library.tsv
   results/[batch]/summary/nanofilt_cutoffs_by_library.tsv
   results/[batch]/summary/qc_flags.tsv
   results/[batch]/summary/seqkit_stats.tsv
 results/[batch]/lengths/              - read-length distributions
 results/emu_runs_16S_[batch]/         - per-sample Emu outputs (16S)
 results/emu_runs_ITS_[batch]/         - per-sample Emu outputs (ITS)
 results/emu_runs_LSU_[batch]/         - per-sample Emu outputs (LSU)

 results/tables_16S_[batch]/           - merged Emu abundance + mapping tables
 results/tables_ITS_[batch]/           - merged Emu abundance + mapping tables
 results/tables_LSU_[batch]/           - merged Emu abundance + mapping tables

 results/plots_16S_[batch]/            - genus stacked-bar plots
 results/plots_ITS_[batch]/            - genus stacked-bar plots
 results/plots_LSU_[batch]/            - genus stacked-bar plots

 metadata/                             - FASTQ manifests and primer lists
   metadata/fastq_meta.[batch].tsv       - sample manifest
   metadata/primers_fwd.list             - forward primers used
   metadata/primers_rev.list             - reverse primers used
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


