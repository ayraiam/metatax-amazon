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

EXAMPLES
---------
# 1) Run both stages (QC + Emu, default settings)
bash workflow/runall.sh

# 2) Run both with custom resources
bash workflow/runall.sh --time 08:00:00 --cpus 8 --mem 32G

# 3) Run only QC
bash workflow/runall.sh --no-emu --time 08:00:00 --cpus 8 --mem 32G

# 4) Run only Emu (assuming QC already done)
bash workflow/runall.sh --no-qc --emu-time 12:00:00 --emu-cpus 16 --emu-mem 64G

# 5) Run both, giving Emu more resources
bash workflow/runall.sh \
  --time 06:00:00 --cpus 8 --mem 32G \
  --emu-time 12:00:00 --emu-cpus 16 --emu-mem 64G

# 6) Run Emu with ITS DB
bash workflow/runall.sh --no-qc --emu-db-its /path/to/emu_unite_its

# 7) Run Emu with LSU DB
bash workflow/runall.sh --no-qc --emu-db-lsu /path/to/emu_silva_lsu
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
