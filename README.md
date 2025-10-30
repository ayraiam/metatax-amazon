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
STRUCTURE [TO BE UPDATED]
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
DEPENDENCIES [TO BE UPDATED]
------------
 - Python ≥ 3.10
 - R ≥ 4.3
 - Conda or Mamba
 - Slurm (or compatible scheduler)
 - NanoPlot · FastQC · MultiQC · SeqKit · Cutadapt · NanoFilt
</pre>

<pre>
USAGE [TO BE UPDATED]
-----
Example (cluster execution):

  bash workflow/runall.sh \
  --partition short \
  --time 04:00:00 \
  --cpus 4 \
  --mem 16G \
  --wd "$PWD"
# optional flags:
#   --primer-fwd SEQ --primer-rev SEQ
#   --seq-summary /path/to/sequencing_summary.txt
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
CONTACT [TO BE UPDATED]
-------
email: ayrabioinf@gmail.com
linkedin: https://www.linkedin.com/company/aryaiam
</pre>

<p align="center"><sub>© 2025 AY:RΔ — data and discovery in flow</sub></p>
