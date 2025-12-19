#!/usr/bin/env bash
# ==========================================================
# Script: workflow/run_downstream_analysis.sh
# Purpose:
#   - Create/activate a clean env for downstream analysis
#   - Ensure ANCOMBC (for ANCOM-BC2) is available
#   - Run downstream_analysis.R
# ==========================================================
set -euo pipefail

# ---------------- User-tunable defaults -------------------
ENV_NAME="${ENV_NAME:-emu-downstream}"
INFILE="${INFILE:-/home/t.sousa/metataxonomy_rds/metatax-amazon/results/tables/abundance_combined.tsv}"
OUTDIR="${OUTDIR:-/home/t.sousa/metataxonomy_rds/metatax-amazon/results/plots}"
BASENAME="${BASENAME:-downstream}"

LOGDIR="${LOGDIR:-logs}"
mkdir -p "$LOGDIR" "$(dirname "$OUTDIR")" "$OUTDIR"

log()  { echo -e ">>> $*" >&2; }
die()  { echo -e "xxx $*" >&2; exit 1; }

ensure_channels() {
  conda config --remove-key channels >/dev/null 2>&1 || true
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda config --add channels defaults
  conda config --set channel_priority strict
}

create_env() {
  if ! command -v conda >/dev/null 2>&1; then
    die "conda not found. Install Miniforge/Conda first."
  fi
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"

  if conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    log "Environment '${ENV_NAME}' exists. Activating..."
    conda activate "${ENV_NAME}"
  else
    log "Creating environment '${ENV_NAME}'..."
    ensure_channels
    if command -v mamba >/dev/null 2>&1; then
      mamba create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        r-base=4.3 \
        r-data.table r-ggplot2 r-vegan r-tidyr r-dplyr r-stringr \
        r-cowplot r-rcolorbrewer r-ggbeeswarm \
        r-devtools \
        -y
    else
      conda create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        r-base=4.3 \
        r-data.table r-ggplot2 r-vegan r-tidyr r-dplyr r-stringr \
        r-cowplot r-rcolorbrewer r-ggbeeswarm \
        r-devtools \
        -y
    fi
    conda activate "${ENV_NAME}"
  fi

  command -v Rscript >/dev/null 2>&1 || die "Rscript not available in env."
  log "Env ready: $(which Rscript)"
}

# ----------------------------------------------------------
# ensure ANCOMBC is installed (GitHub only)
# ----------------------------------------------------------
ensure_ancombc() {
  log "Checking for R package 'ANCOMBC' (ANCOM-BC2)..."

  if Rscript -e 'quit(status = ifelse(requireNamespace("ANCOMBC", quietly=TRUE), 0, 1))'; then
    log "ANCOMBC already installed."
    return 0
  fi

  log "ANCOMBC not found. Will install from GitHub."

  # --- best-effort cleanup of conda/Bioc-installed versions ---
  # (helps avoid mixed installs / stale versions)
  if command -v mamba >/dev/null 2>&1; then
    log "Removing conda package if present: bioconductor-ancombc"
    mamba remove -n "${ENV_NAME}" -y bioconductor-ancombc >/dev/null 2>&1 || true
  else
    log "Removing conda package if present: bioconductor-ancombc"
    conda remove -n "${ENV_NAME}" -y bioconductor-ancombc >/dev/null 2>&1 || true
  fi

  # Also try removing any R-library copy (best effort; don't fail if missing)
  Rscript -e 'try(remove.packages("ANCOMBC"), silent=TRUE); try(remove.packages("ancombc"), silent=TRUE)' || true

  # --- install from GitHub using devtools ---
  log 'Installing ANCOMBC from GitHub: FrederickHuangLin/ANCOMBC'
  Rscript -e '
    if (!requireNamespace("devtools", quietly=TRUE)) {
      install.packages("devtools", repos="https://cloud.r-project.org")
    }
    devtools::install_github("FrederickHuangLin/ANCOMBC", upgrade="never", dependencies=TRUE)
  '

  # Final check
  if ! Rscript -e 'quit(status = ifelse(requireNamespace("ANCOMBC", quietly=TRUE), 0, 1))'; then
    die "Failed to install ANCOMBC from GitHub. Check build deps (compilers/headers) and network."
  fi

  log "ANCOMBC installed via GitHub."
}

run_downstream() {

  export USE_COUNTS_0_4="${USE_COUNTS_0_4:-0}"
  export USE_COUNTS_5="${USE_COUNTS_5:-1}"
  export MODE="${MODE:-16S}"

  log "Running downstream_analysis.R with:"
  log "  infile        = $INFILE"
  log "  outdir        = $OUTDIR"
  log "  basename      = $BASENAME"
  log "  MODE          = $MODE"
  log "  USE_COUNTS_0_4= $USE_COUNTS_0_4"
  log "  USE_COUNTS_5  = $USE_COUNTS_5"

  Rscript workflow/downstream_analysis.R \
    "$INFILE" \
    "$OUTDIR" \
    "$BASENAME"

  log "Plots + tables written to: $OUTDIR"
}

create_env
ensure_ancombc
run_downstream
