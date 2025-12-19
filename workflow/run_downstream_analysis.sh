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
        r-biocmanager \
        -y
    else
      conda create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        r-base=4.3 \
        r-data.table r-ggplot2 r-vegan r-tidyr r-dplyr r-stringr \
        r-cowplot r-rcolorbrewer r-ggbeeswarm \
        r-biocmanager \
        -y
    fi
    conda activate "${ENV_NAME}"
  fi

  command -v Rscript >/dev/null 2>&1 || die "Rscript not available in env."
  log "Env ready: $(which Rscript)"
}

# ----------------------------------------------------------
# ensure ANCOMBC is installed
# ----------------------------------------------------------
ensure_ancombc() {
  log "Checking for R package 'ANCOMBC' (ANCOM-BC2)..."

  # We do a 3-way classification:
  #  - NOT_INSTALLED
  #  - GITHUB
  #  - BIOCONDUCTOR
  #  - OTHER (unknown origin; we'll leave it alone unless you want to force GitHub)
  #
  # Heuristic:
  #  - GitHub installs typically include RemoteType/RemoteRepo fields in package DESCRIPTION.
  #  - Bioconductor installs typically have Repository: Bioconductor and/or biocViews field.

  classify="$(
    Rscript -e '
      status <- "NOT_INSTALLED"
      if (requireNamespace("ANCOMBC", quietly=TRUE)) {
        d <- utils::packageDescription("ANCOMBC")
        repo <- if (!is.null(d[["Repository"]])) as.character(d[["Repository"]]) else ""
        has_bioc_views <- !is.null(d[["biocViews"]]) && nchar(as.character(d[["biocViews"]])) > 0

        has_remote <- any(!is.null(d[[x]]) && nchar(as.character(d[[x]])) > 0
                          for (x in c("RemoteType","RemoteRepo","RemoteUsername","RemoteRef","RemoteSha")))

        # GitHub/remotes/devtools install signature
        if (has_remote) {
          status <- "GITHUB"
        } else if (grepl("Bioconductor", repo, ignore.case=TRUE) || has_bioc_views) {
          status <- "BIOCONDUCTOR"
        } else {
          status <- "OTHER"
        }
      }
      cat(status)
    ' 2>/dev/null
  )"

  log "ANCOMBC status detected: ${classify}"

  if [[ "${classify}" == "GITHUB" ]]; then
    log "ANCOMBC is already from GitHub. Keeping it (no action)."
    # Optional: print version
    Rscript -e 'cat("ANCOMBC version: ", as.character(packageVersion("ANCOMBC")), "\n", sep="")' || true
    return 0
  fi

  if [[ "${classify}" == "BIOCONDUCTOR" ]]; then
    log "ANCOMBC appears to be from Bioconductor. Removing it..."
    Rscript -e 'suppressWarnings(try(remove.packages("ANCOMBC"), silent=TRUE))'
    log "ANCOMBC removed. Installing from GitHub..."
  elif [[ "${classify}" == "NOT_INSTALLED" ]]; then
    log "ANCOMBC not installed. Installing from GitHub..."
  else
    # OTHER
    log "ANCOMBC installed but origin is unclear (OTHER)."
    log "Per your spec: only replace if Bioconductor. Keeping it (no action)."
    # Optional: print version/path
    Rscript -e 'cat("ANCOMBC version: ", as.character(packageVersion("ANCOMBC")), "\n", sep="")' || true
    Rscript -e 'cat("ANCOMBC path: ", find.package("ANCOMBC"), "\n", sep="")' || true
    return 0
  fi

  # Install from GitHub (lightweight installer: remotes)
  # We ensure remotes exists, then install_github.

  log "Ensuring R package 'remotes' is available..."
  Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org")'

  log 'Installing ANCOMBC from GitHub: FrederickHuangLin/ANCOMBC'
  Rscript -e 'remotes::install_github("FrederickHuangLin/ANCOMBC", upgrade="never", dependencies=TRUE)'

  log "Verifying ANCOMBC after GitHub install..."
  if ! Rscript -e 'quit(status = ifelse(requireNamespace("ANCOMBC", quietly=TRUE), 0, 1))'; then
    die "Failed to install ANCOMBC from GitHub."
  fi

  # Re-classify to confirm it’s GitHub now
  classify2="$(
    Rscript -e '
      status <- "NOT_INSTALLED"
      if (requireNamespace("ANCOMBC", quietly=TRUE)) {
        d <- utils::packageDescription("ANCOMBC")
        has_remote <- any(!is.null(d[[x]]) && nchar(as.character(d[[x]])) > 0
                          for (x in c("RemoteType","RemoteRepo","RemoteUsername","RemoteRef","RemoteSha")))
        status <- if (has_remote) "GITHUB" else "NON_GITHUB"
      }
      cat(status)
    ' 2>/dev/null
  )"
  log "Post-install ANCOMBC status: ${classify2}"
  Rscript -e 'cat("ANCOMBC version: ", as.character(packageVersion("ANCOMBC")), "\n", sep="")' || true
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
