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
        "r-base>=4.5" \
        r-data.table r-ggplot2 r-vegan r-tidyr r-dplyr r-stringr \
        r-cowplot r-rcolorbrewer r-ggbeeswarm \
        r-biocmanager \
        -y
    else
      conda create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        "r-base>=4.5" \
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
# ### CHANGED ### Install core Bioconductor deps via conda/mamba
# ----------------------------------------------------------
install_ancombc_deps_conda() {
  log "Ensuring core Bioconductor dependencies (DelayedArray/GenomicRanges/etc) are installed via conda..."

  local pkgs=(
    bioconductor-delayedarray
    bioconductor-genomicranges
    bioconductor-summarizedexperiment
    bioconductor-sparsearray
    bioconductor-delayedmatrixstats
    bioconductor-biostrings
    bioconductor-iranges
    bioconductor-s4vectors
    bioconductor-xvector
    bioconductor-genomeinfodb
    bioconductor-genomeinfodbdata
    r-matrix
  )

  if command -v mamba >/dev/null 2>&1; then
    log "mamba install -n ${ENV_NAME} -c conda-forge -c bioconda ${pkgs[*]}"
    mamba install -n "${ENV_NAME}" -c conda-forge -c bioconda "${pkgs[@]}" -y
  else
    log "conda install -n ${ENV_NAME} -c conda-forge -c bioconda ${pkgs[*]}"
    conda install -n "${ENV_NAME}" -c conda-forge -c bioconda "${pkgs[@]}" -y
  fi

  log "Core Bioconductor deps installed (conda)."
}

# ----------------------------------------------------------
# ensure ANCOMBC is installed (GitHub preferred)
# ----------------------------------------------------------
ensure_ancombc() {
  log "Checking for R package 'ANCOMBC' (ANCOM-BC2)..."

  classify="UNKNOWN"

  # Robust classification (never crash the bash script)
  classify="$(
    Rscript -e '
      out <- "UNKNOWN"
      tryCatch({
        if (!requireNamespace("ANCOMBC", quietly=TRUE)) {
          out <- "NOT_INSTALLED"
        } else {
          d <- utils::packageDescription("ANCOMBC")

          repo <- ""
          if (!is.null(d[["Repository"]])) repo <- as.character(d[["Repository"]])

          has_bioc_views <- FALSE
          if (!is.null(d[["biocViews"]])) {
            bv <- as.character(d[["biocViews"]])
            has_bioc_views <- nchar(bv) > 0
          }

          remote_fields <- c("RemoteType","RemoteRepo","RemoteUsername","RemoteRef","RemoteSha")
          has_remote <- any(sapply(remote_fields, function(x) {
            !is.null(d[[x]]) && nchar(as.character(d[[x]])) > 0
          }))

          if (has_remote) {
            out <- "GITHUB"
          } else if (grepl("Bioconductor", repo, ignore.case=TRUE) || has_bioc_views) {
            out <- "BIOCONDUCTOR"
          } else {
            out <- "OTHER"
          }
        }
      }, error=function(e) {
        out <- "UNKNOWN"
      })
      cat(out)
    ' 2>/dev/null || echo "UNKNOWN"
  )"

  log "ANCOMBC status detected: ${classify}"

  if [[ "${classify}" == "GITHUB" ]]; then
    log "ANCOMBC is already from GitHub. Keeping it (no action)."
    return 0
  fi

  if [[ "${classify}" == "BIOCONDUCTOR" ]]; then
    log "ANCOMBC appears to be from Bioconductor. Removing it..."
    Rscript -e 'suppressWarnings(try(remove.packages("ANCOMBC"), silent=TRUE))' || true
    log "ANCOMBC removed."
  elif [[ "${classify}" == "NOT_INSTALLED" ]]; then
    log "ANCOMBC not installed."
  elif [[ "${classify}" == "OTHER" ]]; then
    log "ANCOMBC installed but origin is OTHER/unknown. Per your rule: only replace if Bioconductor. Keeping it."
    return 0
  else
    log "ANCOMBC classification UNKNOWN. Per your rule: only replace if Bioconductor. Keeping it."
    return 0
  fi

  # Ensure Bioconductor dependency chain is satisfied via conda first
  install_ancombc_deps_conda

  # Install from GitHub WITHOUT pulling huge dependency trees
  log "Ensuring R package 'remotes' is available..."
  Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org")'

  log 'Installing ANCOMBC from GitHub (dependencies=FALSE): FrederickHuangLin/ANCOMBC'
  Rscript -e 'remotes::install_github("FrederickHuangLin/ANCOMBC", upgrade="never", dependencies=FALSE)'

  log "Verifying ANCOMBC after GitHub install..."
  if ! Rscript -e 'quit(status = ifelse(requireNamespace("ANCOMBC", quietly=TRUE), 0, 1))'; then
    die "Failed to install ANCOMBC from GitHub."
  fi

  log "ANCOMBC installed successfully from GitHub."
  Rscript -e 'cat("ANCOMBC version: ", as.character(packageVersion("ANCOMBC")), "\n", sep="")' || true
  Rscript -e 'cat("ANCOMBC path: ", find.package("ANCOMBC"), "\n", sep="")' || true
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
