#!/usr/bin/env bash
# ==========================================================
# Script: workflow/run_downstream_analysis.sh
# Purpose:
#   - Create/activate a clean env for downstream analysis
#   - Ensure ANCOMBC (for ANCOM-BC2) is available (GitHub preferred)
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

# ==========================================================
# Helper: run mamba/conda consistently
# ==========================================================
conda_install() {
  # usage: conda_install <pkgs...>
  if command -v mamba >/dev/null 2>&1; then
    log "mamba install -n ${ENV_NAME} -c conda-forge -c bioconda $*"
    mamba install -n "${ENV_NAME}" -c conda-forge -c bioconda "$@" -y
  else
    log "conda install -n ${ENV_NAME} -c conda-forge -c bioconda $*"
    conda install -n "${ENV_NAME}" -c conda-forge -c bioconda "$@" -y
  fi
}

# ==========================================================
# Robust R>=4.5 assertion + auto-fix (prevents silent R=4.4.3)
# ==========================================================
assert_and_fix_r_version() {
  log "Checking R version inside env..."

  # print which Rscript is being used
  log "Rscript path: $(command -v Rscript || echo 'NOT_FOUND')"

  local rv
  rv="$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || echo "0.0.0")"
  log "R version in env (via Rscript): ${rv}"

  # also print R.home() to confirm it’s from the env
  Rscript -e 'cat("R.home(): ", R.home(), "\n", sep="")' || true

  # If R < 4.5, force-upgrade r-base inside the env (don’t just die)
  if ! Rscript -e 'quit(status=ifelse(getRversion() >= "4.5.0", 0, 1))'; then
    log "R < 4.5.0 detected. Forcing r-base=4.5.* upgrade inside env..."
    conda_install "r-base=4.5.*"

    # re-check after upgrade
    rv="$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || echo "0.0.0")"
    log "R version after forced install: ${rv}"

    if ! Rscript -e 'quit(status=ifelse(getRversion() >= "4.5.0", 0, 1))'; then
      die "Still R < 4.5.0 after forced r-base=4.5.* install. Something is pinning/downgrading the env."
    fi
  fi

  log "R version OK (>= 4.5.0)."
}

# ==========================================================
# install build toolchain + system libs needed to compile Bioconductor from source
# ==========================================================
install_build_deps_conda() {
  log "Ensuring build/system dependencies exist in the env (for compiling Bioconductor packages from source)..."

  local pkgs=(
    make pkg-config cmake autoconf automake libtool
    gcc gxx gfortran llvm-openmp
    zlib bzip2 xz libcurl openssl libxml2 pcre2 icu readline
  )

  # don’t randomly install phyloseq here (unrelated + can affect solver)
  conda_install "${pkgs[@]}"

  log "Build/system deps installed."
}

# ==========================================================
# install heavy CRAN deps via conda to avoid compilation failures (nloptr/lme4/CVXR etc)
# ==========================================================
install_ancombc_cran_deps_conda() {
  log "Installing ANCOMBC CRAN dependencies via conda (reduces compilation problems)..."

  local pkgs=(
    r-cvxr r-desctools r-hmisc r-rdpack r-doparallel r-dorng
    r-energy r-foreach r-gtools r-lme4 r-lmertest r-multcomp r-nloptr
  )

  # don’t install phyloseq as a side-effect
  conda_install "${pkgs[@]}"

  log "CRAN deps (conda) installed."
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
        "r-base=4.5.*" \
        r-data.table r-ggplot2 r-vegan r-tidyr r-dplyr r-stringr \
        r-cowplot r-rcolorbrewer r-ggbeeswarm \
        r-biocmanager r-remotes \
        -y
    else
      conda create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        "r-base=4.5.*" \
        r-data.table r-ggplot2 r-vegan r-tidyr r-dplyr r-stringr \
        r-cowplot r-rcolorbrewer r-ggbeeswarm \
        r-biocmanager r-remotes \
        -y
    fi
    conda activate "${ENV_NAME}"
  fi

  command -v Rscript >/dev/null 2>&1 || die "Rscript not available in env."
  log "Env ready: $(which Rscript)"

  # auto-fix, not just assert
  assert_and_fix_r_version

  # Ensure toolchain exists BEFORE BiocManager tries compiling stuff like XVector
  install_build_deps_conda

  install_phyloseq_conda
}

install_phyloseq_conda() {
  log "Installing phyloseq via conda (bioconductor-phyloseq)..."
  conda_install "bioconductor-phyloseq"
  log "phyloseq installed."
}

# ----------------------------------------------------------
# Install Bioconductor deps via BiocManager (inside R)
# ----------------------------------------------------------
install_bioc_deps_r() {
  log "Installing required Bioconductor dependencies via BiocManager (R-side, matches your R=4.5.*)..."

  Rscript -e '
    options(Ncpus = 1)
    Sys.setenv(MAKEFLAGS = "-j1")

    if (!requireNamespace("BiocManager", quietly=TRUE)) {
      install.packages("BiocManager", repos="https://cloud.r-project.org")
    }

    cat(">>> R: ", as.character(getRversion()), "\n", sep="")
    cat(">>> Bioconductor version: ", as.character(BiocManager::version()), "\n", sep="")
    cat(">>> .libPaths():\n"); print(.libPaths())

    pkgs <- c(
      "XVector",
      "SummarizedExperiment","GenomicRanges","DelayedArray","SparseArray",
      "DelayedMatrixStats","Biostrings","IRanges","S4Vectors",
      "GenomeInfoDb","GenomeInfoDbData"
    )

    BiocManager::install("XVector", update=FALSE, ask=FALSE)
    BiocManager::install(setdiff(pkgs, "XVector"), update=FALSE, ask=FALSE)
  '
  log "Bioconductor dependency install completed (R-side)."
}

# ----------------------------------------------------------
# ensure ANCOMBC is installed (GitHub preferred)
# ----------------------------------------------------------
ensure_ancombc() {
  log "Checking for R package 'ANCOMBC' (ANCOM-BC2)..."

  classify="UNKNOWN"
  classify="$(
    Rscript -e '
      out <- "UNKNOWN"
      tryCatch({
        if (!requireNamespace("ANCOMBC", quietly=TRUE)) {
          out <- "NOT_INSTALLED"
        } else {
          d <- utils::packageDescription("ANCOMBC")
          repo <- if (!is.null(d[["Repository"]])) as.character(d[["Repository"]]) else ""
          has_bioc_views <- !is.null(d[["biocViews"]]) && nchar(as.character(d[["biocViews"]])) > 0
          remote_fields <- c("RemoteType","RemoteRepo","RemoteUsername","RemoteRef","RemoteSha")
          has_remote <- any(sapply(remote_fields, function(x) !is.null(d[[x]]) && nchar(as.character(d[[x]])) > 0))
          if (has_remote) out <- "GITHUB"
          else if (grepl("Bioconductor", repo, ignore.case=TRUE) || has_bioc_views) out <- "BIOCONDUCTOR"
          else out <- "OTHER"
        }
      }, error=function(e) out <- "UNKNOWN")
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
  else
    log "ANCOMBC classification ${classify}. Per your rule: only replace if Bioconductor/Not installed. Keeping it."
    return 0
  fi

  install_ancombc_cran_deps_conda
  install_bioc_deps_r

  log "Installing ANCOMBC from GitHub (dependencies=FALSE, build=FALSE)..."
  Rscript -e 'remotes::install_github(
    "FrederickHuangLin/ANCOMBC",
    upgrade = "never",
    dependencies = FALSE,
    build = FALSE,
    build_vignettes = FALSE
  )'

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
