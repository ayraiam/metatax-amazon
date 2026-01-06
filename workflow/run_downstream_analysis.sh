#!/usr/bin/env bash
# ==========================================================
# Script: workflow/run_downstream_analysis.sh
# Purpose:
#   - Create/activate a clean env for downstream analysis
#   - Ensure phyloseq + microbiome are available (R-side; no conda downgrade)
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
# Strong check: confirm we're using ENV Rscript + R>=4.5
# ==========================================================
assert_r_version_and_path() {
  log "Checking R version + making sure we're using the ENV's Rscript..."

  local expected actual rv
  expected="$(conda info --base)/envs/${ENV_NAME}/bin/Rscript"
  actual="$(command -v Rscript || true)"

  log "Expected env Rscript: ${expected}"
  log "Actual   Rscript on PATH: ${actual}"

  if [[ -z "${actual}" ]] || [[ "${actual}" != "${expected}" ]]; then
    die "Rscript on PATH is not the env Rscript. Something is wrong with conda activation."
  fi

  rv="$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || echo "0.0.0")"
  log "R version in env: ${rv}"

  log "Conda r-base entry:"
  (conda list -n "${ENV_NAME}" r-base || true) >&2

  if ! Rscript -e 'quit(status=ifelse(getRversion() >= "4.5.0", 0, 1))'; then
    die "R < 4.5.0 detected in env. Delete env and recreate."
  fi

  log "R version OK (>= 4.5.0) and PATH is correct."
}

# ==========================================================
# Install toolchain + system libs in env
# ==========================================================
install_build_deps_conda() {
  log "Ensuring build/system dependencies exist in the env (for compiling/loading Bioconductor packages)..."

  local pkgs=(
    make
    pkg-config
    cmake
    autoconf
    automake
    libtool
    gcc
    gxx
    gfortran

    zlib
    bzip2
    xz
    libcurl
    openssl
    libxml2
    pcre2
    icu
    readline

    # HDF5 often needed for rhdf5/biomformat/phyloseq loadability
    hdf5
  )

  if command -v mamba >/dev/null 2>&1; then
    log "mamba install -n ${ENV_NAME} -c conda-forge ${pkgs[*]}"
    mamba install -n "${ENV_NAME}" -c conda-forge "${pkgs[@]}" -y
  else
    log "conda install -n ${ENV_NAME} -c conda-forge ${pkgs[*]}"
    conda install -n "${ENV_NAME}" -c conda-forge "${pkgs[@]}" -y
  fi

  log "Build/system deps installed."
}

# ==========================================================
# install heavy CRAN deps via conda to avoid compilation failures
# ==========================================================
install_ancombc_cran_deps_conda() {
  log "Installing ANCOMBC CRAN dependencies via conda (reduces compilation problems)..."

  local pkgs=(
    r-cvxr
    r-desctools
    r-hmisc
    r-rdpack
    r-doparallel
    r-dorng
    r-energy
    r-foreach
    r-gtools
    r-lme4
    r-lmertest
    r-multcomp
    r-nloptr
  )

  if command -v mamba >/dev/null 2>&1; then
    log "mamba install -n ${ENV_NAME} -c conda-forge ${pkgs[*]}"
    mamba install -n "${ENV_NAME}" -c conda-forge "${pkgs[@]}" -y
  else
    log "conda install -n ${ENV_NAME} -c conda-forge ${pkgs[*]}"
    conda install -n "${ENV_NAME}" -c conda-forge "${pkgs[@]}" -y
  fi

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
        r-biocmanager \
        r-remotes \
        -y
    else
      conda create -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        "r-base=4.5.*" \
        r-data.table r-ggplot2 r-vegan r-tidyr r-dplyr r-stringr \
        r-cowplot r-rcolorbrewer r-ggbeeswarm \
        r-biocmanager \
        r-remotes \
        -y
    fi
    conda activate "${ENV_NAME}"
  fi

  command -v Rscript >/dev/null 2>&1 || die "Rscript not available in env."
  log "Env ready: $(which Rscript)"

  assert_r_version_and_path
  install_build_deps_conda
}

# ----------------------------------------------------------
# Install Bioconductor deps via BiocManager (inside R)
# ----------------------------------------------------------
install_bioc_deps_r() {
  log "Installing required Bioconductor dependencies via BiocManager (R-side; matches your R=4.5.*)..."

  Rscript -e '
    options(Ncpus = 1)
    Sys.setenv(MAKEFLAGS = "-j1")

    if (!requireNamespace("BiocManager", quietly=TRUE)) {
      install.packages("BiocManager", repos="https://cloud.r-project.org")
    }

    pkgs <- c(
      "XVector",
      "IRanges","S4Vectors",
      "GenomeInfoDb","GenomeInfoDbData",
      "GenomicRanges",
      "SparseArray",
      "DelayedArray",
      "SummarizedExperiment",
      "DelayedMatrixStats",
      "Biostrings"
    )

    BiocManager::install(pkgs, update=FALSE, ask=FALSE)
  ' || die "Bioconductor dependency install failed."

  log "Bioconductor dependency install completed (R-side)."
}

# ==========================================================
# Install phyloseq + microbiome via BiocManager and verify load
# ==========================================================
install_phyloseq_and_microbiome_r() {
  log "Ensuring 'phyloseq' + 'microbiome' are installed and loadable (BiocManager; no conda bioconductor-phyloseq)..."
  assert_r_version_and_path

  Rscript -e '
    options(Ncpus = 1)
    Sys.setenv(MAKEFLAGS = "-j1")

    if (!requireNamespace("BiocManager", quietly=TRUE)) {
      install.packages("BiocManager", repos="https://cloud.r-project.org")
    }

    # Install biomformat + rhdf5 stack explicitly (common failure point)
    BiocManager::install(c("rhdf5lib","rhdf5","biomformat"), update=FALSE, ask=FALSE)

    # Install phyloseq
    BiocManager::install("phyloseq", update=FALSE, ask=FALSE)

    ### CHANGED ### Install microbiome (required by ANCOMBC when data is phyloseq)
    BiocManager::install("microbiome", update=FALSE, ask=FALSE)

    # Verify loadability
    ok <- TRUE
    tryCatch({ suppressPackageStartupMessages(library(phyloseq)) }, error=function(e){
      ok <<- FALSE
      cat(">>> ERROR loading phyloseq:\n"); cat(conditionMessage(e), "\n")
    })
    tryCatch({ suppressPackageStartupMessages(library(microbiome)) }, error=function(e){
      ok <<- FALSE
      cat(">>> ERROR loading microbiome:\n"); cat(conditionMessage(e), "\n")
    })

    quit(status=ifelse(ok, 0, 1))
  ' || die "phyloseq/microbiome installation failed (package still not loadable)."

  log "phyloseq + microbiome are installed and loadable."
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
  elif [[ "${classify}" == "OTHER" ]]; then
    log "ANCOMBC installed but origin is OTHER/unknown. Per your rule: only replace if Bioconductor. Keeping it."
    return 0
  else
    log "ANCOMBC classification UNKNOWN. Per your rule: only replace if Bioconductor. Keeping it."
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

# ==========================================================
# Ensure rotl + ape (conda preferred; CRAN fallback)
# ==========================================================
ensure_rotl_and_ape() {
  log "Ensuring 'rotl' + 'ape' are available in env '${ENV_NAME}'..."

  # --- First try conda-forge packages (fast, binary) ---
  local missing_conda=()

  conda list -n "${ENV_NAME}" r-rotl >/dev/null 2>&1 || missing_conda+=("r-rotl")
  conda list -n "${ENV_NAME}" r-ape  >/dev/null 2>&1 || missing_conda+=("r-ape")

  if [[ ${#missing_conda[@]} -gt 0 ]]; then
    log "Conda missing: ${missing_conda[*]} -> installing via conda-forge..."
    if command -v mamba >/dev/null 2>&1; then
      mamba install -n "${ENV_NAME}" -c conda-forge "${missing_conda[@]}" -y
    else
      conda install -n "${ENV_NAME}" -c conda-forge "${missing_conda[@]}" -y
    fi
  else
    log "Conda shows r-rotl and r-ape already installed."
  fi

  # --- Verify from THIS Rscript (critical) ---
  log "Verifying rotl/ape loadability from R (this must succeed)..."
  if ! Rscript -e 'suppressPackageStartupMessages(library(rotl)); suppressPackageStartupMessages(library(ape)); cat("OK\n")'; then
    log "rotl/ape still not loadable from R. Falling back to CRAN install inside env..."

    Rscript -e '
      options(repos = c(CRAN="https://cloud.r-project.org"))
      pkgs <- c("ape","rotl")
      for (p in pkgs) {
        if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
      }
      suppressPackageStartupMessages(library(ape))
      suppressPackageStartupMessages(library(rotl))
      cat("OK\n")
    ' || die "Failed to make rotl/ape loadable (conda + CRAN fallback both failed)."
  fi

  log "rotl + ape are installed and loadable."
}

run_downstream() {
  export USE_COUNTS_0_4="${USE_COUNTS_0_4:-0}"
  export USE_COUNTS_5="${USE_COUNTS_5:-1}"
  export MODE="${MODE:-16S}"

  #Dedicated ANCOM selector; default to estimated_counts
  export USE_COUNTS_ANCOM="${USE_COUNTS_ANCOM:-1}"

  log "Running downstream_analysis.R with:"
  log "  infile         = $INFILE"
  log "  outdir         = $OUTDIR"
  log "  basename       = $BASENAME"
  log "  MODE           = $MODE"
  log "  USE_COUNTS_0_4 = $USE_COUNTS_0_4"
  log "  USE_COUNTS_5   = $USE_COUNTS_5"
  log "  USE_COUNTS_ANCOM = $USE_COUNTS_ANCOM"

  Rscript workflow/downstream_analysis.R \
    "$INFILE" \
    "$OUTDIR" \
    "$BASENAME"

  log "Plots + tables written to: $OUTDIR"
}

create_env

# phyloseq + microbiome must be installed BEFORE downstream step7
install_phyloseq_and_microbiome_r

ensure_rotl_and_ape

ensure_ancombc
run_downstream
