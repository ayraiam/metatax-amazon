#!/usr/bin/env bash
# ==========================================================
# Script: workflow/run_downstream_analysis.sh
# Purpose:
#   - Create/activate a clean env for downstream analysis
#   - Ensure ANCOMBC (for ANCOM-BC2) is available (GitHub preferred)
#   - Run downstream_analysis.R
#
# Key fixes for your cases:
#   (2) Prevent solver downgrading R after creation (hard re-pin r-base=4.5.*)
#   (3) Detect/stop if PATH is not using env's Rscript (prove which R is used)
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
# Robust R>=4.5 assertion + PATH correctness verification
# ==========================================================
assert_r_version_and_path() {
  log "Checking R version + making sure we're using the ENV's Rscript..."

  local expected_rscript
  expected_rscript="$(conda info --base)/envs/${ENV_NAME}/bin/Rscript"
  if [[ -n "${CONDA_PREFIX:-}" ]]; then
    expected_rscript="${CONDA_PREFIX}/bin/Rscript"
  fi
  log "Expected env Rscript: ${expected_rscript}"

  local actual_rscript
  actual_rscript="$(command -v Rscript || true)"
  log "Actual   Rscript on PATH: ${actual_rscript}"

  if [[ -z "${actual_rscript}" ]]; then
    die "Rscript not found on PATH after activation."
  fi

  if [[ "${actual_rscript}" != "${expected_rscript}" ]]; then
    log "DEBUG: CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
    log "DEBUG: PATH=${PATH}"
    die "You are NOT using the env's Rscript (PATH/module issue). Fix your module environment or ensure conda activate works inside the batch job."
  fi

  local rv
  rv="$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || echo "0.0.0")"
  log "R version in env: ${rv}"

  if ! Rscript -e 'quit(status=ifelse(getRversion() >= "4.5.0", 0, 1))'; then
    die "R < 4.5.0 detected in env. Something downgraded R or activation is not real."
  fi

  log "Conda r-base entry:"
  conda list -n "${ENV_NAME}" | awk 'BEGIN{found=0} $1=="r-base"{print; found=1} END{if(!found) print "r-base not found in conda list"}' >&2

  log "R version OK (>= 4.5.0) and PATH is correct."
}

# ==========================================================
# Hard pin R 4.5.* AFTER ANY install (prevents silent downgrade)
# ==========================================================
repin_r_base() {
  log "Re-pinning r-base to 4.5.* to prevent solver downgrades..."

  if command -v mamba >/dev/null 2>&1; then
    # mamba may NOT support --update-specs. Use plain install.
    log "mamba install -n ${ENV_NAME} -c conda-forge -c bioconda 'r-base=4.5.*'"
    mamba install -n "${ENV_NAME}" -c conda-forge -c bioconda "r-base=4.5.*" -y
  else
    # conda DOES support --update-specs; keep it here
    log "conda install -n ${ENV_NAME} -c conda-forge -c bioconda 'r-base=4.5.*' --update-specs"
    conda install -n "${ENV_NAME}" -c conda-forge -c bioconda "r-base=4.5.*" --update-specs -y
  fi

  assert_r_version_and_path
}

# ==========================================================
# install build toolchain + system libs needed to compile Bioconductor from source
# ==========================================================
install_build_deps_conda() {
  log "Ensuring build/system dependencies exist in the env (for compiling Bioconductor packages from source)..."

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
    llvm-openmp
    zlib
    bzip2
    xz
    libcurl
    openssl
    libxml2
    pcre2
    icu
    readline
  )

  if command -v mamba >/dev/null 2>&1; then
    log "mamba install -n ${ENV_NAME} -c conda-forge ${pkgs[*]}"
    mamba install -n "${ENV_NAME}" -c conda-forge "${pkgs[@]}" -y
  else
    log "conda install -n ${ENV_NAME} -c conda-forge ${pkgs[*]}"
    conda install -n "${ENV_NAME}" -c conda-forge "${pkgs[@]}" -y
  fi

  log "Build/system deps installed."
  repin_r_base
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
  repin_r_base
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
  repin_r_base
  install_build_deps_conda
}

# ----------------------------------------------------------
# Install Bioconductor deps via BiocManager (inside R)
# ----------------------------------------------------------
install_bioc_deps_r() {
  log "Installing required Bioconductor dependencies via BiocManager (R-side, matches your R=4.5.*)..."
  assert_r_version_and_path

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
  assert_r_version_and_path
}

# ----------------------------------------------------------
# ensure ANCOMBC is installed (GitHub preferred)
# ----------------------------------------------------------
ensure_ancombc() {
  log "Checking for R package 'ANCOMBC' (ANCOM-BC2)..."
  assert_r_version_and_path

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

  assert_r_version_and_path

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

  assert_r_version_and_path

  Rscript workflow/downstream_analysis.R \
    "$INFILE" \
    "$OUTDIR" \
    "$BASENAME"

  log "Plots + tables written to: $OUTDIR"
}

create_env
ensure_ancombc
run_downstream
