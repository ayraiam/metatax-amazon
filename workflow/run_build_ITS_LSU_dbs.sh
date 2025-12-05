#!/usr/bin/env bash
# ==========================================================
# Script: workflow/run_build_ITS_LSU_dbs.sh
# Purpose: Build Emu ITS (UNITE) and LSU (SILVA) databases
#          for Ascomycota + Basidiomycota, if not present.
# Notes:
#   - Intended to be called from workflow/runall.sh (via srun)
#   - Uses the same "emu-env" conda environment as run_emu_amplicons.sh
#   - if DBs already exist, it just prints a message and exits.
# ==========================================================
set -euo pipefail

# --------- CONFIG FROM ENV OR DEFAULTS --------------------
# ITS UNITE DB target directory (must contain species_taxid.fasta + taxonomy.tsv)
EMU_DB_ITS_DIR="${EMU_DB_ITS_DIR:-refdb/unite_its_2025/its_unite_asco_basi}"

# UNITE FASTA you downloaded
ITS_FASTA="${ITS_FASTA:-refdb/unite_its_2025/sh_general_release_dynamic_19.02.2025.fasta}"

# SILVA LSU DB target directory
EMU_DB_LSU_DIR="${EMU_DB_LSU_DIR:-refdb/silva_lsu_138/lsu_silva_asco_basi}"

# SILVA LSU taxonomy FASTA (with taxonomy in header)
LSU_FASTA="${LSU_FASTA:-refdb/silva_lsu_138/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz}"

# Conda env name for Emu (same as run_emu_amplicons.sh)
ENV_NAME="${ENV_NAME:-emu-env}"

LOGDIR="${LOGDIR:-logs}"
mkdir -p "${LOGDIR}"
RUN_LOG="${LOGDIR}/run_build_ITS_LSU_dbs.$(date +%Y%m%d_%H%M%S).log"
touch "${RUN_LOG}"

log()  { echo -e ">>> $*" | tee -a "${RUN_LOG}" >&2; }
warn() { echo -e "!!! $*" | tee -a "${RUN_LOG}" >&2; }
die()  { echo -e "xxx $*" | tee -a "${RUN_LOG}" >&2; exit 1; }

ensure_channels() {
  conda config --remove-key channels >/dev/null 2>&1 || true
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda config --add channels defaults
  conda config --set channel_priority strict
  conda clean --index-cache -y >/dev/null 2>&1 || true
}

create_env_emu() {
  if ! command -v conda >/dev/null 2>&1; then
    die "conda not found. Install Miniforge/Conda and re-run."
  fi
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"

  if conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    log "Environment '${ENV_NAME}' exists. Activating..."
    conda activate "${ENV_NAME}"
  else
    log "Environment '${ENV_NAME}' not found. Creating..."
    ensure_channels
    if command -v mamba >/dev/null 2>&1; then
      SOLVER=mamba
    else
      SOLVER=conda
    fi
    set +e
    ${SOLVER} create -n "${ENV_NAME}" -c conda-forge -c bioconda \
      python=3.11 emu minimap2 "seqkit>=2.6" "r-base>=4.3" \
      "r-ggplot2>=3.4" "r-data.table" pandas>=1.5 -y
    status=$?
    set -e
    if [ $status -ne 0 ]; then
      die "Failed to create env '${ENV_NAME}' (Emu)."
    fi
    conda activate "${ENV_NAME}"
  fi

  command -v emu >/dev/null 2>&1 || die "emu not available in env '${ENV_NAME}'"
}

build_its_db() {
  if [ -z "${ITS_FASTA}" ] || [ ! -f "${ITS_FASTA}" ]; then
    log "ITS_FASTA not set or file not found; skipping ITS DB build."
    return 0
  fi

  local db_dir="${EMU_DB_ITS_DIR}"
  local db_name
  db_name="$(basename "${db_dir}")"
  local db_parent
  db_parent="$(dirname "${db_dir}")"

  if [ -s "${db_dir}/species_taxid.fasta" ] && [ -s "${db_dir}/taxonomy.tsv" ]; then
    log "ITS DB already present at ${db_dir}; skipping build."
    return 0
  fi

  log "Building ITS UNITE Emu DB at ${db_dir} from ${ITS_FASTA} ..."
  mkdir -p "${db_parent}"
  local prefix="${db_parent}/${db_name}"

  python workflow/build_emu_its_unite.py \
    --fasta "${ITS_FASTA}" \
    --out-prefix "${prefix}"

  pushd "${db_parent}" >/dev/null
    emu build-database "${db_name}" \
      --sequences "${db_name}.fasta" \
      --seq2tax   "${db_name}.seq2tax.map.tsv" \
      --taxonomy-list "${db_name}.taxonomy.tsv"
  popd >/dev/null

  if [ -s "${db_dir}/species_taxid.fasta" ] && [ -s "${db_dir}/taxonomy.tsv" ]; then
    log "ITS DB built successfully at ${db_dir}"
  else
    warn "ITS DB build finished but expected files are missing in ${db_dir}."
  fi
}

build_lsu_db() {
  if [ -z "${LSU_FASTA}" ] || [ ! -f "${LSU_FASTA}" ]; then
    log "LSU_FASTA not set or file not found; skipping LSU DB build."
    return 0
  fi

  local db_dir="${EMU_DB_LSU_DIR}"
  local db_name
  db_name="$(basename "${db_dir}")"
  local db_parent
  db_parent="$(dirname "${db_dir}")"

  if [ -s "${db_dir}/species_taxid.fasta" ] && [ -s "${db_dir}/taxonomy.tsv" ]; then
    log "LSU DB already present at ${db_dir}; skipping build."
    return 0
  fi

  log "Building LSU SILVA Emu DB at ${db_dir} from ${LSU_FASTA} ..."
  mkdir -p "${db_parent}"
  local prefix="${db_parent}/${db_name}"

  python workflow/build_emu_lsu_silva.py \
    --fasta "${LSU_FASTA}" \
    --out-prefix "${prefix}"

  pushd "${db_parent}" >/dev/null
    emu build-database "${db_name}" \
      --sequences "${db_name}.fasta" \
      --seq2tax   "${db_name}.seq2tax.map.tsv" \
      --taxonomy-list "${db_name}.taxonomy.tsv"
  popd >/dev/null

  if [ -s "${db_dir}/species_taxid.fasta" ] && [ -s "${db_dir}/taxonomy.tsv" ]; then
    log "LSU DB built successfully at ${db_dir}"
  else
    warn "LSU DB build finished but expected files are missing in ${db_dir}."
  fi
}

#Calling functions
log "Starting ITS/LSU DB builder..."
create_env_emu
build_its_db
build_lsu_db
log "ITS/LSU DB builder done."
