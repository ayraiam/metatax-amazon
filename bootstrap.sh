#!/usr/bin/env bash
# ==========================================================
# Script: bootstrap.sh
# Purpose: Prepare the environment and directory structure
#          so the MetaTax-Amazon pipelines (libsQC + Emu Amplicons)
#          can run immediately after clone.
# Usage:   bash bootstrap.sh
# ==========================================================

set -euo pipefail

echo ">>> Bootstrapping MetaTax-Amazon project..."

# ----------------------------------------------------------
# 1) Create runtime directories (not tracked by git)
# ----------------------------------------------------------
for dir in data results logs refdb; do
  if [ ! -d "$dir" ]; then
    echo ">>> Creating $dir/"
    mkdir -p "$dir"
  else
    echo ">>> Directory $dir/ already exists."
  fi
done

# ----------------------------------------------------------
# 2) Check for conda/mamba installation
# ----------------------------------------------------------
if command -v mamba >/dev/null 2>&1; then
  MAMBA=mamba
elif command -v conda >/dev/null 2>&1; then
  MAMBA=conda
else
  echo "!!! Neither mamba nor conda found. Please install Miniforge or Mambaforge first."
  exit 1
fi

# ----------------------------------------------------------
# 3) Create the 'libsQC' environment if missing
# ----------------------------------------------------------
if ! $MAMBA env list | grep -qE '^libsQC\s'; then
  if [ -f envs/libsQC.yml ]; then
    echo ">>> Creating environment 'libsQC' from envs/libsQC.yml..."
    $MAMBA env create -f envs/libsQC.yml
  else
    echo ">>> No envs/libsQC.yml found. The run_libsQC.sh script will create the env automatically later."
  fi
else
  echo ">>> Environment 'libsQC' already exists."
fi

# # ----------------------------------------------------------
# # 4) (Optional) Create the 'emu-env' environment if missing
# # ----------------------------------------------------------
# # <<< ADDED (supports Emu pipeline right away)
# if ! $MAMBA env list | grep -qE '^emu-env\s'; then
#   if [ -f envs/emu-env.yml ]; then
#     echo ">>> Creating environment 'emu-env' from envs/emu-env.yml..."
#     $MAMBA env create -f envs/emu-env.yml
#   else
#     echo ">>> No envs/emu-env.yml found. The run_emu_amplicons.sh script will create it automatically later."
#   fi
# else
#   echo ">>> Environment 'emu-env' already exists."
# fi

# ----------------------------------------------------------
# 5) Done
# ----------------------------------------------------------
echo ">>> Bootstrap complete!"
echo ">>> You can now run either:"
echo "      bash workflow/runall.sh           # libsQC (Slurm submission)"
echo "   or"
echo "      bash workflow/run_emu_amplicons.sh # Emu classification stage"
