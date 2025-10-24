#!/usr/bin/env bash
# ==========================================================
# Script: bootstrap.sh
# Purpose: Prepare the environment and directory structure
#          so the libsQC pipeline can run immediately after clone.
# Usage:   bash bootstrap.sh
# ==========================================================

set -euo pipefail

echo ">>> Bootstrapping Meta-Amazon project..."

# ----------------------------------------------------------
# 1) Create runtime directories (not tracked by git)
# ----------------------------------------------------------
for dir in data results logs; do
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

# ----------------------------------------------------------
# 4) Done
# ----------------------------------------------------------
echo ">>> Bootstrap complete!"
echo ">>> You can now run:"
echo "      bash workflow/run_libsQC.sh"
echo "   or"
echo "      bash workflow/runall.sh  (for Slurm)"
