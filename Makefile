SHELL := /usr/bin/env bash

.PHONY: setup run run-local run-emu clean  # <<< ADDED run-emu target

# ------------------------------------------------------------
# 1) Bootstrap project (creates dirs, checks envs)
# ------------------------------------------------------------
setup:
	./bootstrap.sh

# ------------------------------------------------------------
# 2) Run libsQC pipeline via Slurm (recommended for HPC)
# ------------------------------------------------------------
run:
	bash workflow/runall.sh

# ------------------------------------------------------------
# 3) Run libsQC locally (no Slurm) for quick smoke test
# ------------------------------------------------------------
run-local:
	bash workflow/run_libsQC.sh

# ------------------------------------------------------------
# 4) Run Emu Amplicons stage (classification)
# ------------------------------------------------------------
run-emu:                                      # <<< ADDED target for Emu pipeline
	bash workflow/run_emu_amplicons.sh        # <<< ADDED command

# ------------------------------------------------------------
# 5) Clean results and logs (does NOT touch refdb/ or config/)
# ------------------------------------------------------------
clean:
	rm -rf logs/* results/*
