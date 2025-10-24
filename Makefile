SHELL := /usr/bin/env bash

.PHONY: setup run run-local clean

setup:
	./bootstrap.sh

# Run via Slurm wrapper (your current path)
run:
	bash workflow/runall.sh

# Run locally (no Slurm) for a tiny smoke test
run-local:
	bash workflow/run_libsQC.sh

clean:
	rm -rf logs/* results/*
