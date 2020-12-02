#!/bin/bash

# Usage: bash run_pipeline_local.sh

snakemake \
--snakefile Snakefile \
--use-conda \
--cores 2 \
--configfile configs/config.yaml
