#!/bin/bash

# Usage: bash run_pipeline_local.sh

snakemake \
--snakefile Snakefile_strand \
--use-conda \
--conda-frontend mamba \
--cores 2 \
--configfile configs/config_strand.yaml
