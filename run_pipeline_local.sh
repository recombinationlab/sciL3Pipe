#!/bin/bash

# Usage: bash run_pipeline_local.sh

snakemake \
--snakefile Snakefile \
--use-conda \
--conda-frontend mamba \
--cores 2 \
--configfile configs/config_yi331.yaml

# TODO: --conda-prefix /path/to/shared/environment
# TODO: working dir