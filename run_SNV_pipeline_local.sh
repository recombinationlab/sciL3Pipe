#!/bin/bash

# Usage: bash run_pipeline_local.sh

snakemake \
--snakefile Snakefile_SNV \
--use-conda \
--cores 2 \
--configfile configs/config_SNV.yaml
