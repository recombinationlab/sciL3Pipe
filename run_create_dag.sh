#!/bin/bash
# requires graphviz to be installed

snakemake \
--profile profiles/local_profile/ \
--configfile configs/config.yaml \
--dag | dot -Tsvg > dag/dag.svg


