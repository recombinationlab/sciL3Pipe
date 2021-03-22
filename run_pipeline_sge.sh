#!/bin/bash

# Usage: qsub run_pipeline_sge.sh

snakemake \
--snakefile Snakefile \
--use-conda \
--conda-frontend mamba \
--configfile configs/config.yaml \
-j 100 \
--latency-wait 300 \
--cluster-config cluster.yaml \
--cluster "qsub -V -cwd \
-l highp,h_data={cluster.mem},h_rt={cluster.time},nodes={cluster.nodes} \
-pe shared {cluster.cpus} \
-o {cluster.output} \
-e {cluster.error}"

