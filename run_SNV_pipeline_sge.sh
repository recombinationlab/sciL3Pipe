#!/bin/bash

# Usage: qsub run_pipeline_sge.sh

snakemake \
--snakefile Snakefile_SNV \
--use-conda \
--conda-frontend mamba \
--configfile configs/config_SNV.yaml \
-j 100 \
--latency-wait 300 \
--restart-times 1 \
--cluster-config cluster_SNV.yaml \
--cluster "qsub -V -cwd \
-l highp,h_data={cluster.mem},h_rt={cluster.time},nodes={cluster.nodes} \
-pe shared {cluster.cpus} \
-o {cluster.output} \
-e {cluster.error}"

