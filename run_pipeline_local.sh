#!/bin/bash

snakemake \
--configfile configs/config.yaml \
--apptainer-args "-B /mnt/data/Projects_sync/sciL3Pipe:/mnt/data/Projects_sync/sciL3Pipe,/media/chovanec/NVMe/genomes/GRCh38:/media/chovanec/NVMe/genomes/GRCh38" \
--dry-run --debug-dag --notemp 
#--apptainer-args "-B /mnt/d/SynologyDriveCloud/Projects/sciL3Pipe/tests/:/mnt/d/SynologyDriveCloud/Projects/sciL3Pipe/tests/,/mnt/e/genomes/GRCh38/:/genomes/GRCh38/" \


# --apptainer-args
# --singularity-args "-B /mnt/data/Projects_sync/sciL3Pipe:/mnt/data/Projects_sync/sciL3Pipe,/media/chovanec/NVMe/genomes/GRCh38:/media/chovanec/NVMe/genomes/GRCh38" \



