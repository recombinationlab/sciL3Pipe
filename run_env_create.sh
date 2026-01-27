#!/bin/bash

snakemake \
--configfile configs/config.yaml \
--containerize > Dockerfile

# docker build -t scil3pipe -f Dockerfile .
# apptainer build containers/scil3pipe.sif docker-daemon://scil3pipe:latest

