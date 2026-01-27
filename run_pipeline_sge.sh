#!/bin/bash

unset LD_RUN_PATH_modshare
unset LD_RUN_PATH
unset LD_LIBRARY_PATH_modshare

snakemake \
--snakefile Snakefile \
--profile profiles/hoffman2_profile/ \
--configfile configs/config.yaml \
--apptainer-args "-B /u/project/yeastyin/chovanec/:/u/project/yeastyin/chovanec/"
