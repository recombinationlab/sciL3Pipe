
# sciStrandPipe

A Snakemake pipeline for processing and aligning sci-L3-seq data.


## Running pipeline from raw fastq files:

1. (Optional) Create a folder with symbolic links to fastq files
    ```mkdir fastqs```
    ```ln -s <path to raw fastq> <path to fastqs folder>```

2. Create sample.json file that informs the pipeline which fastq files correspond to which sample:
    ```python /u/project/yeastyin/chovanec/pipelines/sciStrandPipe/fastq2json.py --fastq_dir /u/project/yeastyin/chovanec/lib293/fastqs/```

    Will create a file with the following structure:
    ```
    {
    "yi293": {
        "R1": [
            "/u/project/yeastyin/chovanec/lib293/fastqs/yi293_S7_R1_001.fastq.gz"
        ],
        "R2": [
            "/u/project/yeastyin/chovanec/lib293/fastqs/yi293_S7_R2_001.fastq.gz"
        ]
        }
    }
    ```

3. Create a config.yaml file specific to your sample. The config file contains pipeline input parameters. Example config file for sample yi293 is shown below, with parameters likely to change between runs highlighted with `<---`:

```
# BWA or bowtie2
aligner: "bwa"                              <------- Choose which aligner to use
# output directory (include trailing /)
output_dir: "/u/project/yeastyin/chovanec/lib293/" <------- Set output directory
assembly: "hg38" <--- Assembly to align data to (corresponds to fields in indexes)
################################################################################
# SSS extract and split
################################################################################
# fastq file location produced by fastq2json.py
fastqs: "/u/project/yeastyin/chovanec/lib293/samples.json"  <--- Path to json
# samples with corresponding SSS barcode file, must match fastqs
samples: 
    yi293: "barcodes/barcode_sss_lib293.txt" <------ Path to sample SSS barcodes
    <----------------------------------------- Can specify multiple samples here
sss_mismatches: 0
rt_mismatches: 3
rt_primer: GGGATGCAGCTCGCTCCTG
################################################################################
# Tn5 and ligation barcode
################################################################################
# the round1 barcode list (on tn5) for splitting single cells
tn5_barcodes: "barcodes/barcode_tn5.txt"
# the round2 barcode list (by ligation) for splitting single cells
ligation_barcodes: "barcodes/barcode_ligation_104.txt"
tn5_mismatches: 1
ligation_mismatches: 1
# Reads per single cell to keep
cutoff: 10000
################################################################################   
# Indexes
################################################################################
hybrid_reference: True
bwa_index:
    hg38: "/u/project/yeastyin/chovanec/genomes/hybrid/GRCh38.p13_GRCm38.p6.chr"
bowtie2_index:
    hg38: "/u/project/yeastyin/chovanec/genomes/hybrid/GRCh38.p13_GRCm38.p6.chr"
```

4. To run the pipeline on an sge cluster, use the `run_pipeline_sge.sh` script:

```
#!/bin/bash

snakemake \
--snakefile Snakefile \
--use-conda \
--configfile configs/config_lib293.yaml \ <------ Will need to edit path to your config file
-j 100 \
--cluster-config cluster.yaml \
--cluster "qsub -V -cwd \
-l highp,h_data={cluster.mem},h_rt={cluster.time},nodes={cluster.nodes} \
-pe shared {cluster.cpus} \
-o {cluster.output} \
-e {cluster.error}"
```
To run:
```
bash run_pipeline_sge.sh
```

## Running pipeline from trimmed and attached fastq files:

1. Create a directory structure `trimmed/<sample>` what will contain symbolic links to trimmed and attached fastq files. 
    ```mkdir -p fastqs/yi293```
    ```ln -s <path to trimmed attached fastq> <path to fastqs folder>```

    The trimmed files are expected to have a specific formated file name and be gzip'ed (the pipeline will not recognize the files if named differently) e.g.:
    ```
    yi293_GAACCG.R1.trimmed.attached.fastq.gz
    yi293_GAACCG.R2.trimmed.attached.fastq.gz
    ```
    ```
    <sample>_<SSS barcode>.<R1 or R2>.trimmed.attached.fastq.gz
    ```

2. Create a config.yaml file specific to your sample. The config file is the same as above with the exception of omitting the fastqs json parameter (also you will notice the json generation step has not been run):

```
# BWA or bowtie2
aligner: "bwa"                              
# output directory (include trailing /)
output_dir: "/u/project/yeastyin/chovanec/lib293/" 
assembly: "hg38"
################################################################################
# SSS extract and split
################################################################################
# fastq file location produced by fastq2json.py
# fastqs: "/u/project/yeastyin/chovanec/lib293/samples.json"  <--- Comment out 
# samples with corresponding SSS barcode file, must match fastqs
samples: 
    yi293: "barcodes/barcode_sss_lib293.txt"
sss_mismatches: 0
rt_mismatches: 3
rt_primer: GGGATGCAGCTCGCTCCTG
################################################################################
# Tn5 and ligation barcode
################################################################################
# the round1 barcode list (on tn5) for splitting single cells
tn5_barcodes: "barcodes/barcode_tn5.txt"
# the round2 barcode list (by ligation) for splitting single cells
ligation_barcodes: "barcodes/barcode_ligation_104.txt"
tn5_mismatches: 1
ligation_mismatches: 1
# Reads per single cell to keep
cutoff: 10000
################################################################################   
# Indexes
################################################################################
hybrid_reference: True
bwa_index:
    hg38: "/u/project/yeastyin/chovanec/genomes/hybrid/GRCh38.p13_GRCm38.p6.chr"
bowtie2_index:
    hg38: "/u/project/yeastyin/chovanec/genomes/hybrid/GRCh38.p13_GRCm38.p6.chr"
```

3. To run the pipeline on an sge cluster, use the `run_pipeline_sge.sh` script:

```
#!/bin/bash

snakemake \
--snakefile Snakefile \
--use-conda \
--configfile configs/config_lib293.yaml \ <------ Will need to edit path to your config file
-j 100 \
--cluster-config cluster.yaml \
--cluster "qsub -V -cwd \
-l highp,h_data={cluster.mem},h_rt={cluster.time},nodes={cluster.nodes} \
-pe shared {cluster.cpus} \
-o {cluster.output} \
-e {cluster.error}"
```
To run:
```
bash run_pipeline_sge.sh
```

## General troubleshooting

Sometimes snakemake will not recognize output files and produce an error:
```
This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait. completed successfully, but some output files are missing.
```
It is safe to just restart the pipeline with the run script, additionally adding `--latency-wait <seconds>` to the run script may resolve future errors.