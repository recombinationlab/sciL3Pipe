FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="aeb3c9e2c062cd88c3e85be7ce24ddad266575104f18ad478881b4f46c3930a4"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/r_env.yaml
#   prefix: /conda-envs/e2666ede98d6b4f761e29c414afc7b8c
#   name: r_env
#   channels:
#     - bioconda
#     - anaconda
#     - conda-forge
#   dependencies:
#       - r-base=4.3.1
#       - r-argparse=2.2.2
#       - r-ggbeeswarm=0.7.2
#       - r-plotly=4.11.0
#       - r-htmlwidgets=1.6.4
#       - bioconductor-genomicranges=1.52.0
#       - bioconductor-rtracklayer=1.60.0
#       - bioconductor-biocgenerics=0.46.0
#       - bioconductor-iranges=2.34.1
#       - r-dosnow=1.0.20
#       - r-foreach=1.5.2
#       - r-doparallel=1.0.17
#       - r-tidyverse=2.0.0
#       - jq=1.7
RUN mkdir -p /conda-envs/e2666ede98d6b4f761e29c414afc7b8c
COPY workflow/envs/r_env.yaml /conda-envs/e2666ede98d6b4f761e29c414afc7b8c/environment.yaml

# Conda environment:
#   source: workflow/envs/sciL3_env.yaml
#   prefix: /conda-envs/da9669680debec9c1e74030c9717bac4
#   name: sciL3_env
#   channels:
#     - bioconda
#     - anaconda
#     - conda-forge
#   dependencies:
#     - coreutils=9.5
#     - bwa=0.7.18
#     - bowtie2=2.5.4
#     - samtools=1.20
#     - pysam=0.22.1
#     - matplotlib=3.9.2
#     - bedtools=2.31.1
#     - pyfastx=2.1.0
#     - python-levenshtein=0.12.2
#     - natsort=7.1.1
#     - pigz=2.8
#     - cutadapt=4.9
#     - fastqc=0.12.1
#     - seqkit=2.8.2
#     - pandas=2.3.2
RUN mkdir -p /conda-envs/da9669680debec9c1e74030c9717bac4
COPY workflow/envs/sciL3_env.yaml /conda-envs/da9669680debec9c1e74030c9717bac4/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/e2666ede98d6b4f761e29c414afc7b8c --file /conda-envs/e2666ede98d6b4f761e29c414afc7b8c/environment.yaml && \
    conda env create --prefix /conda-envs/da9669680debec9c1e74030c9717bac4 --file /conda-envs/da9669680debec9c1e74030c9717bac4/environment.yaml && \
    conda clean --all -y
