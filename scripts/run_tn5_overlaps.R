#!/usr/bin/env Rscript

if(suppressMessages(!require(argparse, quietly = TRUE))){
    install.packages("argparse", repos = "http://cran.us.r-project.org")
    suppressMessages(library(argparse, quietly = TRUE))
}

if(suppressMessages(!require(devtools, quietly = TRUE))){
    install.packages("devtools", repos = "http://cran.us.r-project.org")
    suppressMessages(library(remotes, quietly = TRUE))
}

if(suppressMessages(!require(BiocManager, quietly = TRUE))){
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
    suppressMessages(library(BiocManager, quietly = TRUE))
}

if(suppressMessages(!require(rtracklayer, quietly = TRUE))){
    BiocManager::install("rtracklayer")
    suppressMessages(library(rtracklayer, quietly = TRUE))
}

parse_arguments <- function(){
    parser <- ArgumentParser(description="Run breakpointR")
    parser$add_argument("-i","--input_bam", metavar="BAM1,BAM2,...", type="character", dest="input_bam",
                        required=TRUE, nargs="+",
                        help="Strand-seq BAM files, one for each single cell")
    parser$add_argument("-o", "--output_folder", dest="output_folder", type="character", required=TRUE,
                        help="Output folder location")
    parser$add_argument("-t", "--threads", dest="threads", type="integer", default = 1,
                        help="Number of threads to use")
    parser$add_argument("-p","--package", dest="package", type="character",
                        help="Path to source of sciStrandR to install if not installed already")

  args <- parser$parse_args()

  return(args)
}



main <- function(){

    args <- parse_arguments()

    require(sciStrandR, quietly = TRUE)
 
    cat("=========================== Running Tn5 overlap finder ===================\n",
        "Number of input BAM(s): ", length(args$input_bam), "\n",
        "Output folder: ", args$output_folder, "\n",
        "Number of threads: ", args$threads, "\n",
        "==========================================================================\n")

    dir.create(file.path(args$output_folder), showWarnings = FALSE)

    bam_counts(args$input_bam, args$output_folder, num_threads = args$threads,
                tn5_overlap=12, ignore_soft_clipping=TRUE,
                mapq=20, offset=3)

    cat("======= Finished =======\n")
  
}

main()
