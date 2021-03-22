#!/usr/bin/env Rscript

if(suppressMessages(!require(argparse, quietly = TRUE))){
  install.packages("argparse")
  suppressMessages(library(argparse, quietly = TRUE))
}

if(suppressMessages(!require(breakpointR, quietly = TRUE))){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("breakpointR")
  suppressMessages(library(breakpointR, quietly = TRUE))
}

parse_arguments <- function(){
  parser <- ArgumentParser(description="Run breakpointR")
  parser$add_argument("-i","--input_folder", metavar="PATH", type="character", dest="input_folder",
                      required=TRUE,
                      help="Strand-seq folder containing BAM files, one for each single cell")
  parser$add_argument("-o", "--output_folder", dest="output_folder", type="character", required=TRUE,
                      help="Output folder location")
  parser$add_argument("-t", "--threads", dest="threads", type="integer",
                      help="Number of threads to use")
  parser$add_argument("--type", dest="type", type="character", required=TRUE,
                      choices=c("haploid", "diploid"),
                      help="Output folder location")
  parser$add_argument("-b","--blacklist", dest="blacklist", type="character",
                      help="Path to blacklist BED file")
  
  args <- parser$parse_args()
  
  return(args)
}


main <- function(){
  
  args <- parse_arguments()
  
  # args <- parse_arguments(c("-i", "/mnt/data/nextseq190419/yi293_GAACCG_1min_UV_with_USER_HAP1/split/subset/filtered_0.5/",
  #                      "-o", "/mnt/data/sci-l3/test/breakpoint_hg19_bwa_b0.05_p2fF_mr10_pt0.33_t10_1Mb/",
  #                      "-t", "4",
  #                      "--type","haploid",
  #                      "-b", "/mnt/data/genomes/hg19-blacklist.v2.ensembl.bed"))
  tryCatch({
  if(args$type == "haploid"){
    
    cat("======= Running BreakpointR with setting optimised for haploid cells =======\n", 
        "Input folder: ", args$input_folder, "\n",
        "Output folder: ", args$output_folder, "\n",
        "windowsize = 1000000\n", 
        "binMethod = size\n",
        "pairedEndReads = TRUE\n", 
        "pair2frgm = FALSE\n", 
        "min.mapq = 20\n", 
        "filtAlt = TRUE\n", 
        "peakTh = 0.33\n", 
        "trim = 10\n",
        "background = 0.05\n", 
        "minReads = 10\n",
        "===========================================================================\n")
    
    breakpointr(inputfolder = args$input_folder,
                outputfolder = args$output_folder,
                windowsize = 1000000, binMethod = "size",
                pairedEndReads = TRUE, pair2frgm = FALSE, min.mapq = 20, filtAlt = TRUE, peakTh = 0.33, trim = 10,
                background = 0.05, minReads = 10, numCPU = args$threads, maskRegions=args$blacklist)
    
  }else if(args$type == "diploid"){
    
    cat("======= Running BreakpointR with setting optimised for diploid cells =======\n", 
        "Input folder: ", args$input_folder, "\n",
        "Output folder: ", args$output_folder, "\n",
        "windowsize = 5000000\n", 
        "binMethod = size\n",
        "pairedEndReads = TRUE\n", 
        "pair2frgm = FALSE\n", 
        "min.mapq = 20\n", 
        "filtAlt = TRUE\n", 
        "peakTh = 0.33\n", 
        "trim = 10\n",
        "background = 0.05\n", 
        "minReads = 50\n",
        "===========================================================================\n")
    
    breakpointr(inputfolder = args$input_folder,
                outputfolder = args$output_folder,
                windowsize = 5000000, binMethod = "size",
                pairedEndReads = TRUE, pair2frgm = FALSE, min.mapq = 20, filtAlt = TRUE, peakTh = 0.33, trim = 10,
                background = 0.05, minReads = 50, numCPU = args$threads, maskRegions=args$blacklist)
  }},
  error=function(cond){
    message("Error during breakpointR execution, some output files may still have been produced.")
  })
}


main()