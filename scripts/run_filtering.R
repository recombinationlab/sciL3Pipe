#!/usr/bin/env Rscript

if(suppressMessages(!require(argparse, quietly = TRUE))){
  install.packages("argparse")
  suppressMessages(library(argparse, quietly = TRUE))
}

# if(suppressMessages(!require(remotes, quietly = TRUE))){
#   install.packages("remotes", repos = "http://cran.us.r-project.org")
#   suppressMessages(library(remotes, quietly = TRUE))
# }

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
  parser$add_argument("-i","--input_folder", metavar="PATH", type="character", dest="input_folder",
                      required=TRUE,
                      help="Strand-seq folder containing BAM files, one for each single cell")
  parser$add_argument("-o", "--output_folder", dest="output_folder", type="character", required=TRUE,
                      help="Output folder location")
  parser$add_argument("-s", "--summary_output", dest="summary_output", type="character", required=TRUE,
                      help="Summary file output path")
  parser$add_argument("-t", "--threads", dest="threads", type="integer", default = 1,
                      help="Number of threads to use")
  parser$add_argument("--type", dest="type", type="character", required=TRUE,
                      choices=c("haploid", "diploid"),
                      help="Output folder location")
  parser$add_argument("-a", "--assembly", dest="assembly", type="character", required=TRUE,
                      choices=c("hg19", "mm10", "hg38"),
                      help="Genome assembly, used to fetch centromere coordinates")
  parser$add_argument("-d", "--distance_cutoff", dest="distance_cutoff", type="integer",
                      default=2000000,
                      help="Double event distance cutoff")
  parser$add_argument("-c", "--centromere_distance", dest="centromere_distance", type="integer",
                      default=3000000,
                      help="Centromere proximal double event distance cutoff")
  parser$add_argument("-p","--package", dest="package", type="character",
                      help="Path to source of breakpointRAddon to install if not installed already")
  parser$add_argument("-f", "--filter", dest="filter", type="character", default=NULL,
                      help="Filter whole chromosome or a single location")
  parser$add_argument("-b", "--bed", dest="bed", type="character",
                      help="BED file with centromere coordinates")

  args <- parser$parse_args()

  return(args)
}



main <- function(){

  args <- parse_arguments()

  # if(suppressMessages(!require(breakpointRAddon, quietly = TRUE))){
  #   remotes::install_local(args$package, dependencies=TRUE, upgrade=FALSE, repos=BiocManager::repositories())
  #   suppressMessages(library(breakpointRAddon, quietly = TRUE))
  # }


# args <- parse_arguments(c("-i", "/mnt/data/tmp/data/",
#                      "-o", "/mnt/data/tmp/data_filtered/",
#                      "-s", "/mnt/data/tmp/summary.csv",
#                      "-t", "1",
#                      "--type","diploid",
#                      "-a", "mm10",
#                      "-d", "5000000",
#                      "-c", "5000000"))

  cat("=========================== Running filtering  ===========================\n",
      "Input folder: ", args$input_folder, "\n",
      "Output folder: ", args$output_folder, "\n",
      "Distance cutoff: ", args$distance_cutoff, "\n",
      "Centromere distance cutoff: ", args$centromere_distance, "\n",
      "Type: ", args$type, "\n",
      "Filtering: ", args$filter, "\n",
      "==========================================================================\n")

  if(args$type == "haploid"){
    gtype <- TRUE
  }else{
    gtype <- FALSE
  }

  if(!is.null(args$bed)){
    cent_coord <- rtracklayer::import.bed(args$bed)
  }else{
    if(args$assembly == "hg19"){
      cent_coord <- get_centromere_coordinates("hg19", ensembl = TRUE)
    }else if(args$assembly == "hg38"){
      cent_coord <- get_centromere_coordinates("hg38", reduce_single = TRUE, ensembl = TRUE)
    }else if(args$assembly == "mm10"){
      cent_coord <- get_centromere_coordinates("mm10", ensembl = FALSE)
    }
  }

  f_summary <- breakpoint_filter(args$input_folder, args$output_folder,
                                distance_cutoff=args$distance_cutoff,
                                cent=cent_coord,
                                filt=args$filt, assembly=args$assembly,
                                centromere_distance_cutoff=args$centromere_distance,
                                num_threads=args$threads,
                                haploid=gtype)

  cat("======= Creating summary =======\n",
      "Input folder: ", args$input_folder, "\n",
      "===========================================================================\n")

  b_summary <- breakpoint_summary(args$output_folder, num_threads=args$threads)
  write.csv(b_summary,
            args$summary_output,
            quote = FALSE, row.names = FALSE)

  cat("======= Finished =======\n")

}

main()
