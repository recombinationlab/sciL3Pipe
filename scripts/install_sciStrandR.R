#!/usr/bin/env Rscript

if(suppressMessages(!require(argparse, quietly = TRUE))){
  install.packages("argparse", repos = "http://cran.us.r-project.org")
  suppressMessages(library(argparse, quietly = TRUE))
}

if(suppressMessages(!require(remotes, quietly = TRUE))){
    install.packages("remotes", repos = "http://cran.us.r-project.org")
    suppressMessages(library(remotes, quietly = TRUE))
}

if(suppressMessages(!require(BiocManager, quietly = TRUE))){
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
  suppressMessages(library(BiocManager, quietly = TRUE))
}


# TODO: work around to DO.db getting installed under R version 3 for some reason
package_list <- data.frame(installed.packages())
if(!(package_list$Built[grepl("DO.db", package_list$Package)] > 4)){
  BiocManager::install("DO.db", force = TRUE)
}




parse_arguments <- function(){
  parser <- ArgumentParser(description="Install sciStrandR")

  parser$add_argument("-p","--package", dest="package", type="character",
                      help="Path to source of sciStrandR to install if not installed already")
  parser$add_argument("-o","--output", dest="output", type="character",
                      help="If successful, will write out a file.")

  args <- parser$parse_args()

  return(args)
}

main <- function(){

    args <- parse_arguments()


    if(suppressMessages(!require(sciStrandR, quietly = TRUE))){
        remotes::install_local(args$package, dependencies=TRUE, upgrade=FALSE, repos=BiocManager::repositories())
        suppressMessages(library(sciStrandR, quietly = TRUE))
    }

    cat("Install finished", file=args$output, sep="\n")

}

main()