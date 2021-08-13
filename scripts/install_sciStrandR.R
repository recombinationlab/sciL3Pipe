#!/usr/bin/env Rscript

if(suppressMessages(!require(argparse, quietly = TRUE))){
  install.packages("argparse")
  suppressMessages(library(argparse, quietly = TRUE))
}

if(suppressMessages(!require(devtools, quietly = TRUE))){
    install.packages("devtools", repos = "http://cran.us.r-project.org")
    suppressMessages(library(remotes, quietly = TRUE))
}


parse_arguments <- function(){
  parser <- ArgumentParser(description="Run breakpointR")

  parser$add_argument("-p","--package", dest="package", type="character",
                      help="Path to source of breakpointRAddon to install if not installed already")

  args <- parser$parse_args()

  return(args)
}

main <- function(){

    args <- parse_arguments()


    if(suppressMessages(!require(breakpointRAddon, quietly = TRUE))){
        devtools::install(args$package, dependencies = TRUE)
        suppressMessages(library(breakpointRAddon, quietly = TRUE))
    }

    cat("Install finished")

}

main()