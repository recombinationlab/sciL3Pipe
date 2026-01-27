#!/usr/bin/env Rscript

# setwd("D://SynologyDriveCloud/Projects/sciL3Pipe/workflow/scripts/")

# if(suppressMessages(!require(argparse, quietly = TRUE))){
#   install.packages("argparse", repos = "http://cran.us.r-project.org")
#   suppressMessages(library(argparse, quietly = TRUE))
# }
# 
# if(suppressMessages(!require(parallel, quietly = TRUE))){
#   install.packages("parallel", repos = "http://cran.us.r-project.org")
#   suppressMessages(library(parallel, quietly = TRUE))
# }
# 
# if(suppressMessages(!require(doSNOW, quietly = TRUE))){
#   install.packages("doSNOW", repos = "http://cran.us.r-project.org")
#   suppressMessages(library(doSNOW, quietly = TRUE))
# }

library(argparse)
library(doSNOW)
library(parallel)
library(ggbeeswarm)
library(dplyr)
library(ggplot2)
# library(plotly)
# library(htmlwidgets)
library(GenomicRanges)
library(base64enc)


parse_arguments <- function(){
  parser <- ArgumentParser(description="Run background estimate")
  parser$add_argument("-i","--input_dir", type="character", dest="input_dir",
                      required=TRUE, nargs="+",
                      help="Strand-seq BAM files, one for each single cell")
  parser$add_argument("-o", "--output_csv", dest="output_file", type="character", 
                      required=TRUE,
                      help="Output csv location")
  parser$add_argument("-o2", "--output_html", dest="output_html", type="character",
                      required=TRUE,
                      help="Output html location")
  parser$add_argument("-t", "--threads", dest="threads", type="integer", default = 1,
                      help="Number of threads to use")
  parser$add_argument("-w", "--bin_width", dest="bin_width", type="integer", default = 10000000,
                      help="Bin width used to estimate background")

  args <- parser$parse_args()
  
  return(args)
}



main <- function(){
  
  args <- parse_arguments()
  
  bkg_df <- run_background_estimation(args$input_dir, num_threads = args$threads, bin_width = args$bin_width)
  
  write.csv(bkg_df, file = args$output_file, row.names = FALSE)
  
  bkg_plot(bkg_df, out_html = args$output_html)
  
}


# bam_path <- "D://SynologyDriveCloud/Projects/sciStrandR/inst/extdata/annotation/yi305_GTCTAT.CCGCGCTTTGATGCA.human.filt.bam"
# 
# data_folder <- "D://SynologyDriveCloud/Projects/sciStrandR/inst/extdata/annotation/"
# 
# bkg_df <- run_background_estimation(data_folder, num_threads = 2)
# bkg_plot(bkg_df, out_html = "test.html")


run_background_estimation <- function(data_folder, num_threads=1,
                                      bin_width=10e+06, limit=NULL,
                                      basename_pattern="^([^\\.]*\\.[^\\.]*).*"){


  # get bam files
  bam_files <- list.files(data_folder, pattern = "*.bam$", full.names = TRUE, recursive = TRUE)
  
  if(!is.null(limit)){
    bam_files <- bam_files[1:limit]
  }

  if(num_threads > 1){

    # progress bar
    pb <- txtProgressBar(max = length(bam_files), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    message("\nUsing ",num_threads," threads")
    cl <- parallel::makeCluster(num_threads)
    doSNOW::registerDoSNOW(cl)

    out_df <- foreach(i=1:length(bam_files), .packages = ("GenomicRanges"),
                      .export = c("readBamFileAsGRanges", "background_estimation_on_bins", "state_assignment_bins"),
                      .options.snow = opts, .combine=rbind) %dopar% {
                        background_estimation_on_bins(bam_files[i], bin_width=bin_width, basename_pattern=basename_pattern)
                      }

    parallel::stopCluster(cl)

  }else{
    out_lst <- list()
    for(data_path in bam_files){

      out_lst[[data_path]] <- background_estimation_on_bins(data_path, bin_width=bin_width, basename_pattern=basename_pattern)

    }
    out_df <- dplyr::bind_rows(out_lst)

  }

  rownames(out_df) <- NULL

  return(out_df)

}


state_assignment_bins <- function(bins, frags){
  
  # get counts for each segment
  ranges_tiles <- bins
  strand(ranges_tiles) <- '+'
  ranges_tiles$Ws <- GenomicRanges::countOverlaps(ranges_tiles, frags)
  strand(ranges_tiles) <- '-'
  ranges_tiles$Cs <- GenomicRanges::countOverlaps(ranges_tiles, frags)
  strand(ranges_tiles) <- '*'
  
  # genotype each segment
  tile_wc_ratio <- (ranges_tiles$Ws-ranges_tiles$Cs)/(ranges_tiles$Ws+ranges_tiles$Cs)
  ranges_tiles$tile_wc_ratio <- tile_wc_ratio
  
  ranges_tiles$states <- "?"
  ranges_tiles$states[tile_wc_ratio > 0.8] <- "ww"
  ranges_tiles$states[tile_wc_ratio < -0.8] <- "cc"
  ranges_tiles$states[tile_wc_ratio < 0.2 & tile_wc_ratio > -0.2] <- "wc"
  # to list all possiblities in table, need to convert to factor
  ranges_tiles$states <- factor(ranges_tiles$states, levels = c("wc", "ww", "cc", "?"))
  
  
  return(ranges_tiles)
}


background_estimation_on_bins <- function(bam_path, bin_width=10e+06, basename_pattern="^([^\\.]*\\.[^\\.]*).*"){

  ID <- sub(basename_pattern, "\\1", basename(bam_path))
  
  reads <- readBamFileAsGRanges(bam_path, pairedEndReads = TRUE, remove.duplicate.reads = FALSE)
  w_reads <- reads[strand(reads) == "+"]
  c_reads <- reads[strand(reads) == "-"]
  
  bins <- tileGenome(seqlengths = seqlengths(reads), tilewidth = bin_width, cut.last.tile.in.chrom=TRUE)
  
  bin_counts <- state_assignment_bins(bins, reads)
  
  ## Estimate background reads
  ww <- bin_counts[bin_counts$states == 'ww']
  cc <- bin_counts[bin_counts$states == 'cc']
  
  if (length(ww) > 0) {
    #add pseudocounts to avoid working with zeros
    ww$Cs <- ww$Cs+1
    ww$Ws <- ww$Ws+1
    bg_estim_ww <-  sum(ww$Cs) / sum(ww$Ws)
  } else {
    bg_estim_ww <- 0
  }
  
  if (length(cc) > 0) {
    #add pseudocounts to avoid working with zeros
    cc$Ws <- cc$Ws+1
    cc$Cs <- cc$Cs+1
    bg_estim_cc <- sum(cc$Ws) / sum(cc$Cs)
  } else {
    bg_estim_cc <- 0
  }
  bg_estim <- mean(c(bg_estim_ww, bg_estim_cc))
  
  
  ## Calculate reads per megabase
  seqlnth <- seqlengths(reads)[!is.na(seqlengths(reads))]
  if (length(seqlnth) > 0) {
    tiles <- unlist(GenomicRanges::tileGenome(seqlnth, tilewidth = 1000000))
    counts <- GenomicRanges::countOverlaps(tiles, reads)
    reads_per_MB <- round(median(counts))
  } else {
    reads_per_MB <- NA
  }
  

  out <- data.frame(ID = ID,
                    background.estimate=bg_estim,
                    med.reads.per.MB=reads_per_MB)

  return(out)

}




bkg_plot <- function(df, out_html){
  
  
  df <- df %>% mutate(library=sapply(strsplit(ID, "\\.|_"), "[[", 1))
  
  p <- df %>% 
    mutate(filter_col = ifelse(background.estimate > 0 & background.estimate <= 0.08, "F", "T")) %>%
    ggplot(aes(x = library, y = background.estimate)) +
    geom_quasirandom(aes(color = filter_col), size = 0.1, dodge.width=NULL) +
    geom_boxplot(width = 0.2, outlier.size = 0.1, outlier.colour = NA, fill = NA, color="black") +
    xlab(NULL) +
    ylab("Background estimate") +
    # scale_y_continuous(expand = expansion(mult = c(0.02, 0)),
    #                limits = c(0, 0.3),
    #                breaks = seq(0,0.3,0.1)) +
    theme_classic() +
    geom_hline(yintercept = 0.08, alpha=0.5) +
    scale_color_manual(values = c("F" = "darkgrey", "T" = "tomato")) +
    theme(strip.background=element_blank(),
          legend.position='none',
          axis.text.x = element_text(angle=45, hjust = 1))
  
  cols <- length(unique(df$library))
  
  tmp_png <- tempfile(fileext = ".png")
  ggsave(tmp_png, plot = p, width = cols, height = 4, dpi = 300)
  
  # read and base64 encode
  b64 <- base64enc::base64encode(tmp_png)
  
  html <- sprintf('<!doctype html>
  <html>
    <head><meta charset="utf-8"><title>Background estimate</title></head>
    <body>
      <h1>Background estimate</h1>
      <img src="data:image/png;base64,%s" alt="ggplot" style="max-width:100%%;height:auto;">
    </body>
  </html>', b64)
  
  writeLines(html, out_html)

  # interactive_p <- ggplotly(p)
  # 
  # saveWidget(interactive_p, out_html, selfcontained = TRUE)

}


#' Import BAM file into GRanges
#'
#' Modified version of the original \code{readBamFileAsGRanges} from breakpointR.
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param bam_path Bamfile with aligned reads.
#' @param bamindex Bam-index file with or without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @param remove.duplicate.reads If DUP flag present in BAM file, a logical indicating whether or not duplicate reads should be kept.
#' @param markDuplicates If DUP flag present in BAM file, a logical indicator whether to add mcol marking which reads are duplicates.
#' @param pair2frgm Set to \code{TRUE} if every paired-end read should be merged into a single fragment.
#' @param filtAlt Set to \code{TRUE} if you want to filter out alternative alignments defined in 'XA' tag.
#' @param useNames Retain read names.
#' @param extendSoftClipped Extend read coordinates by soft clipped length.
#'
#' @return A \code{\link{GRanges}} object.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments first last
#' @author David Porubsky, Aaron Taudt, Ashley Sanders
#' @export
#' @examples
#'## Get an example file
#'exampleFolder <- system.file("extdata", "example_bams", package="breakpointRdata")
#'exampleFile <- list.files(exampleFolder, full.names=TRUE)[1]
#'## Load the file
#'fragments <- readBamFileAsGRanges(exampleFile, pairedEndReads=FALSE, chromosomes='chr22')
#'
readBamFileAsGRanges <- function(bam_path, bamindex=bam_path, chromosomes=NULL, pairedEndReads=FALSE, min.mapq=10,
                                 remove.duplicate.reads=TRUE, markDuplicates=FALSE, pair2frgm=FALSE, filtAlt=FALSE,
                                 useNames=FALSE, extendSoftClipped=FALSE) {
  
  ## Check if bamindex exists
  bamindex.raw <- sub('\\.bai$', '', bamindex)
  bamindex <- paste0(bamindex.raw,'.bai')
  if (!file.exists(bamindex)) {
    bamindex.own <- Rsamtools::indexBam(bam_path)
    warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
    bamindex <- bamindex.own
  }
  
  file.header <- Rsamtools::scanBamHeader(bam_path)[[1]]
  chrom.lengths <- file.header$targets
  chroms.in.data <- names(chrom.lengths)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  
  chroms2use <- intersect(chromosomes, chroms.in.data)
  if (length(chroms2use)==0) {
    chrstring <- paste0(chromosomes, collapse=', ')
    stop('The specified chromosomes ', chrstring, ' do not exist in the data. Please try ', paste(paste0('chr',chromosomes), collapse=', '), ' instead.')
  }
  
  ## Issue warning for non-existent chromosomes
  diff <- setdiff(chromosomes, chroms.in.data)
  if (length(diff)>0) {
    diffs <- paste0(diff, collapse=', ')
    warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
  }
  
  ## Check if bam is truly paired ended in case pairedEndReads set to TRUE
  is.Paired <- Rsamtools::testPairedEndBam(file = bam_path, index = bamindex)
  if (pairedEndReads) {
    if (!is.Paired) {
      warning("You are trying to process single-ended BAM as paired-ended, Please set proper BAM directioanlity!!!")
    }
  } else {
    if (is.Paired) {
      warning("You are trying to process paired-ended BAM as single-ended, Please set proper BAM directioanlity!!!")
    }
  }
  
  ## Import the file into GRanges
  gr <- GenomicRanges::GRanges(seqnames=chroms2use, ranges=IRanges::IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
  if (!remove.duplicate.reads & !markDuplicates) {
    if (pairedEndReads) {
      if (filtAlt) {
        data.raw <- GenomicAlignments::readGAlignmentPairs(bam_path, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq'),
                                                           use.names = useNames)
      } else {
        data.raw <- GenomicAlignments::readGAlignmentPairs(bam_path, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq'),
                                                           use.names = useNames)
      }
    } else {
      if (filtAlt) {
        data.raw <- GenomicAlignments::readGAlignments(bam_path, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq'),
                                                       use.names = useNames)
      } else {
        data.raw <- GenomicAlignments::readGAlignments(bam_path, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq'),
                                                       use.names = useNames)
      }
    }
  } else {
    if (pairedEndReads) {
      if (filtAlt) {
        #NOTE: remove duplicates don't work if there is only one mate of the pair marked as duplicated. Suggestion force proper pairs to do it duplicate removal automatically!!!
        data.raw <- GenomicAlignments::readGAlignmentPairs(bam_path, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what=c('mapq', 'flag')),
                                                           use.names = useNames)
      } else {
        data.raw <- GenomicAlignments::readGAlignmentPairs(bam_path, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=c('mapq', 'flag')),
                                                           use.names = useNames)
      }
    } else {
      if (filtAlt) {
        data.raw <- GenomicAlignments::readGAlignments(bam_path, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=FALSE)),
                                                       use.names = useNames)
      } else {
        data.raw <- GenomicAlignments::readGAlignments(bam_path, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=FALSE)),
                                                       use.names = useNames)
      }
    }
  }
  
  if (pairedEndReads) {
    if (pair2frgm) {
      if (is.null(min.mapq)) { min.mapq <- 0 }
      #only proper pairs can be merged into single fragment (no negative ranges)
      data.prop.pairs <- data.raw[GenomicAlignments::isProperPair(data.raw)]
      
      data.first.ga <- GenomicAlignments::first(data.prop.pairs)
      data.last.ga <- GenomicAlignments::last(data.prop.pairs)
      
      data.first <- as(data.first.ga, 'GRanges')
      data.last <- as(data.last.ga, 'GRanges')
      
      if(extendSoftClipped){
        gr_cigar_first <- GenomicAlignments::cigar(data.first.ga)
        gr_cigar_last <- GenomicAlignments::cigar(data.last.ga)
        
        data.first <- extend_soft_clipped(data.first, gr_cigar_first)
        data.last <- extend_soft_clipped(data.last, gr_cigar_last)
      }
      
      ## filter XA tag
      if (filtAlt) {
        data.first.xa <- is.na(S4Vectors::mcols(data.first)$XA)
        data.last.xa <- is.na(S4Vectors::mcols(data.last)$XA)
        
        mask <- data.first.xa & data.last.xa
        
        data.first <- data.first[mask]
        data.last <- data.last[mask]
      }
      
      
      
      ## Filter by mapping quality and duplicate reads
      if (!is.null(min.mapq)) {
        if (any(is.na(S4Vectors::mcols(data.first)$mapq)) | any(is.na(S4Vectors::mcols(data.last)$mapq))) {
          warning(paste0(bam_path,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
          S4Vectors::mcols(data.first)$mapq[is.na(S4Vectors::mcols(data.first)$mapq)] <- -1
          S4Vectors::mcols(data.last)$mapq[is.na(S4Vectors::mcols(data.last)$mapq)] <- -1
        }
        
        data.first.filt <- S4Vectors::mcols(data.first)$mapq >= min.mapq
        data.last.filt <- S4Vectors::mcols(data.last)$mapq >= min.mapq
        
        mask <- data.first.filt & data.last.filt
        
        data.first <- data.first[mask]
        data.last <- data.last[mask]
        
        if (remove.duplicate.reads) {
          bit.flag <- bitwAnd(1024, S4Vectors::mcols(data.first)$flag)
          data.first.filt <- bit.flag == 0
          bit.flag <- bitwAnd(1024, S4Vectors::mcols(data.last)$flag)
          data.last.filt <- bit.flag == 0
          
          mask <- data.first.filt & data.last.filt
          data.first <- data.first[mask]
          data.last <- data.last[mask]
        }else if(markDuplicates){
          bit.flag <- bitwAnd(1024, S4Vectors::mcols(data.first)$flag)
          data.first.dup <- bit.flag == 1024
          bit.flag <- bitwAnd(1024, S4Vectors::mcols(data.last)$flag)
          data.last.dup <- bit.flag == 1024
          
          duplicates <- data.first.dup | data.last.dup
          S4Vectors::mcols(data.first)$duplicate <- FALSE
          S4Vectors::mcols(data.first)$duplicate[duplicates] <- TRUE
          S4Vectors::mcols(data.last)$duplicate <- FALSE
          S4Vectors::mcols(data.last)$duplicate[duplicates] <- TRUE
          
        }
      }
      
      #split reads by directionality
      data.first.plus <- data.first[strand(data.first) == '+']
      data.first.minus <- data.first[strand(data.first) == '-']
      data.last.plus <- data.last[strand(data.last) == '+']
      data.last.minus <- data.last[strand(data.last) == '-']
      
      #merge pairs into a single range
      frag.plus.mapq <- data.first.plus$mapq + data.last.minus$mapq
      frag.minus.mapq <- data.first.minus$mapq + data.last.plus$mapq
      
      if(markDuplicates){
        data.frag.plus <- GenomicRanges::GRanges(seqnames=seqnames(data.first.plus), ranges=IRanges::IRanges(start=start(data.first.plus), end=end(data.last.minus)), strand=strand(data.first.plus), mapq=frag.plus.mapq, duplicate=data.first.plus$duplicate)
        GenomeInfoDb::seqlengths(data.frag.plus) <- GenomeInfoDb::seqlengths(data.first)
        data.frag.minus <- GenomicRanges::GRanges(seqnames=seqnames(data.first.minus), ranges=IRanges::IRanges(start=start(data.last.plus), end=end(data.first.minus)), strand=strand(data.first.minus), mapq=frag.minus.mapq, duplicate=data.first.minus$duplicate)
        GenomeInfoDb::seqlengths(data.frag.minus) <- GenomeInfoDb::seqlengths(data.first)
        
      }else{
        data.frag.plus <- GenomicRanges::GRanges(seqnames=seqnames(data.first.plus), ranges=IRanges::IRanges(start=start(data.first.plus), end=end(data.last.minus)), strand=strand(data.first.plus), mapq=frag.plus.mapq)
        GenomeInfoDb::seqlengths(data.frag.plus) <- GenomeInfoDb::seqlengths(data.first)
        data.frag.minus <- GenomicRanges::GRanges(seqnames=seqnames(data.first.minus), ranges=IRanges::IRanges(start=start(data.last.plus), end=end(data.first.minus)), strand=strand(data.first.minus), mapq=frag.minus.mapq)
        GenomeInfoDb::seqlengths(data.frag.minus) <- GenomeInfoDb::seqlengths(data.first)
      }
      
      
      if(useNames){
        names(data.frag.plus) <- names(data.first.plus)
        names(data.frag.minus) <- names(data.first.minus)
      }
      
      data <- GenomicRanges::sort(c(data.frag.plus, data.frag.minus), ignore.strand=TRUE)
      
    } else {
      
      data.prop.pairs <- data.raw[GenomicAlignments::isProperPair(data.raw)]
      
      data.ga <- GenomicAlignments::first(data.prop.pairs)
      data <- as(data.ga, 'GRanges') #use only first mate of the read pair in subsequent analysis!!!
      
      if(extendSoftClipped){
        gr_cigar <- GenomicAlignments::cigar(data.ga)
        data <- extend_soft_clipped(data, gr_cigar)
      }
      
      ## Filter by mapping quality
      if (!is.null(min.mapq)) {
        if (any(is.na(S4Vectors::mcols(data)$mapq))) {
          warning(paste0(bam_path,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
          S4Vectors::mcols(data)$mapq[is.na(S4Vectors::mcols(data)$mapq)] <- -1
        }
        data <- data[S4Vectors::mcols(data)$mapq >= min.mapq]
      }
      
      ## Filter XA tag
      if (filtAlt) {
        data <- data[is.na(S4Vectors::mcols(data)$XA)]
      }
      ## Filter out duplicates
      if (remove.duplicate.reads) {
        bit.flag <- bitwAnd(1024, data$flag)
        mask <- bit.flag == 0
        data <- data[mask]
      }else if(markDuplicates){
        bit.flag <- bitwAnd(1024, data$flag)
        duplicates <- bit.flag == 1024
        S4Vectors::mcols(data)$duplicate <- FALSE
        S4Vectors::mcols(data)$duplicate[duplicates] <- TRUE
        
      }
    }
    
  } else {
    
    data <- as(data.raw, 'GRanges')
    
    if(extendSoftClipped){
      gr_cigar <- GenomicAlignments::cigar(data.raw)
      data <- extend_soft_clipped(data, gr_cigar)
    }
    
    ## Filter by mapping quality
    if (!is.null(min.mapq)) {
      if (any(is.na(S4Vectors::mcols(data)$mapq))) {
        warning(paste0(bam_path,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
        S4Vectors::mcols(data)$mapq[is.na(S4Vectors::mcols(data)$mapq)] <- -1
      }
      data <- data[S4Vectors::mcols(data)$mapq >= min.mapq]
    }
    ## Filter XA tag
    if (filtAlt) {
      data <- data[is.na(S4Vectors::mcols(data)$XA)]
    }
  }
  GenomeInfoDb::seqlevels(data) <- GenomeInfoDb::seqlevels(gr)
  return(data)
}




main()
