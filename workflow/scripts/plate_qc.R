

library(argparse)
library(dplyr)
library(ggplot2)
library(base64enc)
library(patchwork)



parse_arguments <- function(){
  parser <- ArgumentParser(description="Run Tn5 and ligation plate QC")
  parser$add_argument("-i","--input", type="character", dest="ss_summary",
                      required=TRUE, nargs="+",
                      help="Path to single_cell_summary_all.txt produced by sciL3Pipe")
  parser$add_argument("-o", "--output_html", dest="output_html", type="character",
                      required=TRUE,
                      help="Output html location")
  args <- parser$parse_args()

  return(args)
}



main <- function(){

  args <- parse_arguments()

  plot_plate_qc(args$ss_summary, out_html = args$output_html)

}


plot_plate_qc <- function(ss_summary, out_html, samples = NULL){

  if(is(ss_summary, "character")){
    ss_summary <- readr::read_table(ss_summary, col_names = FALSE) %>%
      setNames(c("Cell_ID", "Mapped_reads", "Mapped_reads_MAPQgt0", "Unique_reads_by_R1_start_and_end", "ratio")) %>%
      distinct()

  }else if(!is(ss_summary, "data.frame")){
    stop("Incorrect input format for ss_summary.")
  }

  if(!all(names(ss_summary) %in% c("Cell_ID", "Mapped_reads", "Mapped_reads_MAPQgt0", "Unique_reads_by_R1_start_and_end", "ratio"))){
    names(ss_summary) <- c("Cell_ID", "Mapped_reads", "Mapped_reads_MAPQgt0", "Unique_reads_by_R1_start_and_end", "ratio")
  }



  if(!is.null(samples)){
    ss_summary <- ss_summary %>% filter(Cell_ID %in% samples)
    if(nrow(ss_summary) == 0){
      stop("No specified samples in ss_summary.")
    }
  }


  ss_summary <- ss_summary %>% mutate(Tn5=substring(sapply(strsplit(Cell_ID, "\\."), "[[", 2),1,8)) %>%
    mutate(ligation=substring(sapply(strsplit(Cell_ID, "\\."), "[[", 2),9,15)) %>%
    mutate(sss=sapply(strsplit(Cell_ID, "\\.|_"), "[[", 2)) %>%
    mutate(library=sapply(strsplit(Cell_ID, "\\.|_"), "[[", 1)) %>%
    mutate(library_sss=sapply(strsplit(Cell_ID, "\\."), "[[", 1))



  tn5_barcodes <- data.frame(
    stringsAsFactors = FALSE,
    Tn5 = c("TGATATTG", "GATCCCGT","CTCGATTA",
            "CATCAAGG","TCCTTGTG","GGTCATAT","ATCGCGTT",
            "CATGCCCC","GTTACGCG","CCGCGCTT","TCTTAGTG","TCGGCCTA",
            "CTTTCTCT","TCGCGTTT","GTCAGTAG","CCATGGAA",
            "ATGCTGCG","GAGTCTTT","TACGATAT","ACCATTTA","ATCGGGAC",
            "GACGTCGG","CATTGTGT","TTTGACTC")) %>%
    mutate(barcode_num=seq(1,24,1))


  ligation_barcodes <- data.frame(stringsAsFactors = FALSE,
                                  ligation = c("GCCACCT", "GCTATTA", "TATGTTG", "TTAACCG", "GGGTACA", "TTCTATT", "GCTTGAT", "TCCAAGT",
                                               "AGAACTA", "TCGGTTT", "AGAGACT", "CTGTTAA", "AGGTAGT", "GGCTTGG", "CACTGTT", "ACGTCGT",
                                               "TGCTTAA", "TCCATAG", "GTGCCGC", "GCAGGTC", "TGCACCG", "AGAGACT", "ATAAAAG", "GGTCCCA",
                                               "GTCGCAC", "GTAAAGG", "CGACTTG", "CCGCTTA", "ATGGTCA", "CCATCCA", "GGGCGAG", "CTGCATG",
                                               "CCTACAG", "CAGAGGT", "CAAAACG", "GACCTTC", "AGTAGCC", "TAGCCGG", "TAGTCTA", "GTAATTT",
                                               "TCGAGTA", "GGTAGGA", "GACGGGG", "TCGTATC", "TCACAGC", "GCCTATA", "TGCGTCG", "AAATGGA",
                                               "CGGTCTT", "TTACTTA", "GCAGTAG", "TATAAGA", "GTGGGTT", "ATCCGTA", "TGGAATT", "CGGAGAC",
                                               "ACTGCAT", "CAAGCTC", "GTTTCTC", "GGCCAAA", "ACTCGCA", "TGATGCA", "ATATCCC", "GAGGCGA",
                                               "CCTCACT", "AACTGAC", "AGTCTTC", "TAAACTT", "ATACCAT", "CATTGCG", "ACCGCTG", "GATTGTG")) %>%
    mutate(ligation_num=seq(1,72,1))


  ss_summary <- ss_summary %>% left_join(tn5_barcodes, relationship = "many-to-many") %>% left_join(ligation_barcodes, relationship = "many-to-many")

  indx <- expand.grid(seq(1, 12), seq(1, 8))

  plate <- data.frame(well_id = c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12",
                                  "B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12",
                                  "C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12",
                                  "D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12",
                                  "E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12",
                                  "F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12",
                                  "G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12",
                                  "H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12"),
                      Tn5 = c("TGATATTG", "GATCCCGT", "CTCGATTA", "CATCAAGG", "TCCTTGTG", "GGTCATAT", "ATCGCGTT", "CATGCCCC", "GTTACGCG", "CCGCGCTT", "TCTTAGTG", "TCGGCCTA",
                              rep(NA,12),
                              rep(NA,12),
                              "CTTTCTCT", "TCGCGTTT","GTCAGTAG", "CCATGGAA", "ATGCTGCG", "GAGTCTTT", "TACGATAT", "ACCATTTA", "ATCGGGAC", "GACGTCGG", "CATTGTGT", "TTTGACTC",
                              rep(NA,12),
                              rep(NA,12),
                              rep(NA,12),
                              rep(NA,12)),
                      barcode_num = c(seq(1,12,1),
                                      rep(NA,12),
                                      rep(NA,12),
                                      seq(13,24,1),
                                      rep(NA,12),
                                      rep(NA,12),
                                      rep(NA,12),
                                      rep(NA,12))
  ) %>%
    mutate(column = indx$Var1, row = indx$Var2)


  tn5_summary <- ss_summary %>% group_by(Tn5) %>%
    summarise(uniq_r1_start_end = sum(Unique_reads_by_R1_start_and_end),
              mapped_qgt0 = sum(Mapped_reads_MAPQgt0),
              num_cells = n(),
              barcode_num=dplyr::first(barcode_num)) %>%
    mutate(reads_per_cell = uniq_r1_start_end/num_cells,
           ratio = mapped_qgt0/uniq_r1_start_end)


  plate_tn5 <- plate %>% left_join(tn5_summary)

  plate_tn5 <- plate_tn5 %>% mutate(column = as.numeric(column), row = as.numeric(row))


  p_t <-  ggplot(data=plate_tn5, aes(x=column, y=row)) +
    geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2),
               color="grey90", shape=21) +
    scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
    # coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.5, 12.5), ylim=c(0.5, 8.5)) +
    scale_x_continuous(breaks=seq(1, 12)) +
    ylab(NULL) + xlab(NULL)

  print("here")

  p1 <- p_t + geom_point(aes(size = reads_per_cell, colour = ratio)) + labs(title="Tn5 Plate - ratio")
  p2 <- p_t + geom_point(aes(size = num_cells, colour = reads_per_cell)) + labs(title="Tn5 Plate - reads per cell")



  print("here2")


  ligation_summary <- ss_summary %>% group_by(ligation) %>%
    summarise(uniq_r1_start_end = sum(Unique_reads_by_R1_start_and_end),
              mapped_qgt0 = sum(Mapped_reads_MAPQgt0),
              num_cells = n(),
              ligation_num=dplyr::first(ligation_num)) %>%
    mutate(reads_per_cell = uniq_r1_start_end/num_cells,
           ratio = mapped_qgt0/uniq_r1_start_end)


  indx_lig <- expand.grid(rev(seq(1, 8)), seq(1, 12))


  plate_ligation_barcode <- data.frame(well_id = c("H1","G1","F1","E1","D1","C1","B1","A1",
                                                   "H2","G2","F2","E2","D2","C2","B2","A2",
                                                   "H3","G3","F3","E3","D3","C3","B3","A3",
                                                   "H4","G4","F4","E4","D4","C4","B4","A4",
                                                   "H5","G5","F5","E5","D5","C5","B5","A5",
                                                   "H6","G6","F6","E6","D6","C6","B6","A6",
                                                   "H7","G7","F7","E7","D7","C7","B7","A7",
                                                   "H8","G8","F8","E8","D8","C8","B8","A8",
                                                   "H9","G9","F9","E9","D9","C9","B9","A9",
                                                   "H10","G10","F10","E10","D10","C10","B10","A10",
                                                   "H11","G11","F11","E11","D11","C11","B11","A11",
                                                   "H12","G12","F12","E12","D12","C12","B12","A12"),
                                       ligation = c("GCCACCT", "GCTATTA", "TATGTTG", "TTAACCG", "GGGTACA", "TTCTATT", "GCTTGAT", "TCCAAGT",
                                                    "AGAACTA", "TCGGTTT", "AGAGACT", "CTGTTAA", "AGGTAGT", "GGCTTGG", "CACTGTT", "ACGTCGT",
                                                    "TGCTTAA", "TCCATAG", "GTGCCGC", "GCAGGTC", "TGCACCG", "AGAGACT", "ATAAAAG", "GGTCCCA",
                                                    "GTCGCAC", "GTAAAGG", "CGACTTG", "CCGCTTA", "ATGGTCA", "CCATCCA", "GGGCGAG", "CTGCATG",
                                                    "CCTACAG", "CAGAGGT", "CAAAACG", "GACCTTC", "AGTAGCC", "TAGCCGG", "TAGTCTA", "GTAATTT",
                                                    "TCGAGTA", "GGTAGGA", "GACGGGG", "TCGTATC", "TCACAGC", "GCCTATA", "TGCGTCG", "AAATGGA",
                                                    "CGGTCTT", "TTACTTA", "GCAGTAG", "TATAAGA", "GTGGGTT", "ATCCGTA", "TGGAATT", "CGGAGAC",
                                                    "ACTGCAT", "CAAGCTC", "GTTTCTC", "GGCCAAA", "ACTCGCA", "TGATGCA", "ATATCCC", "GAGGCGA",
                                                    "CCTCACT", "AACTGAC", "AGTCTTC", "TAAACTT", "ATACCAT", "CATTGCG", "ACCGCTG", "GATTGTG",
                                                    rep(NA,8),
                                                    rep(NA,8),
                                                    rep(NA,8)),
                                       barcode_num = c(seq(1,72,1),
                                                       rep(NA,8),
                                                       rep(NA,8),
                                                       rep(NA,8))) %>%
    mutate(row = indx_lig$Var1, column = indx_lig$Var2)

  plate_ligation <- plate_ligation_barcode %>% left_join(ligation_summary)

  plate_ligation <- plate_ligation %>% mutate(column = as.numeric(column), row = as.numeric(row))


  p_l <-  ggplot(data=plate_ligation, aes(x=column, y=row)) +
    geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2),
               color="grey90", shape=21) +
    scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
    # coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.5, 12.5), ylim=c(0.5, 8.5)) +
    scale_x_continuous(breaks=seq(1, 12)) +
    ylab(NULL) + xlab(NULL)




  p3 <- p_l + geom_point(aes(size = reads_per_cell, colour = ratio)) + labs(title="Ligation Plate - ratio")
  p4 <- p_l + geom_point(aes(size = num_cells, colour = reads_per_cell)) + labs(title="Ligation Plate - reads per cell")


  p <- p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2)


  tmp_png <- tempfile(fileext = ".png")
  ggsave(tmp_png, plot = p, width = 14, height = 8, dpi = 300)

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

}




main()
