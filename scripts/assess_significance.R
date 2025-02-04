#!/usr/bin/env Rscript

# Parse Args --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3 | length(args) > 3) {
  stop(
    "At least three arguments must be supplied:
    assess_significance.R [input_CNVs] [input_ratio.txt] [output_file]",
    call. = FALSE
  )
}

cnv_path <- args[1]
ratio_path <- args[2]
output_path <- args[3]

# Load Dependencies -------------------------------------------------------
suppressMessages({library(rtracklayer)})

if (require(matrixTests)) {
  message("[INFO] 'matrixTests' package available, using faster function 'matrixTests::col_wilcoxon_twosample' for wilcoxon test.")
  wilcoxon_stat_function <- matrixTests::col_wilcoxon_twosample
} else {
  message("[INFO] 'matrixTests' package available, using slower function 'stats::wilcox.test' for Wilcoxon test.")
  wilcoxon_stat_function <- stats::wilcox.test
}

# Functions ---------------------------------------------------------------
#' Get Expected CNV Column Names
#' @description Deduces the expected column names from `cnv_path` based on the
#' number of columns it has.
#' @param cnv_path path to CNV file.
#'
#' @returns vector of expected column names
get_cnv_colnames <- function(cnv_path) {
  colname_groups = list(
    base = c("chr", "start", "end", "copy number", "status"),
    genotyped = c("genotype", "uncertainty"),
    paired = c("somatic/germline", "percentageOfGermline")
  )
  cnv_col_no <- length(colnames(read.table(cnv_path, header=FALSE)))
  
  if (cnv_col_no == 5) {
    col_names <- c(colname_groups[["base"]])
  } else if (cnv_col_no == 7) {
    col_names <- c(colname_groups[["base"]], colname_groups[["genotyped"]])
  } else if (cnv_col_no == 9) {
    col_names <- c(
      colname_groups[["base"]],
      colname_groups[["genotyped"]],
      colname_groups[["paired"]]
    )
  }
  return(col_names)
}

#' Calculates p-value of cnv_ratios against normal_ratios given
#' @description Calculates the p-value for either Wilcoxon or Kolmogorov-Smirnov
#' tests, returning `NA` if there are any raised errors.
#' @param cnv_ratios the ratios for genomic ranges that intersect with CNVs
#' @param normal_ratios the ratios for genomic ranges that do not intersect CNVs
#' @param stat_function the stat function to use
#' (`matrixTests::col_wilcoxon_twosample|stats::wilcox.test|stats::ks.test`)
#'
#' @returns p-value
get_pval <- function(cnv_ratios, normal_ratios, stat_function) {
  tryCatch(
    {
      result <- suppressWarnings({stat_function(cnv_ratios,normal_ratios)})
      names(result)[names(result) == "p.value"] <- "pvalue"
      pval <- result$pvalue
      return(pval)
    },
    error = function(e) {
      cat("An error occurred, returning 'NA':", conditionMessage(e), "\n")
      return(NA)
    }
  )
}

#' Main body function
#'
#' @param cnv_path path to `*_CNVs` file produce by ControlFREEC
#' @param ratio_path path to `*_ratio.txt` file produce by ControlFREEC 
#' @param output_path output path for annotated CNV file
#' @param wilcoxon_stat_function stat function to use for Wilcoxon test. If
#' `MatrixTests` is installed, use faster `matrixTests::col_wilcoxon_twosample`
#'
#' @returns Writes annotated CNV `DataFrame` with p-values to `output_path`
main <- function(cnv_path, ratio_path, output_path, wilcoxon_stat_function) {
  message("[INFO] Loading data")
  ratio_df <- read.table(ratio_path, header = TRUE)
  ratio_df$Ratio[which(ratio_df$Ratio == -1)] = NA
  cnv_df <- read.table(cnv_path, col.names = get_cnv_colnames(cnv_path))
  
  message("[INFO] Retrieving genomic ranges")
  cnv_bed_gr <- GenomicRanges::GRanges(
    seqnames = cnv_df$chr,
    ranges = IRanges::IRanges(cnv_df$start, cnv_df$end)
  )
  ratio_bed_gr <- GenomicRanges::GRanges(
    seqnames = ratio_df$Chromosome,
    ranges = IRanges::IRanges(ratio_df$Start, ratio_df$Start),
    score = ratio_df$Ratio
  )
  
  normals_gr <- IRanges::subsetByOverlaps(
    ratio_bed_gr,
    IRanges::setdiff(ratio_bed_gr, cnv_bed_gr)
  )
  normal_ratios <- BiocGenerics::score(normals_gr)
  
  message("[INFO] Calculating P-values:")
  pb <- utils::txtProgressBar(
    min = 0,
    max = nrow(cnv_df),
    initial = 0,
    style = 3
  )
  pvalues <- list(wilcox = c(), ks = c())
  
  for (i in 1:nrow(cnv_df)) {
    utils::setTxtProgressBar(pb, i)
    
    cnv_ratios <- BiocGenerics::score(
      IRanges::subsetByOverlaps(ratio_bed_gr, cnv_bed_gr[i])
    )
    
    wilcox_pval <- get_pval(cnv_ratios, normal_ratios, wilcoxon_stat_function)
    ks_pval <- get_pval(cnv_ratios, normal_ratios, stats::ks.test)
    
    pvalues[["wilcox"]] <- c(pvalues[["wilcox"]], wilcox_pval)
    pvalues[["ks"]] <- c(pvalues[["ks"]], ks_pval)
  }
  close(pb)
  
  message("[INFO] Appending P-values to original CNV data")
  cnv_df[["WilcoxonRankSumTestPvalue"]] <- pvalues[["wilcox"]]
  cnv_df[["KolmogorovSmirnovPvalue"]] <- pvalues[["ks"]]
  
  message(paste("[INFO] Writing output to:", output_path))
  write.table(
    cnv_df,
    file = output_path,
    sep="\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# Main --------------------------------------------------------------------
main(cnv_path, ratio_path, output_path, wilcoxon_stat_function)
