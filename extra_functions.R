library(soGGi)
library(GenomicRanges)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)

#' Cut off data at upper and lower quantiles
#' 
#' Replaces all data points beyond a quantile with the value at that quantile, 
#' see: https://en.wikipedia.org/wiki/Winsorising. 
#' Use for removing outliers before plotting data. This function works on vectors 
#' or matrices as it's based on subsetting with `[]`. 
#' 
#' @param vec Vector or matrix to apply winsorisation to.
#' @param lower Lower quantile for winsorisation. Numeric fractional quantile,
#' should be between zero and 1. 
#' @param upper Upper quantile for winsorisation. Numeric fractional quantile,
#' should be between zero and 1. 
#' 
#' @value Winsorised vector or matrix.

winsorise_vec <- function(vec, lower = 0.01, upper = 0.99){
  if (lower < 0 | upper > 1) {
    stop("Lower and upper bounds for winsorisation must be between 0 and 1")
  }
  if (lower > upper) {
    stop("Lower bound is higher than upper bound!")
  } 
  q_l <- quantile(vec, probs = lower)
  q_u <- quantile(vec, probs = upper)
  vec[vec < q_l] <- q_l
  vec[vec > q_u] <- q_u
  return(vec)
}

#' Reshape a ChIPprofile object into a data.frame
#' 
#' This function takes a soGGi ChIPprofile object and reshapes it into a 
#' long-format data.frame. Metadata columns from the rowRanges are added to each
#' assay before assays are combined into a single data.frame with an "assay" 
#' column that uses the names of the assays to identify the data. 
#' 
#' @param object A ChIPprofile object
#' 
#' @value a long-format tbl_df of the genomic data.

reshape_chipprofile <- function(object){
  require(dplyr)
  require(tidyr)
  
  gr <- rowRanges(object)
  obj_assays <- assays(object)
  assay_names <- metadata(object)$names
  if (is.null(assay_names)){
    assay_names <- 1:length(obj_assays)
    }
  names(obj_assays) <- assay_names
  
  obj_reshape_list <- lapply(assay_names, function(n){
    df <- as.data.frame(cbind(mcols(gr), assay = n,
                              as.data.frame(obj_assays[[n]])))
    #get columns with data at each position - reshape these to long format
    g_cols <- names(df)[!(names(df) %in% c(names(mcols(gr)), "assay"))]
    df_reshape <- gather_(df, "xIndex", "score", g_cols)
  })
  
  obj_reshape <- do.call("rbind", obj_reshape_list)
  return(tbl_df(obj_reshape))
}

#' Function to plot a heatmap from a data.frame of genomic data
#' 
#' This function takes a long-format data.frame containing data extracted from a
#' ChIPprofile object and plots it as a heatmap using geom_raster.
#' 
#' @param obj A data.frame or tbl_df of long-form genomic data.
#' @param x Column name to use for x position, default "xIndex"
#' @param y Column name to use for y position, default "giID"
#' @param fill Column name to use for fill of heatmap, default "score
#' @param winsorise Vector of lower and upper bounds for winsorisation using 
#' winsorise_vec, should be numeric values between zero and one. Default = NULL 
#' for no winsorisation.
#' 
#' @value a gpplot plot object


gg_heatmap <- function(obj, x = "xIndex", y = "giID", fill = "score", 
                       winsorise = NULL){
  
  #winsorise score column
  if(!is.null(winsorise)){
    if (length(winsorise)!=2){
      stop("Must supply lower and upper bounds for winsorisation")}
  obj[[fill]] <- winsorise_vec(obj[[fill]], 
                             lower = winsorise[1], 
                             upper = winsorise[2])
  }
  
  p <- ggplot(obj, aes_string(x = x, y = y, fill = fill))
  
  p <- p + geom_raster() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  return(p)
  
}

#' Summarises genomic signal in a data.frame
#' 
#' This function takes a long form data frame of genomic signal and summarises 
#' it by position and by any additional columns you specify.
#' 
#' @param obj_reshape Long-form data.frame of genomic data
#' @param summariseBy Column names to summarise by. Default NULL to only summarise 
#' by position (column "xIndex").
#' @param summary String specifying function to summarise data, e.g. mean, median, 
#' sum. Default "mean". 
#' @param range Method to use to summarise the range of the data. Default NULL. 
#' Currently only quantile ranges are implemented, e.g. to return the upper and 
#' lower quartiles of the data use `c(0.25, 0.75)`. Any percentiles can be used, 
#' e.g. `c(0.01, 0.05, 0.95, 0.99)`. 
#' TO DO: implement confidence intervals, sd, etc...
#' 
#' @value a data.frame of summarised data

summarise_signal <- function(obj_reshape, summariseBy = NULL, summary = "mean",
                             range = NULL){
  #get summary function e.g. mean, median, sum
  #any checks needed here?
  summary_func <- lazyeval::interp(~f(score), f = as.name(summary))
  
  ## ranges
  quantile_list <- NULL
  if (is.numeric(range) & all(range >= 0 & range <= 1)){
    #and length==2?
    message("Treating 'range' argument as if specifying quantiles of data to return...")
    quantile_list <- lapply(range, function(q){
      lazyeval::lazy(quantile(score, probs = q))
    })
    names(quantile_list) <- paste0("q_", range)
  }
  
  args_list <- c(score_summary = summary_func, quantile_list)
  
  #group by position and any other metadata
  grouping <- c(summariseBy, "xIndex")
  obj_reshape <- group_by_(obj_reshape, .dots = grouping)
  
  #summarise
  obj_summary <- summarise_(obj_reshape, 
                            .dots = args_list)
  #rename columns
  obj_summary <- rename_(obj_summary,
                         .dots = setNames("score_summary", paste(summary, "score", sep = "_")))
  return(obj_summary)
}

