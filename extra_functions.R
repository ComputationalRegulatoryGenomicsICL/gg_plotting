library(soGGi)
library(GenomicRanges)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)

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

