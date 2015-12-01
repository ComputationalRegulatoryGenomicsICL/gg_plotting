library(soGGi)
library(GenomicRanges)
library(tidyr)
library(dplyr)
library(ggplot2)

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
  obj_reshape <- as.data.frame(cbind(mcols(rowRanges(object)), 
                                     as.data.frame(assays(object)[[1]])))
  
  #get columns with data at each position - reshape these to long format
  g_cols <- names(obj_reshape)[!(names(obj_reshape) %in% names(mcols(rowRanges(object))))]
  obj_reshape <- gather_(obj_reshape, "xIndex", "score", g_cols)
  return(tbl_df(obj_reshape))
}

gg_heatmap <- function(obj, x = "xIndex", y = "giID", fill = "score", 
                       facetBy = NULL, winsorise = NULL){
  
  #winsorise score column
  if(!is.null(winsorise)){
    if (length(winsorise)!=2){
      stop("Must supply lower and upper bounds for winsorisation")}
  obj[[fill]] <- winsorise_vec(obj[[fill]], 
                             lower = winsorise[1], 
                             upper = winsorise[2])
  }
  
  p <- ggplot(obj, aes_string(x = x, y = y, fill = fill))
  
  if (!is.null(facetBy)){
    p <- p + facet_wrap(formula(paste("~",paste(facetBy,collapse="+"))),
                        scales = "free")
  }
  
  p <- p + geom_raster() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  return(p)
  
}

plot_by_range_group <- function(object, colourBy=NULL, facetBy=NULL){
  require(tidyr)
  require(dplyr)
  obj_reshape <- as.data.frame(cbind(mcols(rowRanges(object)), 
                                     as.data.frame(assays(object)[[1]])))
  
  #get columns with data at each position - reshape these to long format
  g_cols <- names(obj_reshape)[!(names(obj_reshape) %in% names(mcols(rowRanges(object))))]
  obj_reshape <- gather_(obj_reshape, "xIndex", "score", g_cols)
  
  #group and calculate mean per position per group of ranges
  grouping <- list(colourBy, facetBy)
  grouping <- grouping[sapply(grouping, length) > 0]
  
  obj_reshape <- obj_reshape %>% 
    group_by_(.dots = c(grouping, "xIndex")) %>% 
    summarise(mean_score = mean(score)) 
  
  #get obj_reshape axis index numbers and join
  axisIndex <- c(seq(1,(object@params$distanceAround+object@params$distanceAround+1)))
  axisIndex_df <- data.frame(xIndex = unique(obj_reshape$xIndex), axisIndex = axisIndex)
  obj_reshape <- left_join(obj_reshape, axisIndex_df)
  
  #plot
  
  p <- ggplot(obj_reshape, aes_string(x="axisIndex",y="mean_score"))+
    geom_line(alpha = 1, size = 1)+
    xlim(0, 3001) +
    ylab("Score") +
    theme(axis.title.y=element_text(angle=0))
  
  p <- p + scale_x_continuous(breaks = c(1, object@params$distanceAround+1, object@params$distanceAround+1+object@params$distanceAround),
                              labels = c(paste0("-",object@params$distanceAround), "Centre", paste0("+", object@params$distanceAround))) +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
  
  p <- p + aes_string(colour = colourBy) 
  p <- p + facet_wrap(formula(paste("~",paste(facetBy,collapse="+"))))
  return(p)
}
 


#### 
# test <- all_se_signal_list[[5]]
# test_r <- reshape_chipprofile(test)
# gg_heatmap(test_r, winsorise = c(0, 0.99), facet = "class")
# gg_heatmap(test_r, winsorise = c(0, 0.99), facet = "class") + 
#   scale_fill_viridis(option = "magma") 
