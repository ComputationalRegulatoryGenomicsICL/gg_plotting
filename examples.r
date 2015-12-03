##Examples

load("test_plot_image.RData")

test <- object
test_r <- reshape_chipprofile(test)
gg_heatmap(test_r, winsorise = c(0, 0.99), facet = "class")
gg_heatmap(test_r, winsorise = c(0, 0.99), facetBy_grid = c("class", "assay")) + 
  scale_fill_viridis(option = "magma") 

test2 <- c(object, object)
metadata(test2)$names[2] <- "assay2"
test2<- reshape_chipprofile(test2)

gg_heatmap(test2, winsorise = c(0, 0.99), facetBy_grid = c("assay", "class")) + 
    scale_fill_viridis(option = "magma")
gg_heatmap(test2, winsorise = c(0, 0.99), facetBy_grid = c("class", "assay")) + 
    scale_fill_viridis(option = "magma")

gg_heatmap(test2, winsorise = c(0, 0.99)) + 
  scale_fill_viridis(option = "magma")

gg_heatmap(test2, winsorise = c(0, 0.99)) + 
  #scale_fill_viridis(option = "magma")+
  facet_wrap(class ~ assay, scales = "free", nrow = 4) +
  scale_x_discrete("Position", breaks = c("Start.1.1", "End.1"),
                   labels = c("Start.1.1" = "Start", "End.1" = "End"))

grouping <- c("class")
test_avg <- test_r %>% 
  group_by_(.dots = c(grouping, "xIndex")) %>% 
  summarise(mean_score = mean(score), 
            median_score = median(score),
            uq = quantile(score, probs = 0.75), 
            lq = quantile(score, probs = 0.25)) 

ggplot(test_avg, aes(x = xIndex, y = median_score, group = 1)) + 
  geom_ribbon(aes(ymax = uq, ymin = lq), alpha = 0.5, colour = "grey") +
  geom_line() +
  theme_bw()

+ 
  facet_wrap(~class)

test_avg <- summarise_signal(test_r, summariseBy = "class", summary = "median", 
                             range = c(0.25, 0.75))
test_avg$xIndex <- 1:300

ggplot(test_avg, aes(x = xIndex, y = median_score, colour = class, group = 1)) + 
  theme_bw() + 
  facet_wrap(~class) +
  geom_ribbon(aes(ymin = q_0.25, ymax = q_0.75), colour = NA, fill = "grey", alpha = 0.5) + 
  geom_line() +
  coord_cartesian(ylim = c(0, 15))


test_r$xIndex <- rep(1:300, each = 542)
test_r$score <- winsorise_vec(test_r$score, lower = 0, upper = 0.99)
ggplot(test_r, aes(x = xIndex, y = giID, alpha = score)) + 
  geom_point(size = 1, colour = "darkblue") + 
  scale_alpha(range = c(0,1)) + theme_bw()

big_test <- tbl_df(data.frame())
for (i in 1:10){
  tmp <- test_r
  tmp$giID <- paste(tmp$giID, sep = "_", i)
  big_test <- rbind(big_test, tmp)
}

ggplot(big_test, aes(x = xIndex, y = giID, alpha = score)) + 
  geom_point(size = 1, colour = "darkblue") + 
  scale_alpha(range = c(0,1)) + theme_bw()

gg_heatmap(big_test, winsorise = c(0, 0.99))
+ 
  #scale_fill_viridis(option = "magma")+
  #facet_wrap(class ~ assay, scales = "free", nrow = 4) +
  scale_x_discrete("Position", breaks = c("Start.1.1", "End.1"),
                   labels = c("Start.1.1" = "Start", "End.1" = "End"))

big_test_no_zero <- big_test[-which(big_test$score ==0),]
gg_heatmap(big_test_no_zero, winsorise = c(0, 0.99)) + 
  scale_fill_continuous(low = "white", high = "darkblue", na.value = "white")
