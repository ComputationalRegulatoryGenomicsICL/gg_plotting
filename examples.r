##Examples

load("test_plot_image.RData")
test <- object

test_r <- reshape_chipprofile(test)

#basic plot
gg_heatmap(test_r, winsorise = c(0, 0.99)) 
#add faceting - need "free_y" for it to look nice, so only facet_wrap works not facet_grid
gg_heatmap(test_r, winsorise = c(0, 0.99)) + 
  facet_wrap(~class, scales = "free_y")
#more customisation
gg_heatmap(test_r, winsorise = c(0, 0.99)) + 
  facet_wrap(~class + assay, scales = "free_y") + 
  scale_x_discrete("Position", breaks = c("Start.1.1", "End.1"),
                    labels = c("Start.1.1" = "Start", "End.1" = "End")) # I want to improve labelling...

## with multiple assays
test2 <- c(object, object)
metadata(test2)$names[2] <- "assay2"
test2<- reshape_chipprofile(test2)
# facet by assay type - need to manually control column/row number of layout
gg_heatmap(test2, winsorise = c(0, 0.99)) + 
  #scale_fill_viridis(option = "magma") +
  facet_wrap(~class + assay, scales = "free_y", nrow = 4)
   

## summarising data
test_avg <- summarise_signal(test_r, summariseBy = "class", summary = "median", 
                             range = c(0.25, 0.75))


ggplot(test_avg, aes(x = xIndex, y = median_score, colour = class, group = 1)) + 
  theme_bw() + 
  facet_wrap(~class) +
  geom_ribbon(aes(ymin = q_0.25, ymax = q_0.75), colour = NA, fill = "grey", alpha = 0.5) + 
  geom_line() +
  coord_cartesian(ylim = c(0, 15))

## manually summarising data example...
## if you want something other than quantiles
test_avg <- test_r %>% 
  group_by(class, xIndex) %>% 
  summarise(mean_score = mean(score), 
            median_score = median(score),
            uq = quantile(score, probs = 0.75), 
            lq = quantile(score, probs = 0.25),
            sd = sd(score)) 

test_avg$xIndex <- 1:300 #hack to get x axis to look nicer, will fix!

ggplot(test_avg, aes(x = xIndex, y = median_score, group = 1)) + 
  geom_ribbon(aes(ymax = uq, ymin = lq), alpha = 0.5, colour = "grey") +
  geom_line() +
  theme_bw() + 
  facet_wrap(~class, scales = "free_y")

ggplot(test_avg, aes(x = xIndex, y = mean_score, group = 1)) + 
  geom_ribbon(aes(ymax = mean_score + sd, ymin = mean_score-sd), alpha = 0.5, colour = "grey") +
  geom_line() +
  theme_bw() +
  facet_wrap(~class, scales = "free_y")

# plotting as a scatter plot instead
# point size relative to plot size needs to be adjusted
test_r$xIndex <- rep(1:300, each = 542) # hack for pretty x axis again
test_r$score <- winsorise_vec(test_r$score, lower = 0, upper = 0.99)
ggplot(test_r, aes(x = xIndex, y = giID, alpha = score)) + 
  geom_point(size = 1, colour = "darkblue") + 
  scale_alpha(range = c(0,1)) + theme_bw()

## larger dataframe to test performance
big_test <- tbl_df(data.frame())
for (i in 1:10){
  tmp <- test_r
  tmp$giID <- paste(tmp$giID, sep = "_", i)
  big_test <- rbind(big_test, tmp)
}

gg_heatmap(big_test, winsorise = c(0, 0.99)) + 
  #scale_fill_viridis(option = "magma")+
  #facet_wrap(class ~ assay, scales = "free", nrow = 4) +
  scale_x_discrete("Position", breaks = c("Start.1.1", "End.1"),
                   labels = c("Start.1.1" = "Start", "End.1" = "End"))

#remove zeros to check size improvements
# note that you then need to set na.value the same as the zero value when plotting
big_test_no_zero <- big_test[-which(big_test$score ==0),]
gg_heatmap(big_test_no_zero, winsorise = c(0, 0.99)) + 
  scale_fill_continuous(low = "white", high = "darkblue", na.value = "white")
