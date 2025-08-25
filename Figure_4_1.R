####### Figure 4 time series plot
rm(list = ls())

library(rEDM)
library(deSolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(vegan)


setwd("C:\\Users\\Koushik\\Dropbox\\NTU Papers and ESM\\Multivariate data\\lake zurich")
lake_data <- read.csv("gpedm-regime-shifts-main/data/lake_plankton/lake_zurich.csv")
lake_data <- lake_data[13:516,]
lake_data$Date <- as.Date(paste(lake_data$year, lake_data$month, 1, sep = "-"))
lake_data$Time <- seq(1:nrow(lake_data))
lake_data <- lake_data[, c("Time", setdiff(names(lake_data), "Time"))]

data <- lake_data[,c(1,4:18)]



mean_data <- sapply(data[ , 2:14], mean, na.rm = TRUE)
top3_vars <- sort(mean_data, decreasing = TRUE)[1:3]
names(top3_vars)

top3_lake_data <- data[,paste(names(top3_vars))]


#####################################################
# plot time series of top 3

# Create each plot with a vertical line at 1996
gg1 <- ggplot(top3_lake_data, aes(x = lake_data[,19], y = top3_lake_data[,1])) +
  geom_line(color = "steelblue") +
  geom_vline(xintercept = as.Date("1996-01-01"), linetype = "dashed", color = "red", size = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(y = NULL) +
  ggtitle(paste(names(top3_vars[1])))

gg2 <- ggplot(top3_lake_data, aes(x = lake_data[,19], y = top3_lake_data[,2])) +
  geom_line(color = "steelblue") +
  geom_vline(xintercept = as.Date("1996-01-01"), linetype = "dashed", color = "red", size = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylab("Plankton biomass") +
  ggtitle(paste(names(top3_vars[2])))

gg3 <- ggplot(top3_lake_data, aes(x = lake_data[,19], y = top3_lake_data[,3])) +
  geom_line(color = "steelblue") +
  geom_vline(xintercept = as.Date("1996-01-01"), linetype = "dashed", color = "red", size = 1) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle(paste(names(top3_vars[3]))) +
  labs(y = NULL) +
  xlab("Time")

# Combine plots
### plot in a column
win.graph(40,90)
(gg1 / gg2 / gg3) +
  plot_layout(ncol = 1, heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = 'A')
#####################################################################################

########################################
# # combine all figures
# 
# library(gridExtra)
# left_col <- arrangeGrob(gg1, gg2, gg3, ncol = 1)
# mid_col <- arrangeGrob(gg1, gg2, gg3, ncol = 1)
# right_col <- arrangeGrob(gg1, gg2, gg3, ncol = 1)
# 
# grid.arrange(left_col, mid_col, right_col, ncol = 3)
# 



