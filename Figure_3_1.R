####### Figure 3 time series plot
rm(list = ls())

library(rEDM)
library(deSolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(vegan)
library(readxl)

## Data
KTU_biorec_quan=read.csv("Data\\KTU_biorec_quan_D4.csv")
KTU_biorec_quan_D4 <- (KTU_biorec_quan[,c(-1,-2)])
biorec1 <- KTU_biorec_quan_D4

df <- biorec1

# Calculate the number of non-zero entries in each row
non_zero_counts <- rowSums(df != 0)

# Calculate the total number of columns in the data frame
total_columns <- ncol(df)

# Calculate the percentage of non-zero entries for each row
non_zero_percentage <- non_zero_counts / total_columns

# Get the row numbers where non-zero percentage is greater than 90%
row_numbers <- which(non_zero_percentage > 0.90)

# Filter rows with more than 90% non-zero terms
filtered_df <- df[non_zero_percentage > 0.90, ]
filtered_df <- t(filtered_df)

biorec2 <- data.frame(na.omit(filtered_df))

# mean
mean_data <- sapply(biorec2[,1:153], mean, na.rm = TRUE)
top3_vars <- sort(mean_data, decreasing = TRUE)[1:3]
names(top3_vars)

top3_lake_data <- biorec2[,paste(names(top3_vars))]
top3_lake_data$Time <- seq(1:nrow(top3_lake_data))

#####################################################
# plot time series of top 2

# Create each plot with a vertical line at 1996
gg1 <- ggplot(top3_lake_data, aes(x = Time, y = top3_lake_data[,1])) +
  geom_line(color = "steelblue") +
  geom_vline(xintercept = 80, linetype = "dashed", color = "red", size = 1) +
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


gg2 <- ggplot(top3_lake_data, aes(x = Time, y = top3_lake_data[,2])) +
  geom_line(color = "darkgreen") +
  geom_vline(xintercept = 80, linetype = "dashed", color = "red", size = 1) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle(paste(names(top3_vars[2]))) +
  labs(y = NULL) +
  xlab("Time (Days)")

# Combine plots
### plot in a column
win.graph(30,50)
(gg1 / gg2) +
  plot_layout(ncol = 1, heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A')
#####################################################################################

########### for chemical data
D4_reactor_Navie <- read_excel("Data\\UASB_daughter_reactors_environmental_parameter_table_first110days_20241230update_Navie.xlsx")

#D4_reactor_Expertise <- as.data.frame(D4_reactor_Expertise)
D4_reactor_Navie <- as.data.frame(D4_reactor_Navie)

#D4_reactor_Expertise <- na.omit(D4_reactor_Expertise)
D4_reactor_Navie <- na.omit(D4_reactor_Navie)

### gas production multiply 3.1
D4_reactor_Navie[,c(18,26,27)] <- D4_reactor_Navie[,c(18,26,27)]*3.1

mean_data <- sapply(D4_reactor_Navie[,4:27], mean, na.rm = TRUE)
top3_vars <- sort(mean_data, decreasing = TRUE)[1:5]
names(top3_vars)

top3_lake_data <- D4_reactor_Navie[,paste(names(top3_vars))]
top3_lake_data$Time <- seq(1:nrow(top3_lake_data))

#####################################################
# plot time series of top 2 (must be in Expertise dataset, 2 and 5)

# Create each plot with a vertical line at 1996
gg1 <- ggplot(top3_lake_data, aes(x = Time, y = top3_lake_data[,2])) +
  geom_line(color = "steelblue") +
  geom_vline(xintercept = 80, linetype = "dashed", color = "red", size = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(y = NULL) +
  ggtitle(paste(names(top3_vars[2])))

# gg2 <- ggplot(top3_lake_data, aes(x = Time, y = top3_lake_data[,2])) +
#   geom_line(color = "firebrick") +
#   geom_vline(xintercept = 80, linetype = "dashed", color = "red", size = 1) +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x  = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   ) +
#   ylab("Microorganisms") +
#   ggtitle(paste(names(top3_vars[2])))

gg2 <- ggplot(top3_lake_data, aes(x = Time, y = top3_lake_data[,5])) +
  geom_line(color = "darkgreen") +
  geom_vline(xintercept = 80, linetype = "dashed", color = "red", size = 1) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle(paste(names(top3_vars[5]))) +
  labs(y = NULL) +
  xlab("Time (Days)")

# Combine plots
### plot in a column
win.graph(30,50)
(gg1 / gg2) +
  plot_layout(ncol = 1, heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

