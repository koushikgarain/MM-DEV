######## Check sensitive variables, Masuda et al Nat comm 2024 paper
# Clear workspace
rm(list = ls())

# Load libraries
library(tidyverse)
library(zoo)
library(reshape2)
library(R.matlab)
library(ggplot2)

#################################################################################
##### Load data
data_list <- list()

## Bioreactor data
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
time <- seq(1,nrow(biorec2))
biorec3 <- as.data.frame(cbind(time,biorec2))
data_list[[3]] <- biorec3

## Plankton data
lake_data <- read.csv("Data\\lake_zurich.csv")
lake_data <- lake_data[13:516,]
lake_data$Date <- as.Date(paste(lake_data$year, lake_data$month, 1, sep = "-"))
lake_data$Time <- seq(1:nrow(lake_data))
lake_data <- lake_data[, c("Time", setdiff(names(lake_data), "Time"))]
lake_data <- na.omit(lake_data)
data_new <- lake_data[,c(1,4:18)]
data_list[[4]] <- na.omit(data_new)

## Ricker model
kk=read.csv("Data\\500RemoveSeq1_RK21.csv",header=T)
data_list[[2]] <- na.omit(kk)

##lotka volterra 6-dim

# transcritical
mat_file_T <- "Data\\lot_vol_6_1_data.mat"
# # saddle-node
# mat_file_S <- "Data\\lot_vol_6_2_data.mat"

mat_data_T <- readMat(mat_file_T)

# Extract the relevant data (adjust accordingly based on data structure)
your_data_variable_T <- mat_data_T$data  # Replace with the actual variable name
your_data_variable_T <- as.data.frame(your_data_variable_T[15000:21800,1:6])

data_new2 <- na.omit(your_data_variable_T)
time <-  seq(1,nrow(data_new2))
data_list[[1]] <- cbind(time,data_new2)
####################################################################################

### Data  # window length 
L <- 250; L1=6500  # lotka volterra model
X <- data_list[[1]][,-1]
#
L <- 50; L1=400  # ricker model
X <- data_list[[2]][,-1]
#
L <- 25; L1=50  # Bioreactor data
X <- data_list[[3]][,-1]
#
L <- 50; L1=150  # lake zurich data
X <- data_list[[4]][,-1]


# window indices (adjust if you have domain knowledge)
idx1 <- 1:L                     # far from transition
idx2 <- L1:(L1 + L - 1)       # close to transition

X1 <- X[idx1, ]
X2 <- X[idx2, ]

# Covariance
C1 <- cov(X1)
C2 <- cov(X2)

## Objective function
compute_d <- function(C1, C2, S, L) {
  n <- length(S)
  
  C1s <- C1[S, S, drop = FALSE]
  C2s <- C2[S, S, drop = FALSE]
  
  # Mean of averaged variance (Eq. 5)
  mu1 <- sum(diag(C1s)) / n
  mu2 <- sum(diag(C2s)) / n
  
  # Variance of averaged variance (Eq. 6)
  var1 <- 2 * sum(C1s^2) / (n^2 * (L - 1))
  var2 <- 2 * sum(C2s^2) / (n^2 * (L - 1))
  
  d <- (mu2 - mu1) / sqrt(var1 + var2)
  return(d)
}

## Sensetive variable
d_single <- numeric(ncol(X))

for (i in seq_len(ncol(X))) {
  d_single[i] <- compute_d(C1, C2, S = i, L = L)
}

best_var <- which.max(d_single)

cat("Best sentinel variable index:", best_var, "\n")
cat("Maximum d value:", d_single[best_var], "\n")


ranked_vars <- order(d_single, decreasing = TRUE)

M <- 5                    # candidate pool size
candidates <- ranked_vars[1:M]

S <- c(candidates[1])       # start with best variable
d_current <- compute_d(C1, C2, S, L)

for (j in candidates[-1]) {
  d_new <- compute_d(C1, C2, c(S, j), L)
  
  if (d_new > d_current) {
    S <- c(S, j)
    d_current <- d_new
  }
}

#library(zoo)

window_size <- L
step_size   <- 1

ews <- rollapply(
  X,
  width = window_size,
  by = step_size,
  align = "right",
  FUN = function(window_data) {
    C <- cov(window_data[, S, drop = FALSE])
    mean(diag(C))
  },
  by.column = FALSE,
  fill = NA
)

# visualize
win.graph(30,90)
par(mfrow=c(4,1),oma=c(0,0,2,0))


plot(
  ews,
  type = "l",
  col = "blue", las=0, cex.axis=1.75, cex.lab=1.65,
  lwd = 2,
  xlab = "Time",
  ylab = "Variance"
)

## Lotka volterra model
abline(v = 6700, col="red", lwd=1.5, lty=2)
title("Lotka_volterra Model")

## Ricker model
abline(v = 500, col="red", lwd=1.5, lty=2)
title("Ricker Model")

## Bioreactor data
abline(v = 80, col="red", lwd=1.5, lty=2)
title("Bioreactor Data")

## Lake zurich data
abline(v = 250, col="red", lwd=1.5, lty=2)
title("Lake Zurich Data")


# 
# # top variables
# k <- 5
# top_k_vars <- order(d_single, decreasing = TRUE)[1:k]
# top_k_vars
# 
# # Check for MDEV
# new_variable <- X[,candidates]

