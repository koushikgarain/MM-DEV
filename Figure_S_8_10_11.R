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
lake_data2 <- lake_data[1:260,]
data_new2 <- lake_data2[,c(1,4:18)]
data_list[[4]] <- na.omit(data_new)
data_list[[5]] <- na.omit(data_new2)  # data upto 01-12-1999, 260 data points

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
your_data_variable_T <- as.data.frame(your_data_variable_T[15001:21800,1:6])

data_new2 <- na.omit(your_data_variable_T)
time <-  seq(1,nrow(data_new2))
data_list[[1]] <- cbind(time,data_new2)
#emp_time_series <- scale(mat_Isomap_T[1:6800,2:4])

#################################################################################
cov_trace_list <- list()  
cov_det_list <- list()
cov_lamda_list <- list()
cor_abs_list <- list()  
cor_lamda_list <- list()
mar_list <- list()  

for (k in seq(1,5)) {
  data <- data_list[[k]]
  time <- data[,1]
  X <- as.matrix(data[,-1])
  
  T <- nrow(X)
  p <- ncol(X)
  
  X <- scale(X)
  
  ## Rolling window
  #window_size <- 50  # old, for all datasets
  #window_size <- T/2  # new, 50% of time points
  step <- 1
  
  # three different window size: 25%, 50%, 75%
  for (ws in c(round(T/4),T/2,round(3*T/4))) {
    window_size <- ws
  
  windows <- seq(1, T-window_size, by=step)
  
  ## Covariance matrix
  cov_trace <- c()
  cov_det <- c()
  cov_lambda_max <- c()
  
  for(w in windows){
    
    Xw <- X[w:(w+window_size-1), ]
    
    Cmat <- cov(Xw)
    
    eig <- eigen(Cmat)$values
    
    cov_trace <- c(cov_trace, sum(eig))
    cov_det <- c(cov_det, det(Cmat))
    cov_lambda_max <- c(cov_lambda_max, max(eig))
  }
  k_w <- paste("k", k, "_w", ws, sep = "")
  cov_trace_list[[k_w]] <- cov_trace
  cov_det_list[[k_w]] <- cov_det
  cov_lamda_list[[k_w]] <- cov_lambda_max
  
  
  ## Correlation matrix
  mean_abs_corr <- c()
  corr_lambda_max <- c()
  
  for(w in windows){
    
    Xw <- X[w:(w+window_size-1), ]
    
    Rmat <- cor(Xw)
    
    eigR <- eigen(Rmat)$values
    
    mean_abs_corr <- c(mean_abs_corr,
                       mean(abs(Rmat[upper.tri(Rmat)])))
    
    corr_lambda_max <- c(corr_lambda_max, max(eigR))
  }
  cor_abs_list[[k_w]] <- mean_abs_corr
  cor_lamda_list[[k_w]] <- corr_lambda_max
  
  
  ## Avoid singularity in matrix (Bioreactor data)
  mar_spectral_radius <- c()
  
  lambda <- 1e-4
  #window_size <- 50  # old, for all datasets
  #window_size <- T/2  # new, 50% of time points
  
  for(i in windows){
    
    window_data <- X[i:(i + window_size), ]
    
    # X_t and X_{t+1}
    XX <- t(window_data[1:(nrow(window_data)-1), ])
    YY <- t(window_data[2:nrow(window_data), ])
    
    # Ridge MAR estimator
    A <- YY %*% t(XX) %*%
      solve(XX %*% t(XX) + lambda * diag(nrow(XX)))
    
    eigA <- eigen(A, only.values = TRUE)$values
    
    mar_spectral_radius <- c(mar_spectral_radius,
                             max(Mod(eigA)))
  }
  mar_list[[k_w]] <- mar_spectral_radius
  
}
}

######################################################################################
## Plot

### Covariance trace plot with different window sizes

win.graph(50,50)
par(mfrow=c(2,2),oma=c(0,0,2,0))

plot(seq(1700,6800-1),cov_trace_list[[1]],type = "l", col = "blue",ylim = c(0.5,3.6), las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Lotka_volterra Model")
lines(seq(3400,6800-1),cov_trace_list[[2]], col="red", lwd=1.5)
lines(seq(5100,6800-1),cov_trace_list[[3]], col="green", lwd=1.5)
abline(v = 6700, col="red", lwd=1.5, lty=2)
title("Trace of the Covariance Matrix")


plot(seq(125,500-1),cov_trace_list[[4]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Ricker Model")
lines(seq(250,500-1),cov_trace_list[[5]], col="red", lwd=1.5)
lines(seq(375,500-1),cov_trace_list[[6]], col="green", lwd=1.5)
abline(v = 500, col="red", lwd=1.5, lty=2)


plot(seq(28,110-1),cov_trace_list[[7]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Bioreactor Data")
lines(seq(55,110-1),cov_trace_list[[8]], col="red", lwd=1.5)
lines(seq(82,110-1),cov_trace_list[[9]], col="green", lwd=1.5)
abline(v = 80, col="red", lwd=1.5, lty=2)


plot(lake_data2[65:259,19],cov_trace_list[[13]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Lake Zurich Data")
lines(lake_data2[130:259,19],cov_trace_list[[14]], col="red", lwd=1.5)
lines(lake_data2[195:259,19],cov_trace_list[[15]], col="green", lwd=1.5)
abline(v = as.Date("1996-01-01"), col="red", lwd=1.5, lty=2)

############################################################################################
### Covariance plot

win.graph(50,90)
par(mfrow=c(4,2),oma=c(0,0,2,0))

plot(seq(1700,6800-1),cov_trace_list[[1]],type = "l", col = "blue",ylim = c(0.5,3.6), las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Lotka_volterra Model")
lines(seq(3400,6800-1),cov_trace_list[[2]], col="red", lwd=1.5)
lines(seq(5100,6800-1),cov_trace_list[[3]], col="green", lwd=1.5)
abline(v = 6700, col="red", lwd=1.5, lty=2)
title("Trace of the Covariance \n Matrix")

plot(seq(1700,6800-1),cov_lamda_list[[1]],type = "l", col = "blue",ylim = c(0.5,3.6), las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Lotka_volterra Model")
lines(seq(3400,6800-1),cov_lamda_list[[2]], col="red", lwd=1.5)
lines(seq(5100,6800-1),cov_lamda_list[[3]], col="green", lwd=1.5)
abline(v = 6700, col="red", lwd=1.5, lty=2)
title("Maximum Eigenvalue of the \n Covariance Matrix")

plot(seq(125,500-1),cov_trace_list[[4]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Ricker Model")
lines(seq(250,500-1),cov_trace_list[[5]], col="red", lwd=1.5)
lines(seq(375,500-1),cov_trace_list[[6]], col="green", lwd=1.5)
abline(v = 500, col="red", lwd=1.5, lty=2)

plot(seq(125,500-1),cov_lamda_list[[4]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "")
lines(seq(250,500-1),cov_lamda_list[[5]], col="red", lwd=1.5)
lines(seq(375,500-1),cov_lamda_list[[6]], col="green", lwd=1.5)
abline(v = 500, col="red", lwd=1.5, lty=2)

plot(seq(28,110-1),cov_trace_list[[7]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Bioreactor Data")
lines(seq(55,110-1),cov_trace_list[[8]], col="red", lwd=1.5)
lines(seq(82,110-1),cov_trace_list[[9]], col="green", lwd=1.5)
abline(v = 80, col="red", lwd=1.5, lty=2)

plot(seq(28,110-1),cov_lamda_list[[7]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "")
lines(seq(55,110-1),cov_lamda_list[[8]], col="red", lwd=1.5)
lines(seq(82,110-1),cov_lamda_list[[9]], col="green", lwd=1.5)
abline(v = 80, col="red", lwd=1.5, lty=2)

plot(lake_data2[65:259,19],cov_trace_list[[13]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Lake Zurich Data")
lines(lake_data2[130:259,19],cov_trace_list[[14]], col="red", lwd=1.5)
lines(lake_data2[195:259,19],cov_trace_list[[15]], col="green", lwd=1.5)
abline(v = as.Date("1996-01-01"), col="red", lwd=1.5, lty=2)

plot(lake_data2[65:259,19],cov_lamda_list[[13]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "")
lines(lake_data2[130:259,19],cov_lamda_list[[14]], col="red", lwd=1.5)
lines(lake_data2[195:259,19],cov_lamda_list[[15]], col="green", lwd=1.5)
abline(v = as.Date("1996-01-01"), col="red", lwd=1.5, lty=2)



#########################################################################################  
#### Correlation plot

win.graph(50,90)
par(mfrow=c(4,2),oma=c(0,0,2,0))


plot(seq(1700,6800-1),cor_abs_list[[1]],type = "l", col = "blue",ylim = c(0.5,1.2), las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Lotka_volterra Model")
lines(seq(3400,6800-1),cor_abs_list[[2]], col="red", lwd=1.5)
lines(seq(5100,6800-1),cor_abs_list[[3]], col="green", lwd=1.5)
abline(v = 6700, col="red", lwd=1.5, lty=2)
title("Absolute value of the \n Correlation Matrix")

plot(seq(1700,6800-1),cor_lamda_list[[1]],type = "l", col = "blue",ylim = c(4,6.2), las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Lotka_volterra Model")
lines(seq(3400,6800-1),cor_lamda_list[[2]], col="red", lwd=1.5)
lines(seq(5100,6800-1),cor_lamda_list[[3]], col="green", lwd=1.5)
abline(v = 6700, col="red", lwd=1.5, lty=2)
title("Maximum Eigenvalue of the \n Correlation Matrix")

plot(seq(125,500-1),cor_abs_list[[4]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Ricker Model")
lines(seq(250,500-1),cor_abs_list[[5]], col="red", lwd=1.5)
lines(seq(375,500-1),cor_abs_list[[6]], col="green", lwd=1.5)
abline(v = 500, col="red", lwd=1.5, lty=2)

plot(seq(125,500-1),cor_lamda_list[[4]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "")
lines(seq(250,500-1),cor_lamda_list[[5]], col="red", lwd=1.5)
lines(seq(375,500-1),cor_lamda_list[[6]], col="green", lwd=1.5)
abline(v = 500, col="red", lwd=1.5, lty=2)

plot(seq(28,110-1),cor_abs_list[[7]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Bioreactor Data")
lines(seq(55,110-1),cor_abs_list[[8]], col="red", lwd=1.5)
lines(seq(82,110-1),cor_abs_list[[9]], col="green", lwd=1.5)
abline(v = 80, col="red", lwd=1.5, lty=2)

plot(seq(28,110-1),cor_lamda_list[[7]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "")
lines(seq(55,110-1),cor_lamda_list[[8]], col="red", lwd=1.5)
lines(seq(82,110-1),cor_lamda_list[[9]], col="green", lwd=1.5)
abline(v = 80, col="red", lwd=1.5, lty=2)

plot(lake_data2[65:259,19],cor_abs_list[[13]],type = "l", col = "blue",ylim = c(0.15,0.3), las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Lake Zurich Data")
lines(lake_data2[130:259,19],cor_abs_list[[14]], col="red", lwd=1.5)
lines(lake_data2[195:259,19],cor_abs_list[[15]], col="green", lwd=1.5)
abline(v = as.Date("1996-01-01"), col="red", lwd=1.5, lty=2)

plot(lake_data2[65:259,19],cor_lamda_list[[13]],type = "l", col = "blue",ylim = c(3,5), las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "")
lines(lake_data2[130:259,19],cor_lamda_list[[14]], col="red", lwd=1.5)
lines(lake_data2[195:259,19],cor_lamda_list[[15]], col="green", lwd=1.5)
abline(v = as.Date("1996-01-01"), col="red", lwd=1.5, lty=2)


### MAR Coefficient plot

win.graph(50,50)
par(mfrow=c(2,2),oma=c(0,0,2,0))

plot(seq(1700,6800-1),mar_list[[1]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Lotka_volterra Model")
lines(seq(3400,6800-1),mar_list[[2]], col="red", lwd=1.5)
lines(seq(5100,6800-1),mar_list[[3]], col="green", lwd=1.5)
abline(v = 6700, col="red", lwd=1.5, lty=2)
title("MAR Spectral Radius")

plot(seq(125,500-1),mar_list[[4]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Ricker Model")
lines(seq(250,500-1),mar_list[[5]], col="red", lwd=1.5)
lines(seq(375,500-1),mar_list[[6]], col="green", lwd=1.5)
abline(v = 500, col="red", lwd=1.5, lty=2)

plot(seq(28,110-1),mar_list[[7]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Bioreactor Data")
lines(seq(55,110-1),mar_list[[8]], col="red", lwd=1.5)
lines(seq(82,110-1),mar_list[[9]], col="green", lwd=1.5)
abline(v = 80, col="red", lwd=1.5, lty=2)

plot(lake_data2[65:259,19],mar_list[[13]],type = "l", col = "blue", las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Lake Zurich Data")
lines(lake_data2[130:259,19],mar_list[[14]], col="red", lwd=1.5)
lines(lake_data2[195:259,19],mar_list[[15]], col="green", lwd=1.5)
abline(v = as.Date("1996-01-01"), col="red", lwd=1.5, lty=2)

