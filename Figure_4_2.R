######## Figure 4 MM-DEV plot
rm(list = ls())

# DEV analysis
library(rEDM)
library(dplyr)
library(R.matlab)
library(openxlsx)
library(vegan)

setwd("C:\\Users\\Koushik\\Dropbox\\NTU Papers and ESM\\Multivariate data\\lake zurich")
data_new <- read.csv("gpedm-regime-shifts-main/data/lake_plankton/lake_zurich.csv")
lake_data2 <- data_new[13:516,c(1:15)]

data1 <- na.omit(lake_data2)
data2 <- data1[,-1:-2]
# ISOMAP of both chemicals and species

X11 <- scale(data2)
dist_matrix1 <- dist(X11)
isomap_result1 <- isomap(dist_matrix1, k = 5, ndim = 10)
isomap_coordinates1 <- isomap_result1$points

emp_time_series <- isomap_coordinates1[,(1:3)]

# Number of time lags for each variable
E_best<-12
print(E_best)
d1 <- E_best - 3

nolags <-c(9,0,0)   # total must be equal to E_best-total variables, if E_best>total variables

# form a matrix for new data set
new_data <- matrix(NA, nrow = nrow(emp_time_series), ncol = E_best)

# Determine the number of columns to copy from emp_time_series
num_cols <- min(E_best, ncol(emp_time_series))

# Copy data from emp_time_series to new_data
for (i in 1:num_cols) {
  new_data[, i] <- emp_time_series[, i]
}

# If E_best is greater than the number of columns in emp_time_series

if (E_best > ncol(emp_time_series)) {
  
  
  for (i in 1:nolags[1]) {
    new_data[(i+1):nrow(new_data), num_cols + i] <- emp_time_series[1:(nrow(new_data)-i), 1]
  } 
  
  if (nolags[2] > 0) {  
    for (j in 1:nolags[2]) {
      new_data[(j+1):nrow(new_data), num_cols + i+j] <- emp_time_series[1:(nrow(new_data)-j), 2]
    }} 
  
  if (nolags[3] > 0) {  
    for (l in 1:nolags[3]) {
      new_data[(l+1):nrow(new_data), num_cols + i+j+l] <- emp_time_series[1:(nrow(new_data)-l), 3]
    }}   
  
}

# DEV analysis

E <- E_best
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 70
step_size <- 1

theta_best <-1

window_indices <- seq(window_size, NROW(emp_time_series), step_size)
matrix_result <- matrix(NaN, nrow = length(window_indices), ncol = 4)
index <- 0

#create matrix to save data
for(j in window_indices)
{
  index <- index + 1
  rolling_windows <- lapply(1:E, function(i) {
    new_data[(j - window_size + 1):j, i]
  })
  
  n <- window_size
  vectors <- matrix(NaN, nrow = n, ncol = E+1)
  
  vectors[, 1] <- c(1:n)
  for (i in 1:E) {
    rolling_window <- rolling_windows[[i]]
    vectors[, i + 1] <- (rolling_window - mean(rolling_window)) / sd(rolling_window)
  }
  colnames(vectors) <- c("time", paste0("N", 1:E))
  
  # # calculate best theta among each variables
  # theta_best <- NULL  # Initialize variable to store the maximum theta value
  # for (i in 1:E) {
  # smap <- s_map(vectors[, i + 1], E=E, tau=tau, theta=theta, silent=TRUE)
  # best <- order(-smap$rho)[1]
  # theta_best_all <- smap[best,]$theta
  # # Update max_theta if the current theta_best is greater
  # if (is.null(theta_best) || theta_best_all > theta_best) {
  #   theta_best <- theta_best_all
  # }
  # }
  
  if(T) {  # rEDM package ver 1.2.3
    smap_co <- list()
    for(k in 1:E) {
      out_block <- block_lnlp(vectors, norm=2, method="s-map", theta=theta_best, target=k, 
                              stats_only=FALSE, exclusion_radius=0, first_column_time=TRUE, 
                              columns=colnames(vectors), save_smap_coefficients=TRUE, tp=1)
      smap_co[[k]] <- out_block$smap_coefficients[[1]]
    }
  }
  
  # Matrix for saving eigenvalues
  matrix_eigen <- matrix(NA, nrow = NROW(smap_co[[1]]), ncol = 3)
  
  
  for(k in 1:NROW(smap_co[[1]]))
  {
    if(!is.na(smap_co[[1]][k,1]))
    {
      M <- do.call(rbind, lapply(smap_co, function(mat) mat[k, 1:E]))
      M_eigen <- eigen(M)$values
      lambda1 <- M_eigen[order(abs(M_eigen))[E]]
      
      
      matrix_eigen[k,1] <- abs(lambda1)
      matrix_eigen[k,2] <- Re(lambda1)
      matrix_eigen[k,3] <- Im(lambda1)
      
    }
  }
  
  # calculate save data
  matrix_result[index,1] <- j
  matrix_result[index,2] <- mean(matrix_eigen[,1],na.rm=TRUE)
  matrix_result[index,3] <- mean(matrix_eigen[,2],na.rm=TRUE)
  matrix_result[index,4] <- mean(matrix_eigen[,3],na.rm=TRUE)
}

multi_result <- matrix_result

# plot the result 
plot_data_start <- window_size
plot_data_indices <- seq(plot_data_start, nrow(data1), by = step_size)
# Extract corresponding year and month
plot_data <- cbind(data1[plot_data_indices, 1:2],multi_result)

# Convert Year and Month into Date format (assuming day = 1)
plot_data$Date <- as.Date(paste(plot_data$year, plot_data$month, 1, sep = "-"))


#win.graph(40,70)
win.graph(40,60)
par(mfrow=c(2,1),oma=c(0,0,2,0))
#############################################
#### plot dev
plot_data <- plot_data[1:200,]
#par(mar=c(5,5,2,2)+0.1)


mat = cbind(multi_result[,1], multi_result[,2], 1:length(multi_result))
n = 255
data_seq = seq(1, nrow(multi_result), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

plot(plot_data$Date, plot_data[, 4], ylim=c(0.85,1.02),
     type = "p", pch = 20, col = cols, cex = 1.5, lwd = 1.5,
     cex.axis=1.75, cex.lab=1.5, main="", xlab="Time", ylab="|MM-DEV|")

lines(plot_data$Date,plot_data[, 4], col='black', lwd=1)
points(c(0,125000),c(1,1), type="l", col="gray", lwd=3, lty=2)
abline(v = as.Date("1996-06-01"), col = "red", lty = 2, lwd = 3)

mtext("D",side=3, at=as.Date("1983-06-01"), line=0.5,cex=1.5)

## plot dev in complex

#par(mar=c(5,5,2,2)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.65, cex.axis=1.75, type="l", xlab="Re(MM-DEV)", ylab="Im(MM-DEV)", xlim=c(0.6,1), ylim=c(-0.8,0.1))

n <- NROW(multi_result)

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(multi_result[,3], multi_result[,4], lty=1, lwd=1, type="p", cex=2.75, pch=20, col=cols)

#2nd y-axis
#axis(4, las=0, col="black", cex.axis=1.65)
#mtext("Im(DRM-DEV)",side=4,line=3,cex=1.0)

mtext("E",side=3, at=0.6, line=0.5,cex=1.5)

