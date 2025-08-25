###### Figure 2: Model data analysis

rm(list = ls())
### Load SSR ###
library(rEDM)
library(dplyr)
library(R.matlab)
library(openxlsx)
library(vegan)

###
save.plot.file <- T
if(save.plot.file){pdf("figure_2.pdf", width = 9, height = 9)}else{win.graph(80,100)}
par(mfrow=c(3,3),oma=c(0,0,2,0))



#########################################################################################################
### Lotka-Volterra model: trans (Fig. 2A-C) and saddle bifurcation (Fig. 2D-F)

#lotka volterra 6-dim
# transcritical
mat_file_T <- "C:\\Users\\Koushik\\Dropbox\\NTU Papers and ESM\\Multivariate data\\Lot_Vol_6\\lot_vol_6_1_data.mat"
# saddle-node
mat_file_S <- "C:\\Users\\Koushik\\Dropbox\\NTU Papers and ESM\\Multivariate data\\Lot_Vol_6\\lot_vol_6_2_data.mat"

mat_data_T <- readMat(mat_file_T)
mat_data_S <- readMat(mat_file_S)

# Extract the relevant data (adjust accordingly based on data structure)
your_data_variable_T <- mat_data_T$data  # Replace with the actual variable name
your_data_variable_S <- mat_data_S$data

your_data_variable_S <- as.data.frame(your_data_variable_S[15000:25000,1:6])
your_data_variable_T <- as.data.frame(your_data_variable_T[15000:25000,1:6])

# After ISOMAP analysis
setwd("C:\\Users\\Koushik\\Dropbox\\NTU Papers and ESM\\Multivariate data\\Lot_Vol_6\\Figures")
mat_Isomap_T <- read.csv("lot_vol_6_1_Isomap.csv")
mat_Isomap_S <- read.csv("lot_vol_6_2_Isomap.csv")

####### trans bifurcation DRM-DEV analysis (Fig. 2A-C) #####
emp_time_series <- scale(mat_Isomap_T[,2:4])

# calculate E and w from S-map
E_best <-2

# form a matrix for new data set
new_data <- matrix(NA, nrow = nrow(emp_time_series), ncol = E_best)

# Determine the number of columns to copy from emp_time_series
num_cols <- min(E_best, ncol(emp_time_series))

# Copy data from emp_time_series to new_data
for (i in 1:num_cols) {
  new_data[, i] <- emp_time_series[, i]
}

# Calculate multivariate dev

E <- E_best
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 150
step_size <- 25

theta_best <-1

window_indices <- seq(window_size, NROW(new_data), step_size)
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


###### plot ts

par(mar=c(4.1,5.3,1,1)+0.1)

time <-seq(1:nrow(your_data_variable_T))
cols <- rainbow(ncol(your_data_variable_T))
# Plot the first variable
plot(time, your_data_variable_T[, 1], type = "l", col = cols[1], las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Values",
     ylim = c(0,1.3))

# Add the other lines
for (i in 2:ncol(your_data_variable_T)) {
  lines(time, your_data_variable_T[, i], col = cols[i])
}

# # Add legend
# legend("topright", legend = colnames(your_data_variable_S), col = cols, lty = 1, bty = "n")

polygon(c(6500,7500,7500,6500),c(-10,-10,20,20), col="red",density=15, border="white")
points(c(6500,6500),c(-100,200), type="l", col="red", lwd=1.5, lty=2)
mtext("A",side=3, at=0, line=0.5,cex=1.5)

#### plot dev

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0,xlim=c(0,10000), ylim=c(0.63,1.05), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.5, main="", xlab="Time", ylab="|DRM-DEV|")
points(c(0,125000),c(1,1), type="l", col="gray", lwd=3, lty=2)

mat = cbind(multi_result[,1], multi_result[,2], 1:length(multi_result))
n = 255
data_seq = seq(1, nrow(multi_result), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(multi_result[,1], multi_result[,2], type="p", cex=1, pch=20, col=cols)
lines(multi_result[,1],multi_result[,2], col='black', lwd=1)
#lines(multi_result[150:261,1],predict(lm(multi_result[150:261,2]~multi_result[150:261,1])), col='red', lwd=2)
polygon(c(6500,7500,7500,6500),c(-10,-10,20,20), col="red",density=15, border="white")
points(c(6500,6500),c(-100,200), type="l", col="red", lwd=1.5, lty=2)

mtext("B",side=3, at=0, line=0.5,cex=1.5)

## plot dev in complex

par(mar=c(4.1,2.5,1,4.5)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.65, cex.axis=1.75, type="l", xlab="Re(DRM-DEV)", ylab="", xlim=c(0.6,1), ylim=c(-0.2,0.2), yaxt="n")

n <- NROW(multi_result)

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(multi_result[,3], multi_result[,4], lty=1, lwd=1.5, type="p", cex=4.75, pch=20, col=cols)

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DRM-DEV)",side=4,line=3,cex=1.0)

mtext("C",side=3, at=0.6, line=0.5,cex=1.5)

####### saddle bifurcation DRM-DEV analysis (Fig. 2D-F) #############################################################
emp_time_series <- scale(mat_Isomap_S[,2:4])

# calculate E and w from S-map
E_best <-2

# form a matrix for new data set
new_data <- matrix(NA, nrow = nrow(emp_time_series), ncol = E_best)

# Determine the number of columns to copy from emp_time_series
num_cols <- min(E_best, ncol(emp_time_series))

# Copy data from emp_time_series to new_data
for (i in 1:num_cols) {
  new_data[, i] <- emp_time_series[, i]
}

# Calculate multivariate dev

E <- E_best
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 150
step_size <- 25

theta_best <-1

window_indices <- seq(window_size, NROW(new_data), step_size)
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


###### plot ts

par(mar=c(4.1,5.3,1,1)+0.1)

time <-seq(1:nrow(your_data_variable_S))
cols <- rainbow(ncol(your_data_variable_S))
# Plot the first variable
plot(time, your_data_variable_S[, 1], type = "l", col = cols[1], las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Values",
     ylim = c(0,1.3))

# Add the other lines
for (i in 2:ncol(your_data_variable_S)) {
  lines(time, your_data_variable_S[, i], col = cols[i])
}

# # Add legend
# legend("topright", legend = colnames(your_data_variable_S), col = cols, lty = 1, bty = "n")

polygon(c(8900,9900,9900,8900),c(-10,-10,20,20), col="red",density=15, border="white")
points(c(8900,8900),c(-100,200), type="l", col="red", lwd=1.5, lty=2)

mtext("D",side=3, at=0, line=0.5,cex=1.5)

#### plot dev

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0,xlim=c(0,10000), ylim=c(0.63,1.05), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.5, main="", xlab="Time", ylab="|DRM-DEV|")
points(c(0,125000),c(1,1), type="l", col="gray", lwd=3, lty=2)

mat = cbind(multi_result[,1], multi_result[,2], 1:length(multi_result))
n = 255
data_seq = seq(1, nrow(multi_result), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(multi_result[,1], multi_result[,2], type="p", cex=1, pch=20, col=cols)
lines(multi_result[,1],multi_result[,2], col='black', lwd=1)
#lines(multi_result[150:261,1],predict(lm(multi_result[150:261,2]~multi_result[150:261,1])), col='red', lwd=2)
polygon(c(8900,9900,9900,8900),c(-10,-10,20,20), col="red",density=15, border="white")
points(c(8900,8900),c(-100,200), type="l", col="red", lwd=1.5, lty=2)

mtext("E",side=3, at=0, line=0.5,cex=1.5)

## plot dev in complex

par(mar=c(4.1,2.5,1,4.5)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.65, cex.axis=1.75, type="l", xlab="Re(DRM-DEV)", ylab="", xlim=c(0.6,1), ylim=c(-0.2,0.2), yaxt="n")

n <- NROW(multi_result)

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(multi_result[,3], multi_result[,4], lty=1, lwd=1.5, type="p", cex=4.75, pch=20, col=cols)

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DRM-DEV)",side=4,line=3,cex=1.0)

mtext("F",side=3, at=0.6, line=0.5,cex=1.5)


########################################################################################
### Ricker model: DRM-DEV analysis (Fig. 2G-I)
OPath <- 'C:\\Users\\Koushik\\Dropbox\\NTU Papers and ESM\\DEV_Bifurcation-main\\figure_s16_HighDimensionalModels'
setwd(OPath)
kk=read.csv('500RemoveSeq1_RK21.csv',header=T)

# with ISOMAP
# Load a dataset 
X11 <- scale(kk[,-1])

# Calculate pairwise distance matrix
dist_matrix1 <- dist(X11)

# Perform nonlinear PCA using Isomap
isomap_result1 <- isomap(dist_matrix1, k = 10, ndim = 10)

# Extract coordinates of the reduced space
isomap_coordinates1 <- isomap_result1$points

#write.csv(isomap_coordinates1, "lot_vol_6_2_Isomap.csv")
##########

emp_time_series <- isomap_coordinates1[,(1:3)]


# Number of time lags for each variable
E_best<-2

# form a matrix for new data set
new_data <- matrix(NA, nrow = nrow(emp_time_series), ncol = E_best)

# Determine the number of columns to copy from emp_time_series
num_cols <- min(E_best, ncol(emp_time_series))

# Copy data from emp_time_series to new_data
for (i in 1:num_cols) {
  new_data[, i] <- emp_time_series[, i]
}

# Calculate multivariate dev 

E <- E_best
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 150
#step_size <- 25

# window_size <- 125
step_size <- 5

theta_best <-1

window_indices <- seq(window_size, NROW(new_data), step_size)
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

###### plot ts

par(mar=c(4.1,5.3,1,1)+0.1)

time <-seq(1:nrow(kk))
cols <- rainbow(21)
# Plot the first variable
plot(time, kk[, 2], type = "l", col = cols[1], las=0, cex.axis=1.75, cex.lab=1.65, xlab = "Time", ylab = "Values",
     ylim = c(0,1.5))

# Add the other lines
for (i in 3:6) {
  lines(time, kk[, i], col = cols[i])
}

# # Add legend
# legend("topright", legend = colnames(your_data_variable_S), col = cols, lty = 1, bty = "n")

#polygon(c(8900,9900,9900,8900),c(-10,-10,20,20), col="red",density=15, border="white")
points(c(500,500),c(-100,200), type="l", col="red", lwd=1.5, lty=2)

mtext("G",side=3, at=0, line=0.5,cex=1.5)

#### plot dev

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0,xlim=c(0,500), ylim=c(0.3,1.02), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.5, main="", xlab="Time", ylab="|DRM-DEV|")
points(c(0,1250),c(1,1), type="l", col="gray", lwd=3, lty=2)

mat = cbind(multi_result[,1], multi_result[,2], 1:length(multi_result))
n = 255
data_seq = seq(1, nrow(multi_result), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(multi_result[,1], multi_result[,2], type="p", cex=1, pch=20, col=cols)
lines(multi_result[,1],multi_result[,2], col='black', lwd=1)
#lines(multi_result[150:261,1],predict(lm(multi_result[150:261,2]~multi_result[150:261,1])), col='red', lwd=2)
#polygon(c(8900,9900,9900,8900),c(-10,-10,20,20), col="red",density=15, border="white")
points(c(500,500),c(-100,200), type="l", col="red", lwd=1.5, lty=2)

mtext("H",side=3, at=0, line=0.5,cex=1.5)

## plot dev in complex

par(mar=c(4.1,2.5,1,4.5)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.65, cex.axis=1.75, type="l", xlab="Re(DRM-DEV)", ylab="", xlim=c(-1,0.5), ylim=c(-0.2,0.2), yaxt="n")

n <- NROW(multi_result)

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(multi_result[,3], multi_result[,4], lty=1, lwd=1.5, type="p", cex=4.75, pch=20, col=cols)

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DRM-DEV)",side=4,line=3,cex=1.0)

mtext("I",side=3, at=-1, line=0.5,cex=1.5)


if(save.plot.file){dev.off()}
