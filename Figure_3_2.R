############ Figure 3 MM-DEV plot
rm(list = ls())

library(rEDM)
library(deSolve)
library(dplyr)
library(readxl)
library(vegan)
library(ggplot2)

# read the data set 
D4_reactor_Expertise <- read_excel("Data\\UASB_daughter_reactors_environmental_parameter_table_first110days_20241230update_Expertise.xlsx")
D4_reactor_Expertise <- as.data.frame(D4_reactor_Expertise)
D4_reactor_Expertise <- na.omit(D4_reactor_Expertise)

### gas production multiply 3.1
D4_reactor_Expertise[,11] <- D4_reactor_Expertise[,11]*3.1

# KTU data
KTU_biorec_quan=read.csv("Data\\KTU_biorec_quan_D4.csv")
KTU_biorec_quan_D4 <- (KTU_biorec_quan[,c(-1,-2)])

#### Chemical table analysis

# ISOMAP analysis 

X1 <- D4_reactor_Expertise[,4:11]
#X1 <- D4_reactor_Navie[,4:27]
X11 <- scale(X1)

# Calculate pairwise distance matrix

dist_matrix1 <- dist(X11)

# Perform nonlinear PCA using Isomap

isomap_result1 <- isomap(dist_matrix1, k = 5, ndim = 3)

# Extract coordinates of the reduced space

isomap_coordinates1 <- isomap_result1$points
############################

emp_time_series <- isomap_coordinates1[1:90,(1:3)]

# Number of time lags for each variable
E_best<-1
print(E_best)
d1 <- E_best - 3

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
window_size <- 25
step_size <- 1

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

# plot the result 
### plot in a column
win.graph(30,50)
par(mfrow=c(2,1),oma=c(0,0,2,0))

# plot dev

#par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0,xlim=c(0,90), ylim=c(0.2,1.02), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.5, main="", xlab="Time (Days)", ylab="|MM-DEV|")
points(c(0,150),c(1,1), type="l", col="gray", lwd=3, lty=2)

mat = cbind(multi_result[,1], multi_result[,2], 1:length(multi_result))
n = 255
data_seq = seq(1, nrow(multi_result), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(multi_result[,1], multi_result[,2], type="p", cex=2, pch=20, col=cols)
lines(multi_result[,1],multi_result[,2], col='black', lwd=1)
abline(v = 80, col = "red", lty = 2, lwd = 3)

mtext("E",side=3, at=0, line=0.5,cex=1.5)

## plot dev in complex

#par(mar=c(4.1,5.3,2,1)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.65, cex.axis=1.75, type="l", xlab="Re(MM-DEV)", ylab="Im(MM-DEV)", xlim=c(0.6,1), ylim=c(-0.2,0.2))

n <- NROW(multi_result)

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(multi_result[,3], multi_result[,4], lty=1, lwd=1.5, type="p", cex=2.75, pch=20, col=cols)
mtext("F",side=3, at=0.6, line=0.5,cex=1.5)
###########################################################################################
###### KTU table analysis
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

biorec2 <- na.omit(filtered_df)
# ISOMAP analysis 

X1 <- biorec2
X11 <- scale(X1)

# Calculate pairwise distance matrix

dist_matrix1 <- dist(X11)

# Perform nonlinear PCA using Isomap

isomap_result1 <- isomap(dist_matrix1, k = 5, ndim = 3)

# Extract coordinates of the reduced space

isomap_coordinates1 <- isomap_result1$points
############################

emp_time_series <- isomap_coordinates1[1:90,(1:3)]

# Number of time lags for each variable
E_best<-2
print(E_best)
d1 <- E_best - 3

# form a matrix for new data set
new_data <- matrix(NA, nrow = nrow(emp_time_series), ncol = E_best)

# Determine the number of columns to copy from emp_time_series
num_cols <- min(E_best, ncol(emp_time_series))

# Copy data from emp_time_series to new_data
for (i in 1:num_cols) {
  new_data[, i] <- emp_time_series[, i]
}

# Calculate  MM-DEV 

E <- E_best
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 25
step_size <- 1

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

# plot the result
### plot in a column
win.graph(30,50)
par(mfrow=c(2,1),oma=c(0,0,2,0))

# plot dev

#par(mar=c(4.1,5.3,2,1)+0.1)

plot(0,0,xlim=c(0,90), ylim=c(0.7,1.1), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.5, main="", xlab="Time (Days)", ylab="|MM-DEV|")
points(c(0,150),c(1,1), type="l", col="gray", lwd=3, lty=2)

mat = cbind(multi_result[,1], multi_result[,2], 1:length(multi_result))
n = 255
data_seq = seq(1, nrow(multi_result), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(multi_result[,1], multi_result[,2], type="p", cex=2, pch=20, col=cols)
lines(multi_result[,1],multi_result[,2], col='black', lwd=1)
abline(v = 80, col = "red", lty = 2, lwd = 3)

mtext("G",side=3, at=0, line=0.5,cex=1.5)

## plot dev in complex

#par(mar=c(4.1,5.3,2,1)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.65, cex.axis=1.75, type="l", xlab="Re(MM-DEV)", ylab="Im(MM-DEV)", xlim=c(0.6,1), ylim=c(-0.2,0.2))

n <- NROW(multi_result)

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(multi_result[,3], multi_result[,4], lty=1, lwd=1.5, type="p", cex=2.75, pch=20, col=cols)

mtext("H",side=3, at=0.6, line=0.5,cex=1.5)

