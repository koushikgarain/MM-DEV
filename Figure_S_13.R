rm(list = ls())
### Load SSR ###
library(rEDM)
library(dplyr)
library(R.matlab)
library(openxlsx)
library(vegan)

# visualize
win.graph(50,50)
par(mfrow=c(2,2),oma=c(0,0,2,0))

#########################################################################################################
### load Data directly from 'Different_indicators_3_sv'

## mDEV for Lotka-volterra, Bioreactor data, Lake zurich data
new_data <- X[,S]

# Number of time lags for each variable
E_best<-ncol(new_data)

# Calculate multivariate dev 

E <- E_best
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- L
#step_size <- 25

# window_size <- 125
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

#### plot dev

### Lotka volterra model
par(mar=c(4.5,4.5,2,1))
plot(multi_result[,1], multi_result[,2],type = "l",
     col = "blue", las=0, cex.axis=1.25,lwd = 2, cex.lab=1.4, xlab="Time", ylab="|mDEV|")
points(c(0,125000),c(1,1), type="l", col="gray", lwd=3, lty=2)
abline(v = 6700, col="red", lwd=2, lty=2)
title("Lotka-volterra Model")
mtext("A",side=3, at=L, line=0.5,cex=1.5)

### Bioreactor data
par(mar=c(4.5,4.5,2,1))
plot(multi_result[,1], multi_result[,2],type = "l",
     col = "blue", las=0, cex.axis=1.25,lwd = 2, cex.lab=1.4, xlab="Time", ylab="|mDEV|")
points(c(0,125000),c(1,1), type="l", col="gray", lwd=3, lty=2)
abline(v = 80, col="red", lwd=2, lty=2)
title("Bioreactor Data")
mtext("B",side=3, at=L, line=0.5,cex=1.5)

### Lake zurich data
par(mar=c(4.5,4.5,2,1))
plot(lake_data[50:500,19], multi_result[,2],type = "l",
     col = "blue", las=0, cex.axis=1.25,lwd = 2, cex.lab=1.4, xlab="Time", ylab="|mDEV|")
points(c(0,125000),c(1,1), type="l", col="gray", lwd=3, lty=2)
abline(v = as.Date("1996-01-01"), col="red", lwd=2, lty=2)
title("Lake Zurich Data")
mtext("C",side=3, at=as.Date("1982-01-01"), line=0.5,cex=1.5)


#mtext("B",side=3, at=0, line=0.5,cex=1.5)
##############################

## univariate DEV for Ricker model, 
top_ricker_data <- new_data <- X[,S]
top_ricker_data <- scale(top_ricker_data)

result.dev=NULL
raw_time_series <- top_ricker_data
#raw_time_series <- X11[,1]

###set up the parameter ranges for performing cross-validation
E <- seq(1,12,by=1)
tau <- seq(1,1,by=1)
theta <- 0

window_size <- seq(50, 50, by=25)
step_size <- 1

### start algorithm
matrix_results <- data.frame()

for(i in 1:length(window_size))
{
  for(j in 1:length(E))
  {
    for(k in 1:length(tau))
    {
      for(l in 1:length(theta))
      {
        matrix_rho <- c()
        m <- 0
        
        while(m <= length(raw_time_series) - window_size[i] - step_size)
        {
          raw_ts_part <- raw_time_series[(m + 1):(m + window_size[i])]
          time_series <- (raw_ts_part - mean(raw_ts_part, na.rm=TRUE))/sd(raw_ts_part, na.rm=TRUE)
          
          smap <- s_map(time_series, E=E[j], tau=tau[k], theta=theta[l], silent=TRUE)
          matrix_rho <- cbind(matrix_rho, smap$rho)
          
          m <- m + step_size
        }
        
        matrix_results <- rbind(matrix_results, data.frame(
          w = window_size[i],
          E = E[j],
          tau = tau[k],
          theta = theta[l],
          rho = mean(matrix_rho)
        ))
        
      }
    }
  }
}

results <- matrix_results

### sort results
matrix_rho <- c()
for(i in 1:length(window_size))
{
  best_w <- order(-results[results$w == window_size[i],]$rho)[1]
  matrix_rho <- cbind(matrix_rho, results[results$w == window_size[i],]$rho[best_w])
}


d.t=diff(c(matrix_rho));
ind.w1=which((d.t<median(d.t))&c(matrix_rho[-length(window_size)]>0))
if(length(ind.w1)>0){ind.w2=ind.w1[1]}else{ind.w2=which.max(c(matrix_rho))}
# Select the optimal moving window
opt.w=window_size[ind.w2]

results.w=filter(results,w==opt.w)
# Select the optimal parameter sets
param.opt=results.w[which.max(results.w[,'rho']),]


################################################################################
E <- as.numeric(param.opt['E'])
tau <- as.numeric(param.opt['tau'])
theta <- seq(0,2.5,by=0.5)
window_size <- as.numeric(param.opt['w'])
step_size <- 1

### start algorithm
window_indices <- seq(window_size, NROW(raw_time_series), step_size)
matrix_result <- matrix(NaN, nrow = length(window_indices), ncol = 4)
index <- 0

for(j in window_indices)
{
  index <- index + 1
  rolling_window <- raw_time_series[(j-window_size+1):j]
  
  norm_rolling_window <- (rolling_window - mean(rolling_window, na.rm=TRUE))/sd(rolling_window, na.rm=TRUE)
  
  # calculate best theta
  smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta, silent=TRUE)
  best <- order(-smap$rho)[1]
  theta_best <- smap[best,]$theta
  param.opt['theta'] <- theta_best
  
  # calculate eigenvalues for best theta
  smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta_best, silent=TRUE, save_smap_coefficients=TRUE)
  smap_co <- as.matrix(smap$smap_coefficients[[1]])
  colnames(smap_co)=rownames(smap_co)=NULL
  
  matrix_eigen <- matrix(NA, nrow = NROW(smap_co), ncol = 3)
  
  for(k in 1:NROW(smap_co))
  {
    if(!is.na(smap_co[k,1]))
    {
      M <- rbind(c(smap_co[k, 1:E]), cbind(diag(E - 1), rep(0, E - 1)))
      M_eigen <- eigen(M)$values
      lambda1 <- M_eigen[order(abs(M_eigen))[E]]
      
      matrix_eigen[k,1] <- abs(lambda1)
      matrix_eigen[k,2] <- Re(lambda1)
      matrix_eigen[k,3] <- Im(lambda1)
    }
  }
  
  # save results
  matrix_result[index,1] <- j
  matrix_result[index,2] <- mean(matrix_eigen[,1],na.rm=TRUE)
  matrix_result[index,3] <- mean(matrix_eigen[,2],na.rm=TRUE)
  matrix_result[index,4] <- mean(matrix_eigen[,3],na.rm=TRUE)
}

new_result <- matrix_result

## Plot
### Ricker model
par(mar=c(4.5,4.5,2,1))
plot(new_result[,1], new_result[,2],type = "l",
     col = "blue", las=0, cex.axis=1.25,lwd = 2, cex.lab=1.4, xlab="Time", ylab="|DEV|")
points(c(0,125000),c(1,1), type="l", col="gray", lwd=3, lty=2)
abline(v = 500, col="red", lwd=2, lty=2)
title("Ricker Model")
mtext("D",side=3, at=50, line=0.5,cex=1.5)




