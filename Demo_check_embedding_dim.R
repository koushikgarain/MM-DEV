### Demo Code to check Embedding dimension, window size 

rm(list = ls())

library(rEDM)
library(deSolve)
library(dplyr)

# PART 1 #### read the data set ############################################################################
emp_time_series <- read.csv("Data\\lake_zurich.csv")  # Lake Zurich data

# PART 2 ############# Check best E and window size for each variables #####################################
nvs=ncol(emp_time_series)
result.kk=NULL
for(vb.i in 1:nvs){
  raw_time_series <- emp_time_series[,vb.i]
  ###set up the parameter ranges for performing cross-validation
  E <- seq(1,12,by=1)
  tau <- seq(1,1,by=1)
  theta <- 0
  
  window_size <- seq(25, 250, by=25)
  step_size <- 25
  
  
  
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
  
  ### save E and w for each variables
  result.t=unlist(c(vb.i,param.opt))
  names(result.t)=c('VariableNumber',names(param.opt))
  result.kk=rbind(result.kk,result.t)
}

# PART 3 ############# Check best E and window size ####################################################

E_best=first(result.kk[,'E'])       # we choose first embedding dimension
w_best=first(result.kk[,'w'])      # window size of the first variable

