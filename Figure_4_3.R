####### Figure 4 target variables
rm(list = ls())

# Install and load necessary libraries
#install.packages(c("rEDM", "tidyverse", "zoo", "DescTools", "igraph"))
library(rEDM)
library(tidyverse)
library(zoo)
library(igraph)

# --- STEP 1: Load your data ---
# data loading code
setwd("C:\\Users\\Koushik\\Dropbox\\NTU Papers and ESM\\Multivariate data\\lake zurich")
lake_data <- read.csv("gpedm-regime-shifts-main/data/lake_plankton/lake_zurich.csv")
lake_data <- lake_data[13:516,]
lake_data$Date <- as.Date(paste(lake_data$year, lake_data$month, 1, sep = "-"))
lake_data$Time <- seq(1:nrow(lake_data))
lake_data <- lake_data[, c("Time", setdiff(names(lake_data), "Time"))]

data <- lake_data[,c(1,4:18)]

# Time in first column, 13 plankton groups, phosphate, temperature

# # Interpolate missing environmental data
# data$phosphorus <- na.approx(data$phosphorus)
# data$temperature <- na.approx(data$temperature)

# Scale data (important for CCM)
data_scaled <- as.data.frame(scale(data[ , 2:(ncol(data))]))  # excluding Time column

#time_points <- data$Time

# --- STEP 2: Initialize parameters ---
window_size <- 60  # 5 years
step_size <- 1

guilds <- colnames(data_scaled)[1:13]  # Assuming first 13 columns are plankton groups

guild_type <- c("producer", "producer", "consumer", "consumer", "consumer", "producer", 
                "producer", "consumer", "consumer", "consumer", "consumer", "producer", "consumer")

# --- STEP 3: Moving window CCM analysis ---
connectance_list <- c()
interaction_strength_list <- c()
top_down_count <- c()
bottom_up_count <- c()
window_center <- c()

for (start_idx in seq(1, nrow(data_scaled) - window_size, by = step_size)) {
  sub_data <- data_scaled[start_idx:(start_idx+window_size-1), guilds]
  
  link_count <- 0
  strength_sum <- 0
  total_tests <- 0
  top_down <- 0
  bottom_up <- 0
  
  for (i in 1:(ncol(sub_data)-1)) {
    for (j in (i+1):ncol(sub_data)) {
      ccm_result <- ccm(sub_data, E=3, lib_sizes = seq(10, window_size, by=10),
                        lib_column = guilds[i], target_column = guilds[j], silent=TRUE)
      
      if (max(ccm_result$rho, na.rm=TRUE) > 0.2) {
        link_count <- link_count + 1
        strength_sum <- strength_sum + max(ccm_result$rho, na.rm=TRUE)
        
        if (guild_type[i] == "consumer" & guild_type[j] == "producer") {
          top_down <- top_down + 1
        }
        if (guild_type[i] == "producer" & guild_type[j] == "consumer") {
          bottom_up <- bottom_up + 1
        }
      }
      total_tests <- total_tests + 1
    }
  }
  
  connectance_list <- c(connectance_list, link_count / total_tests)
  interaction_strength_list <- c(interaction_strength_list, strength_sum / max(1, link_count))
  top_down_count <- c(top_down_count, top_down)
  bottom_up_count <- c(bottom_up_count, bottom_up)
  window_center <- c(window_center, data$Time[start_idx + window_size/2])
}

transition_data1 <- cbind.data.frame(connectance_list,interaction_strength_list,top_down_count,bottom_up_count,window_center)
transition_data <- cbind (transition_data1,lake_data[(window_size+1):nrow(lake_data),'Date'])

#transition_data <- data.frame(transition_data)
# Plot results ---
#par(mfrow=c(3,1),oma=c(0,0,2,0))

### plot in a column
win.graph(40,90)
par(mfrow=c(3,1),oma=c(0,0,2,0))
###################### Connectance over time
plot(transition_data[,6], transition_data[,1], type='l', col='blue', lwd=2, cex = 1.5,
     cex.axis=1.75, cex.lab=1.5,xlim = c(as.Date("1978-01-01"), as.Date("2019-12-01")),
     ylab="Connectance", xlab="Time", main="")
abline(v=as.Date("1996-01-01"), col="red",lwd=2, lty=2)  # Mark 1996
mtext("F",side=3, at=as.Date("1978-01-01"), line=0.5,cex=1.5)

############# Top-down vs Bottom-up link counts
plot(transition_data[,6], transition_data[,3]-transition_data[,4], type='l', col='blue', lwd=2, cex = 1.5,
     cex.axis=1.75, cex.lab=1.5,xlim = c(as.Date("1978-01-01"), as.Date("2019-12-01")),
     ylab="TD - BU Count", xlab="Time", main="")
#abline(h=0, lty=2)
abline(v=as.Date("1996-01-01"), col="red", lwd=2, lty=2)  # Mark 1996
mtext("G",side=3, at=as.Date("1978-01-01"), line=0.5,cex=1.5)

######### Phoporus time series
data <- read.csv("gpedm-regime-shifts-main/data/lake_plankton/lake_zurich.csv")

# Convert Year and Month into Date format (assuming day = 1)
new_data <- data[13:516,]
new_data$Date <- as.Date(paste(new_data$year, new_data$month, 1, sep = "-"))
 
# Plot time series for the variable
plot(new_data[[18]],new_data[[17]], type='l', col='blue', lwd=2,cex = 1.5,
     cex.axis=1.75, cex.lab=1.5, xlab = "Time", ylab = 'Phosphorus')
#lines(new_data[[18]],new_data[[17]])
abline(v=as.Date("1996-01-01"), col="red", lwd=2,lty=2)  # Mark 1996
mtext("H",side=3, at=as.Date("1978-01-01"), line=0.5,cex=1.5)



