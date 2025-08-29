########## Figure 3 target variables

rm(list = ls())

library(rEDM)
library(deSolve)
library(dplyr)
library(readxl)
library(vegan)
library(ggplot2)

# load the data
D4_reactor_Navie <- read_excel("Data\\UASB_daughter_reactors_environmental_parameter_table_first110days_20241230update_Navie.xlsx")
D4_reactor_Navie <- as.data.frame(D4_reactor_Navie)
D4_reactor_Navie <- na.omit(D4_reactor_Navie)

EC_data  <- read_excel("Data\\UASB_daughter4_reactors.xlsx")
EC_data  <- as.data.frame(EC_data)

### plot in a column
win.graph(30,50)
par(mfrow=c(2,1),oma=c(0,0,2,0))

# Plot time series for the COD-removal
plot(D4_reactor_Navie[[2]],D4_reactor_Navie[[16]], type='l', col='blue', lwd=2,cex = 1.5,
     cex.axis=1.75, cex.lab=1.5, xlab = "Time (Days)", ylab = paste(colnames(D4_reactor_Navie[16])))
#lines(new_data[[18]],new_data[[17]])
abline(v=80, col="red", lwd=2,lty=2)  # Mark 80
mtext("I",side=3, at=0, line=0.5,cex=1.5)

# Plot time series for the EC
plot(EC_data[[2]],EC_data[[7]], type='l', col='darkgreen', lwd=2,cex = 1.5,
     cex.axis=1.75, cex.lab=1.5, xlab = "Time (Days)", ylab = "Electrical Conductivity (EC)")
#lines(new_data[[18]],new_data[[17]])
abline(v=80, col="red", lwd=2,lty=2)  # Mark 80
mtext("J",side=3, at=0, line=0.5,cex=1.5)


