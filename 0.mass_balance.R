## this file was used to validate the PBPK model through the mass balance checking,
## the model was built based on the mrgsolve package, the model code was saved in the mousePBPK.code file,
## the mass balance was calculated by summing the mass of the organs and the mass of the blood, 
## then comparing the total mass with the initial mass of the dose


## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(gridExtra)

Obs.A1 <- read.csv(file ="dataset/tk/mouse/gold/R_input_mouse_study1_13nm_short.csv")  


#---------------------Build mrgsolve-based PBPK Model-------
mod <- mcode ("mouse_PBPK", mousePBPK.code)

set.seed(5)

#--------------------1. model mass balance checking---------------

## Define the exposure scenario

BW           = 0.02                              ## kg, body weight
tinterval    = 0.1                               ## hr, Time interval
TDoses       = 1                                 ## Dose times, only one dose
PDOSE        = 0.85                              ## mg/kg-day, Single dose
DOSE         = PDOSE*BW                          ## mg, amount of iv dose
ex.iv<- ev(ID=1, amt= DOSE,                      ## Set up the exposure events
           ii=tinterval, addl=TDoses-1, 
           cmt="MBV", replicate = FALSE) 

## Set up the exposure time
tsamp=tgrid(0,max(Obs.A1$Time),tinterval)     ## Simulation time 24*7 hours (180 days)

## calculate the deposition volume
out <- 
  mod %>% 
  update(atol = 1E-80,maxsteps = 5000000) %>%
  mrgsim_d(data = ex.iv, tgrid=tsamp)

## save the calculated into data frame
out <- data.frame(Time=out$time, 
                  Qbal = out$Q,
                  CL=out$Liver,
                  CLt = out$Liver_t,
                  CK = out$Kidney,
                  CKt = out$Kidney_t,
                  CS = out$Spleen,
                  CSt = out$Spleen_t,
                  CR=out$Rest,
                  Clung = out$Lung,
                  CLungt = out$Lung_t,
                  mtot = out$M_tot,
                  iv = out$iv
)

mean(out$mtot) # check the mass balance 

# first result plotting
plot (y=out$mtot, x=out$Time,xlab = "Time", main = paste("Lung Tissue concentration"),
      type='l')

nm_size = "13nm"

plot (y=out$CL, x=out$Time,ylab = "Concentration in Liver (ng/g)",
      xlab="Time (h)", type='l')
plot (y=out$CLt, x=out$Time,xlab = "Time", main = paste("Liver Tissue concentration",nm_size),
      type='l',log="x")
plot (y=out$CSt, x=out$Time,xlab = "Time", main = paste("Spleen Tissue concentration",nm_size),
      type='l',log="x",xlim=c(0.1,1000),ylim=c(0,50000))
plot (y=out$CKt, x=out$Time,xlab = "Time", main = paste("Kidney Tissue concentration",nm_size),
      type='l',log="x",xlim=c(0.1,1000),ylim=c(0,1500))
plot (y=out$CLungt, x=out$Time,xlab = "Time", main = paste("Lung Tissue concentration",nm_size),
      type='l',log="x",xlim=c(0.1,1000),ylim=c(0,1500))

