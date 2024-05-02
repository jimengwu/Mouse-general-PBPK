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

source("helper_functions.R")
source("Mouse_PBPK.R")
#---------------------------13 nm full dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_13nm.csv")  
# LONG TERM (WITHIN 24h)
Obs.A1 <- Obs.A1[1:7,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/13nm_nsc_GI/'


#----------------------ZnO-----------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/ZnO/1-ZnO.csv") 
Obs.A1 <- Obs.A1[1:3,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/ZnO/1_GI/'



#------------------Study2: 13 nm high dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/31501470_Gold dextran.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/31501470_nsc_GI/'
#---------------------------Study2: 13nm low dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/31501470_Gold dextran_2.csv")  
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/31501470_2_nsc_GI/'


#--------------------------Study 3: 61.2nm------------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/29677597_PEG coated Au NPs 61.2nm.csv") 
Obs.A1 <- Obs.A1[1:8,]
PDOSE = Obs.A1$Dose..mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/29677597_61.2_nsc_GI/'
#--------------------------Study 3: 24.3nm------------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/29677597_PEG coated Au NPs 24.3nm.csv") 
Obs.A1 <- Obs.A1[1:8,]
PDOSE = Obs.A1$Dose..mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/29677597_24.3_nsc_GI/'


#---------------------------100 nm short dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_100nm.csv")  
# SHORT TERM (WITHIN 24h)
Obs.A1 <- Obs.A1[1:3,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","Clung")] #100nm only three organs
folder = 'plots/100nm_short_nsc_GI/'



#------------------GO: Study1: 20 nm high dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/GO/21162527_125I-NGS-PEG 10-30nm.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/GO/2116257_GI/'



#------------------TiO2: Study1: 20 nm high dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/TiO2/1-TiO2.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/TiO2/1_GI/'


#-------------------------done-------



#---------------------to be improved--------












#---------------------------13 nm short dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_13nm.csv")  
# SHORT TERM (WITHIN 24h)
Obs.A1 <- Obs.A1[1:3,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/13nm_short_nsc_2/'


#--------------------------Study 3: 6.2nm short------------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/29677597_PEG_coated_AuNPs_6.2nm.csv") 
Obs.A1 <- Obs.A1[1:5,] # due to the data missing. the time range is limited here
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/29677597_6.2_short_nsc/'






#---------------------------4 nm short dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_4nm.csv")  
# SHORT TERM (WITHIN 24h)
Obs.A1 <- Obs.A1[1:3,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/4nm_short_nsc/'




#--------------------------Study 3: 6.2nm time------------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/29677597_PEG_coated_AuNPs_6.2nm.csv") 
Obs.A1 <- Obs.A1[1:8,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/29677597_6.2_nsc/'


#--------------------------Study 3: 42.5nm------------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/29677597_PEG coated Au NPs 42.5nm.csv") 
Obs.A1 <- Obs.A1[1:8,]
PDOSE = Obs.A1$Dose..mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/29677597_42.5_nsc_GI/'



#---------------------------4 nm full dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_4nm.csv")  
# SHORT TERM (WITHIN 24h)
Obs.A1 <- Obs.A1[1:7,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/4nm_nsc_GI/'


#---------------------------100 nm full dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_100nm.csv")  
Obs.A1 <- Obs.A1[1:7,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","Clung")] #100nm only three organs
folder = 'plots/100nm_nsc_GI/'



#---------------------Build mrgsolve-based PBPK Model-------
mod <- mcode ("mouse_PBPK", mousePBPK.code)

set.seed(5)
## 2.1 initial parameters-----------
params.init <- log(c(
  K_release_Liver = 0.001,  # h-1
  K_max_Liver = 20,         # h-1
  A_cap_liver = 1000,
  #K_50_Liver = 48,          # h
  #n_Liver = 5,              # Unitless
  K_release_GI = 0.001,     # h-1
  K_max_GI = 0.075,         # h-1
  K_GI_b = 4e-5,
  #K_50_GI = 24,             # h
  #n_GI = 5,                 # Unitless
  K_release_Spleen = 0.001, # h-1
  K_max_Spleen = 40,        # h-1
  #A_cap_Spleen = 100,
  #K_50_Spleen = 48,
  #n_Spleen = 5,
  K_release_Kidney = 0.0004, # h-1
  K_max_Kidney = 0.075,
  #A_cap_kidney = 100,
  #K_50_Kidney = 24,
  #n_Kidney = 5,
  K_release_Lung = 0.003,   # h-1
  K_max_Lung = 0.075,
  #K_50_Lung = 24,
  #n_Lung = 5,               
  P_Liver  = 0.08,
  P_Brain  = 0.15,
  P_Kidney  = 0.15,
  P_Spleen  = 0.15,
  P_Lung  = 0.15,
  P_Rest  = 0.15,
  P_GI = 0.15,
  DLC_Liver = 0.001,
  DLC_Brain = 0.000001,
  DLC_Kidney = 0.001,
  DLC_Spleen = 0.03,
  DLC_Lung = 0.001,
  DLC_Rest = 0.000001,
  DLC_GI = 0.001,
  Kbile = 0.00003,       # Biliary clearance (L/hr)
  Kurine = 0.000003,   # Urine clearance (L/hr)
  Kfecal = 0.000003
))


tstep = min(1,min(Obs.A1$Time))


pred.mouse <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)
  
  ## Define the exposure scenario
  
  BW           = 0.02                              ## kg, body weight
  tinterval    = 1                                 ## hr, Time interval
  TDoses       = 1                                 ## Dose times, only one dose
  #PDOSE        = 0.85                              ## mg/kg-day, Single dose
  DOSE         = PDOSE*BW                          ## mg, amount of iv dose
  ex.iv<- ev(ID=1, amt= DOSE,                  ## Set up the exposure events
             ii=tinterval, addl=TDoses-1, 
             cmt="MBV", replicate = FALSE) 
  
  ## Set up the exposure time
  tsamp=tgrid(0,max(Obs.A1$Time),tstep)     ## Simulation time 24*7 hours (180 days)
  
  ## calculate the deposition volume
  out <- 
    mod %>% 
    param(pars) %>%
    ##Req(Liver,M_tot,MBV)%>%d
    update(atol=1e-50,maxsteps = 500000000) %>%
    mrgsim_d(data = ex.iv, tgrid=tsamp)
  
  ## save the calculated into data frame
  out <- data.frame(Time=out$time, 
                    CL=out$Liver_t,                    
                    CS = out$Spleen_t,
                    CK = out$Kidney_t,
                    Clung = out$Lung_t)
  
  return(out)
}

#-----oral
pred.mouse <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)
  
  ## Define the exposure scenario
  
  BW           = 0.02                              ## kg, body weight
  tinterval    = 1                                 ## hr, Time interval
  TDoses       = 1                                 ## Dose times, only one dose
  #PDOSE        = 0.85                              ## mg/kg-day, Single dose
  DOSE         = PDOSE*BW                          ## mg, amount of iv dose
  ex.iv<- ev(ID=1, amt= DOSE,                  ## Set up the exposure events
             ii=tinterval, addl=TDoses-1, 
             cmt="M_GI_lumen", replicate = FALSE) 
  
  ## Set up the exposure time
  tsamp=tgrid(0,max(Obs.A1$Time),tstep)     ## Simulation time 24*7 hours (180 days)
  
  ## calculate the deposition volume
  out <- 
    mod %>% 
    param(pars) %>%
    ##Req(Liver,M_tot,MBV)%>%d
    update(atol=1e-50,maxsteps = 500000000) %>%
    mrgsim_d(data = ex.iv, tgrid=tsamp)
  
  ## save the calculated into data frame
  out <- data.frame(Time=out$time, 
                    CL=out$Liver_t,
                    CK = out$Kidney_t,
                    CS = out$Spleen_t,
                    Clung = out$Lung_t)
  
  return(out)
}


########################################## 4. Model Evalution with MCMC  ###########################################
#################################################################################################################

paras = read.csv(paste0(folder,'mod_fit/params_fitted_v2.csv'))
params.mod = as.list(paras)$x
names(params.mod)<- names(params.init)

Fitted_output.A1 = pred.mouse(par=params.mod)
#------Input parameters set----------

# Define parameter names
sig_names <- c(
  "sig2",
  "sig_K_release_Liver", "sig_K_max_Liver", 
  "sig_A_cap_Liver", 
  "sig_K_release_GI", "sig_K_max_GI", "sig_K_GI_b", 
  "sig_K_release_Spleen", "sig_K_max_Spleen", 
  "sig_K_release_Kidney", "sig_K_max_Kidney", 
  "sig_K_release_Lung", "sig_K_max_Lung", 
  "sig_P_Liver", "sig_P_Brain", "sig_P_Kidney", 
  "sig_P_Spleen", "sig_P_Lung", "sig_P_Rest","sig_P_GI",
  "sig_DLC_Liver", "sig_DLC_Brain", "sig_DLC_Kidney", 
  "sig_DLC_Spleen", "sig_DLC_Lung", "sig_DLC_Rest","sig_DLC_GI",
  "sig_Kbile", "sig_Kurine","sig_Kfecal"
)

sig_population = 0.5

# Create a named vector with default values
sig_values <- rep(sig_population, length(sig_names))
names(sig_values) <- sig_names

# Update specific values
sig_pred = 0.3
sig_values["sig2"] <- sig_pred

# Create sig_list
sig_list <- log(sig_values)

theta.MCMC<- c(params.mod, sig_list)

# 13nm
if (!file.exists(paste0(folder,"MCMC/"))) {
  dir.create(paste0(folder,"MCMC/"), recursive = TRUE)
}


saveRDS(theta.MCMC,file =paste0(folder,'MCMC/theta.rds'))

which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG


##-------------------------prior distribution plotting--------------

# Mean and SD
mean.prior = exp(theta.MCMC[-which_sig])
sd.prior   = mean.prior*sig_population

a.prior = qnorm(0.025, mean=mean.prior,sd = sd.prior)
b.prior = qnorm(0.975, mean=mean.prior,sd=sd.prior)
c.prior = qnorm(0.5, mean=mean.prior, sd=sd.prior)

# Initialize lists to store samples for each parameter
param_samples <- list()

# Number of points to use in density estimation
n_points <- 1000  # Adjust as needed for smoother or rougher curves


for (param_index in seq_along(mean.prior)) {
  # Draw samples from the normal distribution for the current parameter
  param_samples[[param_index]] <- rnorm(n_points, mean = mean.prior[[param_index]], sd = sd.prior[[param_index]])
}

library(ggplot2)

# Initialize an empty list to store density estimation dataframes for each parameter
density_data <- list()

# Generate density estimation and store it for each parameter
for (param_index in seq_along(mean.prior)) {
  # Generate density estimation for the current parameter
  density_estimation <- density(param_samples[[param_index]], n = n_points)
  
  # Create a dataframe for the density estimation data
  density_df <- data.frame(x = density_estimation$x, density = density_estimation$y, 
                           Par = names(mean.prior)[param_index])
  
  # Store the density dataframe
  density_data[[param_index]] <- density_df
}

# Combine all density dataframes into one dataframe
combined_density_df <- do.call(rbind, density_data)


# Calculate mean and standard deviation for each parameter
pri_mean_sd_df <- data.frame(Par = names(mean.prior),
                             Mean = unlist(mean.prior),
                             SD = unlist(sd.prior))


# Plot density curves with facets for each parameter
pri_den_plot<- ggplot(combined_density_df, aes(x = x, y = density, fill = "Density")) +
  geom_area(color="black",size=1) +
  #geom_area(fill = "gray70") +  # Fill the area under the curve with gray color
  facet_wrap(~ Par, scales = "free", ncol = 6) +
  labs(title = "Density Distribution by Parameter") +
  theme_minimal()+
  theme(panel.grid = element_blank(),  # Remove grid lines
        panel.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 30, hjust = 1)) +  # Set background to transparent
  geom_text(data = pri_mean_sd_df, aes(label = paste0("Mean:", sprintf("%.1e", Mean), 
                                                      "\nSD:", sprintf("%.1e", SD))), 
            x = Inf, y = -Inf, hjust = 1, vjust = -1, size = 3.5,color="black") +  # Add mean and sd as text annotations
  scale_fill_manual(values = "gray")


# Save the plot
ggsave(paste0(folder,"mc_sens/pri_density_distribution_plot.png"), pri_den_plot,
       width = 30, height = 20, units = "cm")

# Save to CSV
write.csv(pri_mean_sd_df, file = paste0(folder,"mc_sens/pri_paras_stats.csv"), row.names = FALSE)
##---------------------------prior distribution plotting--------------


#--------------Maximum likelihood estimation (MLE) function for MCMC----
mcmc.fun <- function (pars){
  
  out <- pred.mouse(pars[-which_sig])
  out = out[which(out$Time %in% Obs.A1$Time),]
  # making sure the predicted value and observed value matching with each other as
  # the time step in predicting function might not be small enough to cover the time step in the 
  # observed time step
  Obs.A1 = Obs.A1[which(Obs.A1$Time %in% out$Time),] 
  
  ## log-transformed prediction
  log.yhat.CL     <- log(out$CL)
  #log.yhat.CK     <- log(out$CK)
  log.yhat.CS     <- log(out$CS)
  log.yhat.Clung  <- log(out$Clung)
  
  
  ## log-transformed experimental data
  
  log.y.CL        <- log(Obs.A1$CL)
  #log.y.CK        <- log(Obs.A1$CK)
  log.y.CS        <- log(Obs.A1$CS)
  log.y.Clung     <- log(Obs.A1$Clung)
  
  
  ## The method of Maximum likelihood
  #log.yhat        <- c(log.yhat.CL,log.yhat.CK,log.yhat.CS,log.yhat.Clung)
  #log.y           <- c(log.y.CL,log.y.CK,log.y.CS,log.y.Clung)
  
  log.yhat        <- c(log.yhat.CL,log.yhat.CS,log.yhat.Clung)
  log.y           <- c(log.y.CL,log.y.CS,log.y.Clung)
  

  
  #todo：100nm
  #log.yhat        <- c(log.yhat.CL)
  #log.y           <- c(log.y.CL)
  
  non_nan_indices <- which(!is.na(log.y))
  log.yhat        <- log.yhat[non_nan_indices]
  log.y           <- log.y[non_nan_indices]
  
  sig2            <- as.numeric((exp(pars[which_sig][1])))
  
  log_likelihood  <- -2*sum((dnorm (log.y,
                                    mean = log.yhat,
                                    sd   = sqrt(sig2), 
                                    log=TRUE)))
  
  return(log_likelihood)
  
}


## Define the Prior distributions: either normal or log normal distribution
## normal distribution
Prior <- function(pars) {
  
  ## Population level
  # The likelihood for population mean (parameters)
  # Calculate likelihoods of each parameters; P(u|mean,CV)
  pars.data = exp(pars[-which_sig]) # parameter value
  
  mean           = exp(theta.MCMC[-which_sig]) # parameter mean in the distribution
  CV             = 0.3     # Coefficient of variation; Default value of 0.5 in all parameters (Bois,2000; Bois et al., 1996)
  
  prior_pars     = dtruncnorm(pars.data, 
                              a = qnorm(0.025, mean = mean, sd = mean*CV), 
                              b = qnorm(0.975, mean = mean, sd = mean*CV), 
                              mean = mean, sd = mean*CV ) 
  
  # The likelihood for population variance; P(sigmal^2|sigmal0^2)
  
  CV.sig         = exp(theta.MCMC[which_sig])[2:length(theta.MCMC[which_sig])]               # Sigmal0
  
  CU             = 1                                              # Coefficient of uncertainty (CU) (Hack et al., 2006)
  alpha          = (2+1)/(CU^2)                                   # Shape parameter  of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  beta           = (alpha-1)*CV.sig^2                             # Scale parameter  of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  
  # Calculate likelihoods of model error (sig2) and population variance (sig) parameters
  
  sig  <- as.numeric (exp(pars[which_sig][2:length(pars.data [which_sig])]))   # Coefficient of variation from population variance; sigmal0
  prior_sig      = dinvgamma (sig, shape = alpha , scale = beta)  # prior distribution for population variance; sigma2
  
  
  ## individual level; P(theta|u,sigmal^2)
  mean_i         = prior_pars
  sd_i           = sqrt(prior_sig)
  prior_pars_i   = dtruncnorm (prior_pars, 
                               a = qnorm(0.025, mean = mean_i, sd = sd_i), 
                               b = qnorm(0.975, mean = mean_i, sd = sd_i), 
                               mean = mean_i, sd = sd_i) 
  
  # model residuals
  
  sig2 <- as.numeric (exp(pars[which_sig][1]))                    # error variances from model residual
  prior_sig2     = dunif (sig2, min = 0.01, max = 3)            # error variances, Lower and upper boundary from Chiu et al., 2009; Chiu et al., 2014)   
  
  
  # log-transformed (log-likelihoods of each parameters)
  log.pri.pars   = log (prior_pars)
  log.pri.sig    = log (prior_sig)
  log.pri.pars.i = log (prior_pars_i)
  log.pri.sig2   = log (prior_sig2)
  
  # maximum likelihood estimation (MLE): negative log-likelihood function, (-2 times sum of log-likelihoods)
  #MLE =  -2*sum(log.pri.pars, log.pri.sig2)  
  MLE =  -2*sum(log.pri.pars, log.pri.sig , log.pri.pars.i,log.pri.sig2)  
  
  return(MLE)
}



#################### 5. MCMC simulation with parallel computing ############################
detectCores()                                ## check the cores
cl<- makeCluster(detectCores())              ## use all cores in our system     
registerDoParallel(cl)                       ## registers a cluster of all the cores on our system

# start time
tstr<-Sys.time()
#niter = 120000
#burninlength  = 60000
#outputlength  = 60000
niter = 200
burninlength  = 100
outputlength  = 100


niter = 300000
burninlength  = 150000
outputlength  = 150000


# parallel
system.time(
  MCMC <- foreach( i = 1:4, .packages = c('mrgsolve','magrittr','FME',
                                          'truncnorm','EnvStats',
                                          'invgamma','dplyr')) %dopar% {
                                            mod <- mod
                                            modMCMC(f             = mcmc.fun, 
                                                    p             = theta.MCMC, 
                                                    niter = niter,
                                                    jump          = 0.01,             ## jump function generation new parameters distribution using covrate matrix
                                                    prior         = Prior,            ## prior function
                                                    updatecov     = 100,               ## adaptive Metropolis
                                                    #wvar0         = 0.01,             ## "weight" for the initial model variance
                                                    ntrydr        = 5,                ## delayed Rejection
                                                    burninlength  = burninlength,    ## number of initial iterations to be removed from output.
                                                    outputlength  = outputlength,     ## number of output iterations  
                                                    verbose=1
                                            )                    
                                            
                                          }
)
tend<-Sys.time()

#end time
print(tend-tstr)

stopCluster(cl)   


## Performance four chains to check the convergences
MC.mouse.1 = as.mcmc (MCMC[[1]]$pars) # first  chain
MC.mouse.2 = as.mcmc (MCMC[[2]]$pars) # second chain
MC.mouse.3 = as.mcmc (MCMC[[3]]$pars) # third  chain
MC.mouse.4 = as.mcmc (MCMC[[4]]$pars) # fourth chain

## combine all chains
combinedchains = mcmc.list(MC.mouse.1,MC.mouse.2,MC.mouse.3,MC.mouse.4) 

gelman.diag(combinedchains)

mcmc_trace (
  combinedchains,
  pars =names(theta.MCMC[1:length(theta.MCMC [-which_sig])]),
  size = 0.5,
  facet_args = list(nrow = 5)) +
  ggplot2::scale_color_brewer()


# Create the directory if it doesn't exist
if (!file.exists(paste0(folder,"MCMC/"))) {
  dir.create(paste0(folder,"MCMC/"), recursive = TRUE)
}

# output the MCMC results

write.csv(MC.mouse.1,file=paste0(folder,"MCMC/mouse.chain1.csv"))
saveRDS(MCMC,file =paste0(folder,'MCMC/mouse.MCMC.rds'))
saveRDS(combinedchains,file=paste0(folder,'MCMC/mouse.MCMC.comb.rds'))

#--------------------------------------plots-------------------------
#--------------------1. Densities of posterior parameter uncertainty distributions of the population mean (μ). ----------


Mouse.MCMC        <- readRDS(file = paste0(folder,"MCMC/mouse.MCMC.rds"))[[1]] # use only first chain

## Save the posterior parameters (95% CI)
quan.mouse = exp(summary(MC.mouse.1)$quantiles)  
write.csv(quan.mouse,file=paste0(folder,"MCMC/mouse.summary_pos.csv"))
# Calculate mean and standard deviation for each parameter
post_mean_sd_df <- summary(MC.mouse.1)$statistics

## loading the theta names
theta.names       <- names(readRDS(file = paste0(folder,"MCMC/theta.rds")))

theta             <- readRDS(file = paste0(folder,"MCMC/theta.rds"))

## Sampling from posterior parameters to generate the posterior distributions 

## Mouse posteiror distributions
M.Mouse  <- exp(Mouse.MCMC$pars) %>% apply(2,mean)
SD.Mouse <- exp(Mouse.MCMC$pars) %>% apply(2,sd)

# Initialize lists to store samples for each parameter
post_param_samples <- list()

# Number of points to use in density estimation
n_points <- 5000  # Adjust as needed for smoother or rougher curves


for (param_index in seq_along(M.Mouse)) {
  # Draw samples from the normal distribution for the current parameter
  post_param_samples[[param_index]] <- rnorm(n_points, 
                                             mean = M.Mouse[[param_index]], 
                                             sd = SD.Mouse[[param_index]])
}

# Initialize an empty list to store density estimation dataframes for each parameter
post_density_data <- list()

# Generate density estimation and store it for each parameter
for (param_index in seq_along(M.Mouse)) {
  # Generate density estimation for the current parameter
  post_density_estimation <- density(post_param_samples[[param_index]], n = n_points)
  
  # Create a dataframe for the density estimation data
  post_density_df <- data.frame(x = post_density_estimation$x, 
                                density = post_density_estimation$y, 
                                Par = names(M.Mouse)[param_index])
  
  # Store the density dataframe
  post_density_data[[param_index]] <-post_density_df
}

# Combine all density dataframes into one dataframe
combined_post_density_df <- do.call(rbind, post_density_data)

# Calculate mean and standard deviation for each parameter
post_mean_sd_df <- data.frame(Par = names(M.Mouse),
                              Mean = unlist(M.Mouse),
                              SD = unlist(SD.Mouse))


# Plot density curves with facets for each parameter the nrow and ncol might need to be changed

library(ggplot2)


post_den_plot <- ggplot(combined_post_density_df %>%
                          filter(!grepl("sig", Par)), aes(x = x, y = density)) +
  geom_line(color="#87BFD6") +
  facet_wrap(~ Par, scales = "free", nrow = 7, ncol = 5) +
  labs(title = "Posterior Density Distribution by Parameter") +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid lines
        panel.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 30, hjust = 1)) + #"#f7f7f7"
  geom_text(data = post_mean_sd_df %>%
              filter(!grepl("sig", Par)), 
            aes(label = paste0("Mean:", sprintf("%.2e", Mean), "\nSD:", sprintf("%.2e", SD))), 
            x = Inf, y = Inf, hjust = 1, vjust = 1.5, size = 3.5, color = "black")  # Add mean and sd as text annotations

#print(post_den_plot)


prior_post_den_plot <- ggplot() +
  geom_area(data = combined_density_df, aes(x = x, y = density, fill = "Prior"), color = "black", alpha = 0.7) +
  geom_area(data = combined_post_density_df %>% filter(!grepl("sig", Par)), 
            aes(x = x, y = density, fill = "Posterior"), color = "black", alpha = 0.7) +
  facet_wrap(~ Par, scales = "free", nrow = 7, ncol = 6) +
  labs(title = "Prior Density and Posterior Density Distribution by Parameter") +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid lines
        panel.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 30, hjust = 1)) + #"#f7f7f7"
  geom_text(data = post_mean_sd_df %>% filter(!grepl("sig", Par)), 
            aes(label = paste0("Mean:", sprintf("%.1e", Mean), "\nSD:", sprintf("%.1e", SD))), 
            x = Inf, y = Inf, hjust = 1, vjust = 1.5, size = 3, color = "#D17E5E")  + # Add text annotations for combined_post_density_df
  geom_text(data = pri_mean_sd_df %>% filter(!grepl("sig", Par)), 
            aes(label = paste0("Mean:", sprintf("%.1e", Mean), "\nSD:", sprintf("%.1e", SD))), 
            x = Inf, y = Inf, hjust = 1, vjust = 3, size = 3, color = "#9EC7C5") + # Add text annotations for combined_density_df
  scale_fill_manual(name = "Distribution", values = c("Prior" = "#9EC7C5", "Posterior" = "#D17E5E")) + # Add legend for distribution
  guides(fill = guide_legend(title = "Distribution")) # Add legend title

#print(prior_post_den_plot)



# Save the plot
ggsave(paste0(folder,"mc_sens/post_density_distribution_plot.png"), post_den_plot, 
       width = 30, height = 16, units = "cm")

ggsave(paste0(folder,"mc_sens/prior_post_density_distribution_plot.png"), prior_post_den_plot, 
       width = 30, height = 16, units = "cm")


# Save to CSV
write.csv(post_mean_sd_df, file = paste0(folder,"mc_sens/post_paras_stats.csv"), row.names = FALSE)
#-----------------------------2. get sensitivity range plots-------------

Mouse.MCMC        <- readRDS(file = paste0(folder,"MCMC/mouse.MCMC.rds"))[[1]] # only first

#Mouse.MCMC <- MCMC[[1]]
Newtime.r   = pred.mouse(theta.MCMC)$Time  

# this is the new time variable, now it has been changed to sample per day.
nrwo.r = length (Newtime.r)


MC.mouse.a.CL    = matrix(nrow = nrwo.r, ncol = outputlength/100)
MC.mouse.a.CK    = matrix(nrow = nrwo.r, ncol = outputlength/100)
MC.mouse.a.Clung = matrix(nrow = nrwo.r, ncol = outputlength/100)
MC.mouse.a.CS    = matrix(nrow = nrwo.r, ncol = outputlength/100)

for(i in 1:(outputlength/100)){
  
  j = i *100  # sample parameter set once every ten sets, 
  #so you will have 5000 sets from 50000 total sets
  
  pars.mouse             = Mouse.MCMC$pars    [j,]     
  MCdata               = pred.mouse (pars.mouse)
  MC.mouse.a.CL[,i]   = MCdata$CL
  MC.mouse.a.CK[,i]   = MCdata$CK
  MC.mouse.a.Clung[,i]= MCdata$Clung
  MC.mouse.a.CS[,i]= MCdata$CS
  cat("iteration = ", i , "\n") # Shows the progress of iterations, so you can see the number of iterations that has been completed, and how many left.
}

M.mouse.a.CL  = MC.mouse.a.CL %>% apply(1,mean) # get the mean value of MC.rat.a.CA
SD.mouse.a.CL = MC.mouse.a.CL %>% apply(1,sd)
M.mouse.a.CK  = MC.mouse.a.CK %>% apply(1,mean) # get the mean value of MC.rat.a.CA
SD.mouse.a.CK = MC.mouse.a.CK %>% apply(1,sd)
M.mouse.a.Clung  = MC.mouse.a.Clung %>% apply(1,mean) # get the mean value of MC.rat.a.CA
SD.mouse.a.Clung = MC.mouse.a.Clung %>% apply(1,sd)
M.mouse.a.CS  = MC.mouse.a.CS %>% apply(1,mean) # get the mean value of MC.rat.a.CA
SD.mouse.a.CS = MC.mouse.a.CS %>% apply(1,sd)



MC.mouse.CL.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.mouse.a.CL, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)
MC.mouse.CK.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.mouse.a.CK, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)
MC.mouse.Clung.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.mouse.a.Clung, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)


MC.mouse.CS.plot <- cbind(
  Time = Newtime.r, 
  as.data.frame(t(apply(MC.mouse.a.CS, 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.1, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.9, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)



p.r.L <- 
  ggplot() + 
  geom_ribbon(data = MC.mouse.CL.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est, color = "95% CI"), 
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(data = MC.mouse.CL.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3, color = "50% CI"), 
              fill = "mediumblue", alpha = 0.3) +
  geom_line(data = MC.mouse.CL.plot, aes(x = Time, y = median_est, color = "Median"), 
            linewidth = rel(1)) +
  geom_point(data = Obs.A1, aes(x = Time, y = CL, color = "Observed"), 
             shape = 1, fill = "white", size = 3, stroke = 2) +
  scale_color_manual(values = c("95% CI" = alpha("lightblue", 0.3), 
                                "50% CI" = alpha("mediumblue", 0.3), 
                                "Median" = "blue", 
                                "Observed" = "red"), 
                     labels = c("95% CI", "50% CI", "Median", "Observed")) + 
  ylab("Concentration in Liver (ng/g)") + xlab("Time (h)") +
  theme_minimal()+
  theme(axis.title = element_text(size = 22,color="black"),  # Adjust the size of x-axis label
        axis.text = element_text(size = 20,color="black"),
        panel.grid.major = element_blank(),    # Remove major grid lines
        panel.grid.minor = element_blank(),     # Remove minor grid lines
        axis.line = element_line(color = "black")         # Show x and y axes
  )

p.r.L

#ggsave(paste0(folder,"MCMC/MCMC_fit_liver.png"),width = 10, height = 8)


p.r.K <- 
  ggplot() + 
  
  geom_ribbon(data = MC.mouse.CK.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est, color = "95% CI"), 
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(data = MC.mouse.CK.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3, color = "50% CI"), 
              fill = "mediumblue", alpha = 0.3) +
  geom_line(data = MC.mouse.CK.plot, aes(x = Time, y = median_est, color = "Median"), 
            size = rel(1)) +
  geom_point(data = Obs.A1, aes(x = Time, y = CK, color = "Observed"), 
             shape = 1, fill = "white", size = 3, stroke = 2) +
  scale_color_manual(values = c("95% CI" = alpha("lightblue", 0.3), 
                                "50% CI" = alpha("mediumblue", 0.3), 
                                "Median" = "blue", 
                                "Observed" = "red"), 
                     labels = c("95% CI", "50% CI", "Median", "Observed")) + 
  
  ylab("Concentration in Kidney (ng/g)") + xlab("Time (h)") +
  theme_minimal()+
  theme(axis.title = element_text(size = 22,color="black"),  # Adjust the size of x-axis label
        axis.text = element_text(size = 20,color="black"),
        panel.grid.major = element_blank(),    # Remove major grid lines
        panel.grid.minor = element_blank(),     # Remove minor grid lines
        axis.line = element_line(color = "black")         # Show x and y axes
  )

#p.r.K
#ggsave(paste0(folder,"MCMC/MCMC_fit_Kidney.png"),width = 10, height = 8)


p.r.lung <- 
  ggplot() + 
  
  geom_ribbon(data = MC.mouse.Clung.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est, color = "95% CI"), 
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(data = MC.mouse.Clung.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3, color = "50% CI"), 
              fill = "mediumblue", alpha = 0.3) +
  geom_line(data = MC.mouse.Clung.plot, aes(x = Time, y = median_est, color = "Median"), 
            size = rel(1)) +
  geom_point(data = Obs.A1, aes(x = Time, y = Clung, color = "Observed"), 
             shape = 1, fill = "white", size = 3, stroke = 2) +
  scale_color_manual(values = c("95% CI" = alpha("lightblue", 0.3), 
                                "50% CI" = alpha("mediumblue", 0.3), 
                                "Median" = "blue", 
                                "Observed" = "red"), 
                     labels = c("95% CI", "50% CI", "Median", "Observed")) + 
  
  ylab("Concentration in Lung (ng/g)") + xlab("Time (h)") +
  theme_minimal()+
  theme(axis.title = element_text(size = 22,color="black"),  # Adjust the size of x-axis label
        axis.text = element_text(size = 20,color="black"),
        panel.grid.major = element_blank(),    # Remove major grid lines
        panel.grid.minor = element_blank(),     # Remove minor grid lines
        axis.line = element_line(color = "black")         # Show x and y axes
  )

#p.r.lung
#ggsave(paste0(folder,"MCMC/MCMC_fit_lung.png"),width = 10, height = 8)


p.r.S <- 
  ggplot() + 
  geom_ribbon(data = MC.mouse.CS.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est, color = "95% CI"), 
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(data = MC.mouse.CS.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3, color = "50% CI"), 
              fill = "mediumblue", alpha = 0.3) +
  geom_line(data = MC.mouse.CS.plot, aes(x = Time, y = median_est, color = "Median"), 
            size = rel(1)) +
  geom_point(data = Obs.A1, aes(x = Time, y = CS, color = "Observed"), 
             shape = 1, fill = "white", size = 3, stroke = 2) +
  scale_color_manual(values = c("95% CI" = alpha("lightblue", 0.3), 
                                "50% CI" = alpha("mediumblue", 0.3), 
                                "Median" = "blue", 
                                "Observed" = "red"), 
                     labels = c("95% CI", "50% CI", "Median", "Observed")) + 
  ylab("Concentration in Spleen (ng/g)") + 
  xlab("Time (h)") +
  theme_minimal() +
  theme(axis.title = element_text(size = 22, color = "black"),  
        axis.text = element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(),    
        panel.grid.minor = element_blank(),     
        axis.line = element_line(color = "black")) +
  guides(color = guide_legend(title = "Legend"))  # Add legend title

#p.r.S

#ggsave(paste0(folder,"MCMC/MCMC_fit_spleen.png"),width = 10, height = 8)


library(gridExtra)
# Arrange the four plots together
combined_MCMC_plot <- grid.arrange(p.r.L,  p.r.lung, p.r.S,
                              ncol = 2, nrow = 2)
# Arrange the four plots together
combined_MCMC_plot <- grid.arrange(p.r.L, p.r.K, p.r.lung, p.r.S,
                              ncol = 2, nrow = 2)

# Save the combined plot
ggsave(paste0(folder, "MCMC/mcmc_fit_combined.png"), combined_MCMC_plot, width = 14, height = 10)





#-----------------------------



#---------------------3. trace plot and probability density function plot----------
## Trace plot using bayes plot
## Convergences plot
color_scheme_set("blue")

combinedchains        <- readRDS(file = paste0(folder,'MCMC/mouse.MCMC.comb.rds')) # use only first chain
R <-gelman.diag(combinedchains) # Gel man convergence diagnosis

mcmc_trace (
  combinedchains,
  pars =names(theta.MCMC[1:length(pars.mouse [-which_sig])]),
  size = 0.5,
  facet_args = list(nrow = 5)) +
  ggplot2::scale_color_brewer() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust angle of x-axis labels
ggsave(paste0(folder,"mc_sens/trace_plot.png"))

# todo probabilistic distribution of this difference with Fig3 
mcmc_dens_overlay(
  combinedchains,
  pars = names(theta.MCMC[1:length(pars.mouse [-which_sig])]),
  size=0.5,
  facet_args = list(nrow=5) + ggplot2::scale_color_brewer()
)
ggsave(paste0(folder,"mc_sens/prob_chains.png"))

write.csv(R[1],paste0(folder,"mc_sens/r.csv"))

