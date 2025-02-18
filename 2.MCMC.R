

# ---------------------------------------------------------
# Author: Jimeng Wu
# Email: jimeng.wu@empa.ch
# Date: 2025-02-18
# Description: MCMC analysis to generate the posterior parameter distribution from the prior parameter distribution
# ---------------------------------------------------------


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
source("dataset_info.R")

ls_np_name = c("Au: Study1_12nm_0.85mg/kg","Au: Study1_23nm_0.85mg/kg","Au: Study1_100nm_0.85mg/kg",
               "Au: Study2_34.6nm_3mg/kg","Au: Study2_55.5nm_3mg/kg","Au: Study2_77.1nm_3mg/kg",
               "Au: Study2_82.6nm_3mg/kg","Au: Study3_27.6nm_4.26mg/kg","Au: Study3_27.6nm_0.85mg/kg",
               "Si: Study1_20nm_10mg/kg","Si: Study1_80nm_10mg/kg","GO: Study1_20nm_20mg/kg",
               "GO: Study2_243nm_1mg/kg", "GO: Study2_914nm_1mg/kg_all","GO: Study2_914nm_1mg/kg_w/o_CS",
               "TiO2: Study1_385nm_10mg/kg","TiO2: Study2_220nm_60mg/kg",
               "FeO: Study1_29nm_5mg/kg","FeO: Study2_41nm_4mg/kg")


Obs.df = read_observation_data(np.name)$Obs.df
folder = read_observation_data(np.name)$folder
pathway = read_observation_data(np.name)$pathway
PDOSE = read_observation_data(np.name)$PDOSE

#---------------------Build mrgsolve-based PBPK Model-------
mod <- mcode ("mouse_PBPK", mousePBPK.code)

set.seed(5)

if (np.name == "FeO: Study2_41nm_4mg/kg") {
  tstep <- 0.5/60 #2.5/60 in the time step
} else if (np.name %in% c("GO: Study2_243nm_1mg/kg", "GO: Study2_914nm_1mg/kg_all",
                          "GO: Study2_914nm_1mg/kg_w/o_CS")) {
  # Assign the appropriate value for tstep
  tstep <- 1/60   # Fill in the appropriate value here 2/60 & 5/60 in the time points
} else {
  tstep <- min(1, min(Obs.df$Time))
}

if (pathway == "intraperitoneal injection") {
  pred.mouse <- pred.mouse.oral
} else if (pathway == "intravenous injection") {
  pred.mouse <- pred.mouse.iv
}

########################################## 4. Model Evalution with MCMC  ###########################################
#################################################################################################################

#------------------generate prior parameter distribution----------
#------------------Input parameters set----------

paras = read.csv(paste0(folder,'mod_fit/params_fitted.csv'))

params.mod = as.list(paras)$x

names(params.mod)<- names(params.init)

#-------------sig value for parameters and errors between obs and pred

# Define parameter names

sig_population = 0.5

# Create a named vector with default values
sig_values <- rep(sig_population, length(sig_names))
names(sig_values) = sig_names

# error with prediction and observation's deviation variance
sig_pred = 0.3
sig_values["sig2"] <- sig_pred

# Create sig_list
sig_list <- log(sig_values)

theta.MCMC<- c(params.mod, sig_list)


if (!file.exists(paste0(folder,"MCMC/"))) {
  dir.create(paste0(folder,"MCMC/"), recursive = TRUE)
}

#saveRDS(theta.MCMC,file =paste0(folder,'MCMC/theta.rds'))

##-------------------------prior distribution plotting--------------
which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG

# Mean and SD
mean.prior = exp(theta.MCMC[-which_sig])
sd.prior   = mean.prior*sig_population



#-------------- for generating the plot, draw samples from the settings---------

# Initialize lists to store samples for each parameter
param_samples <- list()
# Number of points to use in density estimation
n_points <- 10000  # Adjust as needed for smoother or rougher curves

CU             = 1                                              # Coefficient of uncertainty (CU) (Hack et al., 2006)
alpha          = (2+1)/(CU^2)                                   # Shape parameter  of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
beta           = (alpha-1)*sig_population^2                             # Scale parameter  of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F

for (param_index in seq_along(mean.prior)) {

  # Draw samples from the normal distribution for the current parameter
  param_samples[[param_index]] <- rnorm(n_points,
                                        mean = mean.prior[[param_index]], 
                                        sd =  sd.prior[[param_index]])
}



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
combined_prior_density_df <- do.call(rbind, density_data)



if (!file.exists(paste0(folder,"mc_sens/"))) {
  dir.create(paste0(folder,"mc_sens/"), recursive = TRUE)
}


# Calculate mean and standard deviation for each parameter
pri_mean_sd_df <- data.frame(Par = names(mean.prior),
                             Mean = unlist(mean.prior),
                             SD = unlist(sd.prior))


# Save to CSV
#write.csv(pri_mean_sd_df, file = paste0(folder,"mc_sens/pri_paras_stats.csv"), row.names = FALSE)


#saveRDS(combined_prior_density_df, file =paste0(folder,"mc_sens/pri_paras_density.rds"))




#--------------Maximum likelihood estimation (MLE) function for MCMC----
mcmc.fun <- function (pars){
  
  out <- pred.mouse(pars[-which_sig])

  out = out[which(out$Time %in% Obs.df$Time),]
  # making sure the predicted value and observed value matching with each other as
  # the time step in predicting function might not be small enough to cover the time step in the 
  # observed time step
  #Obs.df = Obs.df[which(Obs.df$Time %in% out$Time),] 
  col_order <- names(Obs.df)
  out = out[col_order]
  
  ## log-transformed prediction
  cols_to_log <- setdiff(names(out), "Time")
  out[cols_to_log] <- lapply(out[cols_to_log], function(x) ifelse(x != 0, log(x), NA))
  log.yhat <- unname(unlist(out[cols_to_log]))
  
  
  ## log-transformed experimental data
  Obs.df[cols_to_log] <- lapply(Obs.df[cols_to_log], function(x) ifelse(x != 0, log(x), NA))
  log.y <- unname(unlist(Obs.df[cols_to_log]))
  

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
sig_mean_population = 0.5
Prior <- function(pars) {
  
  ## Population level
  # The likelihood for population mean (parameters)
  # Calculate likelihoods of each parameters; P(u|mean,CV)
  pars.data = exp(pars[-which_sig]) # parameter value
  
  mean           = exp(theta.MCMC[-which_sig]) # parameter mean in the distribution
  CV             = sig_mean_population     # Coefficient of variation; Default value of 0.5 in all parameters (Bois,2000; Bois et al., 1996)
  
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
  prior_sig      = dinvgamma (sig, shape = alpha, scale = beta)  # prior distribution for population variance; sigma2
  
  
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
  MLE =  -2*sum(log.pri.pars, log.pri.sig , log.pri.pars.i,log.pri.sig2)  
  
  return(MLE)
}



#################### 5. MCMC simulation with parallel computing ############################
detectCores()                                ## check the cores
cl<- makeCluster(detectCores())              ## use all cores in our system     
registerDoParallel(cl)                       ## registers a cluster of all the cores on our system

# start time
tstr<-Sys.time()


niter = 600000
burninlength  = 300000
outputlength  = 300000


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


# Create the directory if it doesn't exist
if (!file.exists(paste0(folder,"MCMC/"))) {
  dir.create(paste0(folder,"MCMC/"), recursive = TRUE)
}

# output the MCMC results

#write.csv(MC.mouse.1,file=paste0(folder,"MCMC/mouse.chain1.csv"))
#saveRDS(MCMC,file =paste0(folder,'MCMC/mouse.MCMC.rds'))
#saveRDS(combinedchains,file=paste0(folder,'MCMC/mouse.MCMC.comb.rds'))


