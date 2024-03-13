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
folder = 'plots/100nm_short_v2/'
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_100nm_short.csv")  
set.seed(12)
#--------------------1. model mass balance checking---------------

## Build mrgsolve-based PBPK Model
mod_100nm <- mcode ("mouse_PBPK", mousePBPK_100nm.code)

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
tsamp=tgrid(0,tinterval*(TDoses-1)+24*7,tinterval)     ## Simulation time 24*7 hours (180 days)

## calculate the deposition volume
out <- 
  mod_100nm %>% 
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


nm_size = "100nm"

plot (y=out$CL, x=out$Time,xlab = "Time", main = paste("Liver concentration",nm_size),
      type='l',log="x")
plot (y=out$CLt, x=out$Time,xlab = "Time", main = paste("Liver Tissue concentration",nm_size),
      type='l',ylim=c(0,10000))
plot (y=out$CSt, x=out$Time,xlab = "Time", main = paste("Spleen Tissue concentration",nm_size),
      type='l',ylim=c(0,10000))
plot (y=out$CKt, x=out$Time,xlab = "Time", main = paste("Kidney Tissue concentration",nm_size),
      type='l',ylim=c(0,800))
plot (y=out$CLungt, x=out$Time,xlab = "Time", main = paste("Lung Tissue concentration",nm_size),
      type='l',ylim=c(0,1500))


plot.l=
  ggplot() +
  geom_line(data  = out, aes(Time,CLt), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CL), size=2.5) + 
  ylab("Concentration in Liver (ng/g)") +
  xlab("Time (h)")
plot.l
ggsave(paste(folder,"mb_liver.png"))


plot.s=
  ggplot() +
  geom_line(data  = out, aes(Time,CSt), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CS), size=2.5) + 
  ylab("Concentration in Spleen (ng/g)") +
  xlab("Time (h)")
plot.s
ggsave(paste(folder,"mb_spleen.png"))

plot.k=
  ggplot() +
  geom_line(data  = out, aes(Time,CKt), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CK), size=2.5) + 
  ylab("Concentration in Kidney (ng/g)") +
  xlab("Time (h)")
plot.k
ggsave(paste(folder,"mb_kidney.png"))


plot.lu=
  ggplot() +
  geom_line(data  = out, aes(Time,CLungt), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, Clung), size=2.5) + 
  ylab("Concentration in Lung (ng/g)") +
  xlab("Time (h)")
plot.lu
ggsave(paste(folder,"mb_lung.png"))


#--------------------2. initial parameter sensitivity analysis------------
## A1 data set = iv single dose of 0.85 mg/kg; Matrix: plasma; 13nm. 
## https://linkinghub.elsevier.com/retrieve/pii/S0041-008X(10)00072-4
## A2 data set same as above, but with 100nm

pred.mouse <- function(pars) {

    ## Get out of log domain
    pars %<>% lapply(exp)

    ## Define the exposure scenario

    BW           = 0.02                              ## kg, body weight
    tinterval    = 0.1                                 ## hr, Time interval
    TDoses       = 1                                 ## Dose times, only one dose
    PDOSE        = 0.85                              ## mg/kg-day, Single dose
    DOSE         = PDOSE*BW                          ## mg, amount of iv dose
    ex.iv<- ev(ID=1, amt= DOSE,                  ## Set up the exposure events
            ii=tinterval, addl=TDoses-1, 
            cmt="MBV", replicate = FALSE) 

    ## Set up the exposure time
    tsamp=tgrid(0,tinterval*(TDoses-1)+24*7,tinterval)     ## Simulation time 24*7 hours (180 days)

    ## calculate the deposition volume
    out <- 
    mod_100nm %>% 
    param(pars) %>%
    ##Req(Liver,M_tot,MBV)%>%
    update(atol = 1E-80,maxsteps = 5000000) %>%
    mrgsim_d(data = ex.iv, tgrid=tsamp)
    
    ## save the calculated into data frame
    out <- data.frame(Time=out$time, 
                      CL=out$Liver_t,
                      CK = out$Kidney_t,
                      CS = out$Spleen_t,
                      Clung = out$Lung_t)

    return(out)
}

## 2.1 initial parameters-----------

# for 100nm 
params.init <- log(c(
  K_release_Liver = 0.0075,  # h-1
  K_max_Liver = 4,         # h-1
  K_50_Liver = 24,          # h
  n_Liver = 0.1,              # Unitless
  #K_release_GI = 0.003,     # h-1
  #K_max_GI = 0.075,         # h-1
  #K_50_GI = 24,             # h
  #n_GI = 5,                 # Unitless
  K_release_Spleen = 0.003, # h-1
  K_max_Spleen = 10,        # h-1
  K_50_Spleen = 24,
  n_Spleen = 0.1,
  K_release_Kidney = 0.01, # h-1
  K_max_Kidney = 0.1,
  K_50_Kidney = 24,
  n_Kidney = 0.1,
  K_release_Lung = 0.005,   # h-1
  K_max_Lung = 0.1,
  K_50_Lung = 24,
  n_Lung = 0.1,               
  P_Liver  = 0.08,
  P_Brain  = 0.15,
  P_Kidney  = 0.15,
  P_Spleen  = 0.15,
  P_Lung  = 0.15,
  P_Rest  = 0.15,
  DLC_Liver = 0.001,
  DLC_Brain = 0.000001,
  DLC_Kidney = 0.001,
  DLC_Spleen = 0.001,
  DLC_Lung = 0.001,
  #DLC_GI = 0.001,
  DLC_Rest = 0.000001,
  Kbile = 0.0012,       # Biliary clearance (L/hr)
  Kurine = 0.00012     # Urine clearance (L/hr)
))


### Local sensitivity analysis
# Choose the sensitive parameters in the model

## Sensitivity function (FME)
## Check the sensitive parameters in the model ,sensvar = c("CL","CK","CS","CB","Clung")
Sens <- sensFun(func = pred.mouse, parms = params.init,varscale = 1)

head(Sens)

png(paste(folder,"pars_sens.png"),width=2500,height=2800,res=300)
plot(Sens)
dev.off()

plot(Sens, type = "b", pch = 15:19, col = 2:6, 
     main = "Sensitivity all vars")

require(reshape2)
# for making separated plot
Sens_1 <- melt(Sens , id.vars = c('x',"var"), variable.name = 'series')
ggplot(subset(Sens_1,series==c("K_release_Liver","K_max_Liver","K_50_Liver",
                               "n_Liver","K_release_Spleen")), 
       aes(x,value)) + geom_line(aes(colour = series))

df_Sens=summary(Sens)
write.csv(df_Sens,file=paste(folder,"pars_sens.csv"))
png(paste(folder,"pars_sens_summary.png"),width=2500,height=3000,res=300)
plot(df_Sens)
dev.off()

png(paste(folder,"pars_corr.png"),width=2500,height=1800,res=300)
pairs(Sens)
dev.off()

# ALMOST EVERY PARAMETERS ARE CORRELATED WITH EACH OTHER




## --------2.2 selected sensitive parameters-----------------------
### set up sensitive or necessary parameters as model input
# for 100nm 
params2fit <- log(c(
  K_release_Spleen = 0.003, # h-1
  K_50_Spleen = 24,
  K_release_Lung = 0.005,   # h-1
  DLC_Kidney = 0.001,
  DLC_Lung = 0.001
))

## Cost function (FME)
## Estimate the model residual by modCost function
MCcost<-function (pars, obs){
  out<- pred.mouse(pars)
  cost<- modCost(model=out,obs=obs,weight='std',x="Time")
  return(cost)
}





#------------------3. Fitting with A1 dataset using modFit function-----------------------------
Obs.A1 <- Obs.A1[ ,!(colnames(Obs.A1) == "CK")]
Fit.Result.A1<- modFit(f=MCcost, p=params2fit, obs=Obs.A1, method ="bobyqa", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"

summary(Fit.Result.A1)                           ## Summary of fit
#exp(Fit.Result$par)                          ## Get the arithmetic value out of the log domain

res.A1=MCcost(Fit.Result.A1$par, obs=Obs.A1)$residuals$res      ## Check the residual for each time points
sum(res.A1^2)                                    ## Total residuals 

# Calculated the model output with fitted parameters
Fitted_output.A1 = pred.mouse(par=Fit.Result.A1$par)


plot.A1=
  ggplot() +
  geom_line(data  = Fitted_output.A1, aes(Time,CL), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CL), size=2.5) + 
  ylab("Concentration in Liver (ng/g)") +
  xlab("Time (h)")
plot.A1
ggsave(paste(folder,"init_fit_liver.png"))


plot.A1_Lung=
  ggplot() +
  geom_line(data  = Fitted_output.A1, aes(Time,Clung), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, Clung), size=2.5) + 
  ylab("Concentration in Lung (ng/g) ") +
  xlab("Time (h)")
plot.A1_Lung
ggsave(paste(folder,"init_fit_lung.png"))


#plot.A1_Kidney=
#  ggplot() +
#  geom_line(data  = Fitted_output.A1, aes(Time,CK), col="firebrick", lwd=2)+
#  geom_point(data = Obs.A1    , aes(Time, CK), size=2.5) + 
#  ylab("Concentration in Kidney (ng/g)") +
#  xlab("Time (h)")
#plot.A1_Kidney
#ggsave(paste(folder,"init_fit_kidney.png"))

plot.A1_Spleen=
  ggplot() +
  geom_line(data  = Fitted_output.A1, aes(Time,CS), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CS), size=2.5) + 
  ylab("Concentration in Spleen (ng/g)") +
  xlab("Time (h)")
plot.A1_Spleen
ggsave(paste(folder,"init_fit_spleen.png"))

#plot.A2_GI=
#  ggplot() +
#  geom_line(data  = Fitted_output.A2, aes(Time,CG), col="firebrick", lwd=2)
#plot.A2_GI





########################################## 4. Model Evalution with MCMC  ###########################################
## Oral.obs.A-C:  Cynomolgus Monkeys oral daily dose to 0.03,0.15,0.75 mg/kg for 182 days and monitored 1 year  #
## B1: 0.03 mg/kg-d                                                                                             #
## Matrix: plasma, data from Seacat et al. (2002)                                                               #                                                 
#################################################################################################################


#------Input parameters set----------

theta.MCMC<-log(c(
  
  K_release_Liver = 0.0075,  # h-1
  K_max_Liver = 4,         # h-1
  K_50_Liver = 24,          # h
  n_Liver = 0.1,              # Unitless
  K_release_Spleen = as.double(exp(Fit.Result.A1$par[1])), # h-1
  K_max_Spleen = 10,        # h-1
  K_50_Spleen = as.double(exp(Fit.Result.A1$par[2])),
  n_Spleen = 0.1,
  K_release_Kidney = 0.01, # h-1
  K_max_Kidney = 0.1,
  K_50_Kidney = 24,
  n_Kidney = 0.1,
  K_release_Lung = as.double(exp(Fit.Result.A1$par[3])),   # h-1
  K_max_Lung = 0.1,
  K_50_Lung = 24,
  n_Lung = 0.1,               
  P_Liver  = 0.08,
  P_Brain  = 0.15,
  P_Kidney  = 0.15,
  P_Spleen  = 0.15,
  P_Lung  = 0.15,
  P_Rest  = 0.15,
  DLC_Liver = 0.001,
  DLC_Brain = 0.000001,
  DLC_Kidney = as.double(exp(Fit.Result.A1$par[4])),
  DLC_Spleen = 0.001,
  DLC_Lung = as.double(exp(Fit.Result.A1$par[5])),
  #DLC_GI = 0.001,
  DLC_Rest = 0.000001,
  Kbile = 0.0012,       # Biliary clearance (L/hr)
  Kurine = 0.00012,     # Urine clearance (L/hr)
  sig2                  = 0.5, ## Model error (residuals); 
  ## mostly between 0.3 and 0.5 (corresponding to 
  ## coefficients of variation of about 30-50%); Bois et al. (1998, 2000)
  
  ## population variance; equal to the CV of parameters (this study assued the cv of all parametesr is 0.3)
  sig_K_release_Liver = 0.3,  # h-1
  sig_K_max_Liver = 0.3,         # h-1
  sig_K_50_Liver = 0.3,          # h
  sig_n_Liver = 0.3,              # Unitless
  sig_K_release_Spleen = 0.3, # h-1
  sig_K_max_Spleen = 0.3,        # h-1
  sig_K_50_Spleen = 0.3,
  sig_n_Spleen = 0.3,
  sig_K_release_Kidney = 0.3, # h-1
  sig_K_max_Kidney = 0.3,
  sig_K_50_Kidney = 0.3,
  sig_n_Kidney = 0.3,
  sig_K_release_Lung = 0.3,   # h-1
  sig_K_max_Lung = 0.3,
  sig_K_50_Lung = 0.3,
  sig_n_Lung = 0.3,               
  sig_P_Liver  = 0.3,
  sig_P_Brain  = 0.3,
  sig_P_Kidney  = 0.3,
  sig_P_Spleen  = 0.3,
  sig_P_Lung  = 0.3,
  sig_P_Rest  = 0.3,
  sig_DLC_Liver = 0.3,
  sig_DLC_Brain = 0.3,
  sig_DLC_Kidney = 0.3,
  sig_DLC_Spleen = 0.3,
  sig_DLC_Lung = 0.3,
  sig_DLC_Rest = 0.3,
  sig_Kbile = 0.3,       # Biliary clearance (L/hr)
  sig_Kurine = 0.3     # Urine clearance (L/hr)
))


saveRDS(theta.MCMC,file =paste(folder,'theta.rds'))

which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG


#--------------Maximum likelihood estimation (MLE) function for MCMC----
mcmc.fun <- function (pars, pred=FALSE){
  
  ## Get out of log domain
  pars.data <- lapply(pars [-which_sig],exp)
  
  
  # Exposure scenario for single oral dose of 0.85 mg/kg
  
  
  BW           = 0.02                              ## kg, body weight
  tinterval    = 0.1                                 ## hr, Time interval
  TDoses       = 1                                 ## Dose times, only one dose
  PDOSE        = 0.85                              ## mg/kg-day, Single dose
  DOSE         = PDOSE*BW                          ## mg, amount of iv dose
  ex.iv<- ev(ID=1, amt= DOSE,                  ## Set up the exposure events
             ii=tinterval, addl=TDoses-1, 
             cmt="MBV", replicate = FALSE) 
  
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*7,tinterval)     ## Simulation time 24*7 hours (180 days)
  
  
  ## Get a prediction
  out <- 
    mod_100nm %>%
    param(pars.data) %>%
    update(atol = 1E-8, maxsteps= 5000) %>%
    mrgsim_d(data = ex.iv, tgrid=tsamp)
  
  out <- data.frame(Time=out$time, 
                    CL=out$Liver_t,
                    #CK = out$Kidney_t,
                    CS = out$Spleen_t,
                    Clung = out$Lung_t)
  
  if (pred) return (out)
  
  out = out[which(out$Time %in% Obs.A1$Time),]
  
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
  log.yhat        <- c(log.yhat.CL,log.yhat.CS,log.yhat.Clung)
  log.y           <- c(log.y.CL,log.y.CS,log.y.Clung)
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
  pars.data = exp(pars[-which_sig])
  sig  <- as.numeric (exp(pars[which_sig][2:31]))                 # Coefficient of variation from population variance; sigmal0
  sig2 <- as.numeric (exp(pars[which_sig][1]))                    # error variances from model residual
  
  mean           = exp(theta.MCMC[-which_sig])
  CV             = 0.5                                            # Coefficient of variation; Default value of 0.5 in all parameters (Bois,2000; Bois et al., 1996)
  sd             = mean*CV
  
  # Calculate likelihoods of each parameters; P(u|M,S)
  prior_pars     = dtruncnorm(pars.data, 
                              a = qnorm(0.025, mean = mean, sd = sd), 
                              b = qnorm(0.975, mean = mean, sd = sd), 
                              mean = mean, sd = sd ) 
  
  # The likelihood for population variance; P(sigmal^2|sigmal0^2)
  CU             = 1                                              # Coefficient of uncertainty (CU) (Hack et al., 2006)
  CV.sig         = exp(theta.MCMC[which_sig])[2:31]               # Sigmal0
  alpha          = (2+1)/(CU^2)                                   # Shape parameter  of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  beta           = (alpha-1)*CV.sig^2                             # Scale parameter  of gamma distribution; Appendix Table A-7 from EPA (2011) :EPA/635/R-09/011F
  
  # Calculate likelihoods of model error (sig2) and population variance (sig) parameters
  prior_sig      = dinvgamma (sig, shape = alpha , scale = beta)  # prior distribution for population variance; sigma2
  prior_sig2     = dunif (sig2, min = 0.01, max = 3.3)            # error variances, Lower and upper boundary from Chiu et al., 2009; Chiu et al., 2014)   
  
  ## individual level; P(theta|u,sigmal^2)
  mean_i         = prior_pars
  sd_i           = sqrt(prior_sig)
  prior_pars_i   = dtruncnorm (prior_pars, 
                               a = qnorm(0.025, mean = mean_i, sd = sd_i), 
                               b = qnorm(0.975, mean = mean_i, sd = sd_i), 
                               mean = mean_i, sd = sd_i) 
  
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

# parallel
system.time(
  MCMC <- foreach( i = 1:4, .packages = c('mrgsolve','magrittr','FME',
                                          'truncnorm','EnvStats',
                                          'invgamma','dplyr')) %dopar% {
                                            mod_100nm <- mcode ("mouse_PBPK", mousePBPK_100nm.code)
                                            modMCMC(f             = mcmc.fun, 
                                                    p             = theta.MCMC, 
                                                    niter         = 500000,           ## iteration number 
                                                    jump          = 0.01,             ## jump function generation new parameters distribution using covrate matrix
                                                    prior         = Prior,            ## prior function
                                                    updatecov     = 50,               ## adaptive Metropolis
                                                    var0          = NULL,             ## initial model variance;
                                                    wvar0         = 0.01,             ## "weight" for the initial model variance
                                                    ntrydr        = 2,                ## delayed Rejection
                                                    burninlength  = 250000,           ## number of initial iterations to be removed from output.
                                                    outputlength  = 50000,
                                                    verbose=5)            ## number of output iterations           
                                            
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

## Save the posterior parameters (95% CI)
quan.mouse = exp(summary(MC.mouse.1)$quantiles)  


# output the MCMC results
write.csv(quan.mouse,file=paste(folder,"mouse.summary_pos.csv"))
write.csv(MC.mouse.1,file=paste(folder,"mouse.pos.csv"))
saveRDS(MCMC,file =paste(folder,'mouse.MCMC.rds'))
saveRDS(combinedchains,file=paste(folder,'mouse.comb.rds'))




## Plot using MCMC posterior parameters
## Model validation using one of the best posterior parameter combination
Sim.fit.MCMC.B = mcmc.fun (par = MCMC[[1]]$bestpar,pred=TRUE)

df.sim.MCMC.B = cbind.data.frame (Time=Sim.fit.MCMC.B$Time, 
                                  CL=Sim.fit.MCMC.B$CL,
                                  CK=Sim.fit.MCMC.B$CK)

plot.BL =
  ggplot() +
  geom_line(data = df.sim.MCMC.B,aes(Time,CL), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1 ,aes(Time, CL),size=2.5) + ylab("Concentration") 
plot.BL
plot.BK =
  ggplot() +
  geom_line(data = df.sim.MCMC.B,aes(Time,CK), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1 ,aes(Time, CK),size=2.5) + ylab("Concentration") 
plot.BK

#--------------------------------------plots-------------------------
#--------------------1. Densities of posterior parameter uncertainty distributions of the population mean (μ). ----------

Mouse.MCMC        <- readRDS(file = paste(folder,"mouse.MCMC.rds"))[[1]] # use only first chain


## loading the theta names
theta.names       <- names(readRDS(file = paste(folder,"theta.rds")))
## Sampling from posterior parameters to generate the posterior distributions 
## Mouse posteiror distributions
M.Mouse  <- exp(Mouse.MCMC$pars) %>% apply(2,mean)
SD.Mouse <- exp(Mouse.MCMC$pars) %>% apply(2,sd)
dis.Mouse <- matrix(ncol = 50000, nrow = 30)
rownames(dis.Mouse) <- names(M.Mouse)[1:30]

for (i in 1:50000){
  for (j in 1:30){
    dis.Mouse[j,i] <- rnorm (1, mean = M.Mouse[j], sd= SD.Mouse[j])
  }
}

dis.Mouse<-melt(dis.Mouse)
names(dis.Mouse)=c("Par","Species","Value")
dis.Mouse$Species = c("Mouse")

## Summary the mouse, rat, monkey, and human posterior distribution in the "Pos.mean" object
Pos.mean <- rbind.data.frame (dis.Mouse)
Pos.mean$log.value <- log(Pos.mean$Value)

## Fig. 3: Densities of posterior parameter uncertainty distributions
## P1: All species gather together in one plot
Pos.mean$ggjoy = c("A") # Create a new column ggjoy in Pos.mean table. Without this, 
# the plot p1 will have four layers for four species overlapped in each panel, which does not look good.

# With this, y axis is A for all species in one layer in one panel.

Pos.mean$Par = factor(Pos.mean$Par, levels = theta.names[1:30])
library(ggjoy)     # Used to create Figure 3
p1 = 
  ggplot (Pos.mean, aes(x = as.numeric(log.value), y = as.factor(ggjoy),fill = Species)) + # fill, fill one color for each species
  geom_joy (scale = 8, size = 0.25, rel_min_height = 0.01, alpha = 0.4) + # over size of the entire ggjoy plot
  #geom_joy (data=filter(Pos.mean, Species == "Monkey"),scale = 28, size = 0.25, rel_min_height = 0.01, fill="red", alpha=0.3) +
  #geom_joy (data=filter(Pos.mean, Species == "Mouse"), scale = 8, size = 0.25, rel_min_height = 0.01, fill="green", alpha=0.3) +
  #geom_joy (data=filter(Pos.mean, Species == "Rat"),   scale = 28, size = 0.25, rel_min_height = 0.01, fill="pink", alpha=0.3)+
  scale_fill_manual(values = c("pink", "grey", "blue","red")) + 
  scale_colour_manual(values = c("pink", "grey", "blue","red"))+ 
  facet_wrap ( ~ Par, nrow = 6, ncol= 5,scale="free_x") # arrange the plot by parameter name, free_x means that the x scale is not fixed, so automatically updated.

p1=
  p1 + theme_joy()+ # change the theme of the plot p1 (incorporate p1 into the theme)
  theme (
    plot.background         = element_blank(),
    text                    = element_text (family = "Times",face="bold",size = 18),
    #panel.border            = element_rect (colour = "black", fill=NA, size=2),
    #panel.background        = element_blank(),
    panel.background        = element_rect (fill   = "#f7f7f7"),
    #panel.grid.major.y      = element_line (size   = 0.5, colour = "grey"),
    #panel.grid.minor.y      = element_blank(),
    panel.grid.major.x      = element_blank(), 
    #panel.grid.minor.x      = element_blank(), 
    #axis.text               = element_text (size   = 18, colour = "black",face="bold"),  
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    axis.ticks.x            = element_line (size   = 1.5, colour = "black"), 
    axis.ticks.y            = element_blank(), 
    axis.text.y             = element_blank(),
    strip.text              = element_text (size   = 18),
    legend.title            = element_text (size   = 18, face="bold"),
    legend.justification    = "none",
    legend.position         = "none",
    legend.text             = element_text (size = 18,face="bold")) + 
  labs (x = "", y = "") 

# Save the image by using the Export button, adjust the height and width, choose image format, preview
ggsave("Fig.3.tiff",scale = 1,
       plot = p1,
       path = folder,
       width = 25, height = 20, units = "cm",dpi=320)


dev.off()



#-----------------------------2. get sensitivity range plots-------------
Mouse.MCMC        <- readRDS(file = paste(folder,"mouse.MCMC.rds")) # use only first chain

Newtime.r   = pred.mouse(theta.MCMC)$Time  
# this is the new time variable, now it has been changed to sample per day.
nrwo.r = length (Newtime.r)

MC.mouse.a.CL    = matrix(nrow = nrwo.r, ncol = 5000)
MC.mouse.a.CK    = matrix(nrow = nrwo.r, ncol = 5000)
MC.mouse.a.Clung    = matrix(nrow = nrwo.r, ncol = 5000)
MC.mouse.a.CS    = matrix(nrow = nrwo.r, ncol = 5000)

for(i in 1:5000){
  
  j = i *10  
  pars.mouse             = MCMC[[1]]$pars    [j,]     # sample parameter set once every ten sets, so you will have 5000 sets from 50000 total sets
  
  MCdata               = pred.mouse (pars.mouse)
  MC.mouse.a.CL[,i]= MCdata$CL
  MC.mouse.a.CK[,i]= MCdata$CK
  MC.mouse.a.Clung[,i]= MCdata$Clung
  MC.mouse.a.CS[,i]= MCdata$CS
  #cat("iteration = ", i , "\n") # Shows the progress of iterations, so you can see the number of iterations that has been completed, and how many left.
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
# for plotting the kidney results
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_100nm_short.csv")  

p.r.L <- 
  ggplot() + 
  geom_ribbon(data = MC.mouse.CL.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="yellowgreen", alpha=0.3) +
  geom_ribbon(data = MC.mouse.CL.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="green4", alpha = 0.3) +
  geom_line(data= MC.mouse.CL.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "green3") +

  geom_point(data=Obs.A1, aes(x=Time, y= CL), shape = 1, colour = "black", 
             fill = "white", size = 1, stroke = 2) +
  ylab("Concentration in Liver (ng/g)") + xlab("Time (h)")
p.r.L
ggsave(paste(folder,"MCMC_fit_liver.png"))

p.r.L

p.r.K <- 
  ggplot() + 
  geom_ribbon(data = MC.mouse.CK.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="yellowgreen", alpha=0.3) +
  geom_ribbon(data = MC.mouse.CK.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="green4", alpha = 0.3) +
  geom_line(data= MC.mouse.CK.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "green3") +
  geom_point(data=Obs.A1, aes(x=Time, y= CK), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  ylab("Concentration in Kidney (ng/g)") + xlab("Time (h)")
p.r.K
ggsave(paste(folder,"MCMC_fit_Kidney.png"))


p.r.lung <- 
  ggplot() + 
  geom_ribbon(data = MC.mouse.Clung.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="yellowgreen", alpha=0.3) +
  geom_ribbon(data = MC.mouse.Clung.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="green4", alpha = 0.3) +
  geom_line(data= MC.mouse.Clung.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "lightgreen") +
  geom_point(data=Obs.A1, aes(x=Time, y= Clung), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  ylab("Concentration in Lung (ng/g)") + xlab("Time (h)")
p.r.lung
ggsave(paste(folder,"MCMC_fit_lung.png"))


p.r.S <- 
  ggplot() + 
  geom_ribbon(data = MC.mouse.CS.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), 
              fill="yellowgreen", alpha=0.3) +
  geom_ribbon(data = MC.mouse.CS.plot, aes(x = Time, ymin = ci_q1, ymax = ci_q3), 
              fill="green4", alpha = 0.3) +
  geom_line(data= MC.mouse.CS.plot, aes(x = Time, y = median_est), 
            size = rel(1), colour = "lightgreen") +
  geom_point(data=Obs.A1, aes(x=Time, y= CS), shape = 1, colour = "black", 
             fill = "white", size = 3, stroke = 2) +
  ylab("Concentration in Spleen (ng/g)") + xlab("Time (h)")
p.r.S
ggsave(paste(folder,"MCMC_fit_spleen.png"))


#---------------------3. trace plot and probability density function plot----------
## Trace plot using bayes plot
## Convergences plot
color_scheme_set("blue")

combinedchains        <- readRDS(file = paste(folder,"mouse.comb.rds")) # use only first chain
R <-gelman.diag(combinedchains) # Gel man convergence diagnosis

mcmc_trace (
  combinedchains,
  pars =names(theta.MCMC[1:30]),
  size = 0.5,
  facet_args = list(nrow = 5)) +
  ggplot2::scale_color_brewer()
ggsave(paste(folder,"trace_plot.png"))


mcmc_dens_overlay(
  combinedchains,
  pars = names(theta.MCMC[1:30]),
  size=0.5,
  facet_args = list(nrow=5) + ggplot2::scale_color_brewer()
)
ggsave(paste(folder,"prob_chains.png"))

write.csv(R[1],paste(folder,"r.csv"))


#----------------------------------------old---------------------------------------------
#-------------------------------------Fitting with A1 dataset---------------
Obs.A1 <- read.csv(file ="C:/switchdriver/dataset/tk/mouse/R_input_mouse_study1_4nm.csv")



Fit.Result.A1<- modFit(f=MCcost, p=params2fit, obs=Obs.A1, method ="Nelder-Mead", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"

summary(Fit.Result.A1)                           ## Summary of fit
#exp(Fit.Result$par)                          ## Get the arithmetic value out of the log domain
res.A1=MCcost(Fit.Result.A1$par, obs=Obs.A1)$residuals$res      ## Check the residual for each time points
sum(res.A1^2)                                    ## Total residuals 

# Calculated the model output with fitted parameters

Fitted_output.A1 = pred.mouse(par=Fit.Result.A1$par)
                                
plot.A1=
    ggplot() +
    geom_line(data  = Fitted_output.A1, aes(Time,CL), col="firebrick", lwd=2)+
    geom_point(data = Obs.A1    , aes(Time, CL), size=2.5) + ylab("Concentration") 
plot.A1

plot.A1_Lung=
  ggplot() +
  geom_line(data  = Fitted_output.A1, aes(Time,Clung), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, Clung), size=2.5) + ylab("Concentration") 
plot.A1_Lung

plot.A1_Kidney=
  ggplot() +
  geom_line(data  = Fitted_output.A1, aes(Time,CK), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CK), size=2.5) + ylab("Concentration") 
plot.A1_Kidney

plot.A1_Spleen=
  ggplot() +
  geom_line(data  = Fitted_output.A1, aes(Time,CS), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CS), size=2.5) + ylab("Concentration") 
plot.A1_Spleen


plot.A1_Brain=
  ggplot() +
  geom_line(data  = Fitted_output.A1, aes(Time,CB), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CB), size=2.5) + ylab("Concentration") 
plot.A1_Brain

plot.A1_GI=
  ggplot() +
  geom_line(data  = Fitted_output.A1, aes(Time,CG), col="firebrick", lwd=2)
plot.A1_GI

#-------------------fitting with only liver data-------------------
Obs.A1_L <- data.frame(Time=Obs.A1['Time'], 
                       CL=Obs.A1['CL'])

Fit.Result.A1_L<- modFit(f=MCcost, p=params2fit, obs=Obs.A1_L, method ="Nelder-Mead", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"

summary(Fit.Result.A1_L)                           ## Summary of fit
#exp(Fit.Result$par)                          ## Get the arithmetic value out of the log domain
res.A1_L=MCcost(Fit.Result.A1_L$par, obs=Obs.A1)$residuals$res      ## Check the residual for each time points
sum(res.A1_L^2)                                    ## Total residuals 

# Calculated the model output with fitted parameters

Fitted_output.A1_L = pred.mouse(par=Fit.Result.A1_L$par)

plot_L.A1=
  ggplot() +
  geom_line(data  = Fitted_output.A1_L, aes(Time,CL), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CL), size=2.5) + ylab("Concentration") 
plot_L.A1

plot_L.A1_Lung=
  ggplot() +
  geom_line(data  = Fitted_output.A1_L, aes(Time,Clung), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, Clung), size=2.5) + ylab("Concentration") 
plot.A1_Lung

plot_L.A1_Kidney=
  ggplot() +
  geom_line(data  = Fitted_output.A1_L, aes(Time,CK), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CK), size=2.5) + ylab("Concentration") 
plot_L.A1_Kidney

plot_L.A1_Spleen=
  ggplot() +
  geom_line(data  = Fitted_output.A1_L, aes(Time,CS), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CS), size=2.5) + ylab("Concentration") 
plot_L.A1_Spleen


plot_L.A1_Brain=
  ggplot() +
  geom_line(data  = Fitted_output.A1_L, aes(Time,CB), col="firebrick", lwd=2)+
  geom_point(data = Obs.A1    , aes(Time, CB), size=2.5) + ylab("Concentration") 
plot.A1_Brain

plot_L.A1_GI=
  ggplot() +
  geom_line(data  = Fitted_output.A1_L, aes(Time,CG), col="firebrick", lwd=2)
plot_L.A1_GI





#-------------------fitting with only A2 liver data-------------------
Obs.A2_L <- data.frame(Time=Obs.A2['Time'], 
                       CL=Obs.A2['CL'])

Fit.Result.A2_L<- modFit(f=MCcost, p=params2fit, obs=Obs.A2_L, method ="Nelder-Mead", 
                         control = nls.lm.control(nprint=1)) #"Nelder-Mead"

summary(Fit.Result.A2_L)                           ## Summary of fit
#exp(Fit.Result$par)                          ## Get the arithmetic value out of the log domain
res.A2_L=MCcost(Fit.Result.A2_L$par, obs=Obs.A2)$residuals$res      ## Check the residual for each time points
sum(res.A2_L^2)                                    ## Total residuals 

# Calculated the model output with fitted parameters

Fitted_output.A2_L = pred.mouse(par=Fit.Result.A2_L$par)

plot_L.A2=
  ggplot() +
  geom_line(data  = Fitted_output.A2_L, aes(Time,CL), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CL), size=2.5) + ylab("Concentration") 
plot_L.A2

plot_L.A2_Lung=
  ggplot() +
  geom_line(data  = Fitted_output.A2_L, aes(Time,Clung), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, Clung), size=2.5) + ylab("Concentration") 
plot_L.A2_Lung

plot_L.A2_Kidney=
  ggplot() +
  geom_line(data  = Fitted_output.A2_L, aes(Time,CK), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CK), size=2.5) + ylab("Concentration") 
plot_L.A2_Kidney

plot_L.A2_Spleen=
  ggplot() +
  geom_line(data  = Fitted_output.A2_L, aes(Time,CS), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CS), size=2.5) + ylab("Concentration") 
plot_L.A2_Spleen


plot_L.A2_Brain=
  ggplot() +
  geom_line(data  = Fitted_output.A2_L, aes(Time,CB), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CB), size=2.5) + ylab("Concentration") 
plot.A2_Brain

plot_L.A2_GI=
  ggplot() +
  geom_line(data  = Fitted_output.A2_L, aes(Time,CG), col="firebrick", lwd=2)
plot_L.A2_GI





#-------------------fitting with only spleen data-------------------
Obs.A2_S <- data.frame(Time=Obs.A2['Time'], 
                       CL=Obs.A2['CS'])

Fit.Result.A2_S<- modFit(f=MCcost, p=params2fit, obs=Obs.A2_S, method ="Nelder-Mead", 
                         control = nls.lm.control(nprint=1)) #"Nelder-Mead"

summary(Fit.Result.A2_S)                           ## Summary of fit
#exp(Fit.Result$par)                          ## Get the arithmetic value out of the log domain
res.A2_S=MCcost(Fit.Result.A2_S$par, obs=Obs.A2_S)$residuals$res      ## Check the residual for each time points
sum(res.A2_S^2)                                    ## Total residuals 

# Calculated the model output with fitted parameters

Fitted_output.A2_S = pred.mouse(par=Fit.Result.A2_S$par)

plot_S.A2=
  ggplot() +
  geom_line(data  = Fitted_output.A2_S, aes(Time,CL), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CL), size=2.5) + ylab("Concentration") 
plot_S.A2

plot_S.A2_Lung=
  ggplot() +
  geom_line(data  = Fitted_output.A2_S, aes(Time,Clung), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, Clung), size=2.5) + ylab("Concentration") 
plot.A2_Lung

plot_S.A2_Kidney=
  ggplot() +
  geom_line(data  = Fitted_output.A2_S, aes(Time,CK), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CK), size=2.5) + ylab("Concentration") 
plot_S.A2_Kidney

plot_S.A2_Spleen=
  ggplot() +
  geom_line(data  = Fitted_output.A2_S, aes(Time,CS), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CS), size=2.5) + ylab("Concentration") 
plot_S.A2_Spleen


plot_S.A2_Brain=
  ggplot() +
  geom_line(data  = Fitted_output.A2_S, aes(Time,CB), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CB), size=2.5) + ylab("Concentration") 
plot.A2_Brain

plot_S.A2_GI=
  ggplot() +
  geom_line(data  = Fitted_output.A2_S, aes(Time,CG), col="firebrick", lwd=2)
plot_S.A2_GI







#-------------------fitting with only brain data-------------------
# only brain has the different shape with experiment data, but the results of other organ 
# becomes worse when we change to use only brain data
Obs.A2_B <- data.frame(Time=Obs.A2['Time'], 
                       CL=Obs.A2['CB'])

Fit.Result.A2_B<- modFit(f=MCcost, p=params2fit, obs=Obs.A2_B, method ="Nelder-Mead", 
                         control = nls.lm.control(nprint=1)) #"Nelder-Mead"

summary(Fit.Result.A2_B)                           ## Summary of fit
#exp(Fit.Result$par)                          ## Get the arithmetic value out of the log domain
res.A2_B=MCcost(Fit.Result.A2_B$par, obs=Obs.A2_B)$residuals$res      ## Check the residual for each time points
sum(res.A2_B^2)                                    ## Total residuals 

# Calculated the model output with fitted parameters

Fitted_output.A2_B = pred.mouse(par=Fit.Result.A2_B$par)

plot_B.A2=
  ggplot() +
  geom_line(data  = Fitted_output.A2_B, aes(Time,CL), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CL), size=2.5) + ylab("Concentration") 
plot_B.A2

plot_B.A2_Lung=
  ggplot() +
  geom_line(data  = Fitted_output.A2_B, aes(Time,Clung), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, Clung), size=2.5) + ylab("Concentration") 
plot.A2_Lung

plot_B.A2_Kidney=
  ggplot() +
  geom_line(data  = Fitted_output.A2_B, aes(Time,CK), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CK), size=2.5) + ylab("Concentration") 
plot_B.A2_Kidney

plot_B.A2_Spleen=
  ggplot() +
  geom_line(data  = Fitted_output.A2_B, aes(Time,CS), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CS), size=2.5) + ylab("Concentration") 
plot_B.A2_Spleen


plot_B.A2_Brain=
  ggplot() +
  geom_line(data  = Fitted_output.A2_B, aes(Time,CB), col="firebrick", lwd=2)+
  geom_point(data = Obs.A2    , aes(Time, CB), size=2.5) + ylab("Concentration") 
plot_B.A2_Brain

plot_B.A2_GI=
  ggplot() +
  geom_line(data  = Fitted_output.A2_B, aes(Time,CG), col="firebrick", lwd=2)
plot_B.A2_GI







#---------------------------Fitting with A3 dataset-----------------------------
Obs.A3 <- read.csv(file ="C:/switchdriver/dataset/tk/mouse/R_input_mouse_study1_100nm.csv")   

Fit.Result.A3<- modFit(f=MCcost, p=params2fit, obs=Obs.A3, method ="Nelder-Mead", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"

summary(Fit.Result.A3)                           ## Summary of fit
#exp(Fit.Result$par)                          ## Get the arithmetic value out of the log domain
res.A3=MCcost(Fit.Result.A3$par, obs=Obs.A3)$residuals$res      ## Check the residual for each time points
sum(res.A3^2)                                    ## Total residuals 

# Calculated the model output with fitted parameters

Fitted_output.A3 = pred.mouse(par=Fit.Result.A3$par)

plot.A3=
  ggplot() +
  geom_line(data  = Fitted_output.A3, aes(Time,CL), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3    , aes(Time, CL), size=2.5) + ylab("Concentration") 
plot.A3

plot.A3_Lung=
  ggplot() +
  geom_line(data  = Fitted_output.A3, aes(Time,Clung), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3    , aes(Time, Clung), size=2.5) + ylab("Concentration") 
plot.A3_Lung

plot.A3_Kidney=
  ggplot() +
  geom_line(data  = Fitted_output.A3, aes(Time,CK), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3    , aes(Time, CK), size=2.5) + ylab("Concentration") 
plot.A3_Kidney

plot.A3_Spleen=
  ggplot() +
  geom_line(data  = Fitted_output.A3, aes(Time,CS), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3    , aes(Time, CS), size=2.5) + ylab("Concentration") 
plot.A3_Spleen

plot.A3_Brain=
  ggplot() +
  geom_line(data  = Fitted_output.A3, aes(Time,CB), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3    , aes(Time, CB), size=2.5) + ylab("Concentration") 
plot.A3_Brain

plot.A3_GI=
  ggplot() +
  geom_line(data  = Fitted_output.A3, aes(Time,CG), col="firebrick", lwd=2)
plot.A3_GI
#-------------------fitting with only liver data-------------------
Obs.A3_L <- data.frame(Time=Obs.A3['Time'], 
                       CL=Obs.A3['CL'])

Fit.Result.A3_L<- modFit(f=MCcost, p=params2fit, obs=Obs.A3_L, method ="Nelder-Mead", 
                         control = nls.lm.control(nprint=1)) #"Nelder-Mead"

summary(Fit.Result.A3_L)                           ## Summary of fit
#exp(Fit.Result$par)                          ## Get the arithmetic value out of the log domain
res.A3_L=MCcost(Fit.Result.A3_L$par, obs=Obs.A3)$residuals$res      ## Check the residual for each time points
sum(res.A3_L^2)                                    ## Total residuals 

# Calculated the model output with fitted parameters

Fitted_output.A3_L = pred.mouse(par=Fit.Result.A3_L$par)

plot_L.A3=
  ggplot() +
  geom_line(data  = Fitted_output.A3_L, aes(Time,CL), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3    , aes(Time, CL), size=2.5) + ylab("Concentration") 
plot_L.A3

plot_L.A3_Lung=
  ggplot() +
  geom_line(data  = Fitted_output.A3_L, aes(Time,Clung), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3    , aes(Time, Clung), size=2.5) + ylab("Concentration") 
plot.A3_Lung

plot_L.A3_Kidney=
  ggplot() +
  geom_line(data  = Fitted_output.A3_L, aes(Time,CK), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3    , aes(Time, CK), size=2.5) + ylab("Concentration") 
plot_L.A3_Kidney

plot_L.A3_Spleen=
  ggplot() +
  geom_line(data  = Fitted_output.A3_L, aes(Time,CS), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3    , aes(Time, CS), size=2.5) + ylab("Concentration") 
plot_L.A3_Spleen


plot_L.A3_Brain=
  ggplot() +
  geom_line(data  = Fitted_output.A3_L, aes(Time,CB), col="firebrick", lwd=2)+
  geom_point(data = Obs.A3    , aes(Time, CB), size=2.5) + ylab("Concentration") 
plot.A3_Brain

plot_L.A3_GI=
  ggplot() +
  geom_line(data  = Fitted_output.A3_L, aes(Time,CG), col="firebrick", lwd=2)
plot_L.A3_GI

#------------------------------------------------------------------------------------
