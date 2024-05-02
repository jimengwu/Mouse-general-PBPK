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

# Function to create each individual plot
create_plot <- function(plot_data, y, organ, title) {
  ggplot(data = plot_data, aes(x = Time)) +
    geom_line(aes(y = !!rlang::sym(y)), col = "firebrick", lwd = 1.5) +
    geom_point(data = Obs.A1, aes(y = !!rlang::sym(y)), size = 2, col = "black") +
    labs(title = title, y = paste("Concentration in", organ, "(ng/g)"), x = "Time (h)") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
}
#----------------------------------------------------
#---------------------1.Input Dataset-------------------
##----------------------ZnO-----------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/ZnO/1-ZnO.csv") 
Obs.A1 <- Obs.A1[1:3,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/ZnO/1_GI/'

##---------------------------13 nm full dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_13nm.csv")  
# LONG TERM (WITHIN 24h)
Obs.A1 <- Obs.A1[1:7,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/13nm_nsc_GI/'


##------------------Study2: 13 nm high dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/31501470_Gold dextran.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/31501470_nsc_GI/'



##---------------------------Study2: 13nm low dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/31501470_Gold dextran_2.csv")  
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/31501470_2_nsc_GI/'

##--------------------------Study 3: 61.2nm------------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/29677597_PEG coated Au NPs 61.2nm.csv") 
Obs.A1 <- Obs.A1[1:8,]
PDOSE = Obs.A1$Dose..mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/29677597_61.2_nsc_GI/'


##--------------------------Study 3: 24.3nm------------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/29677597_PEG coated Au NPs 24.3nm.csv") 
Obs.A1 <- Obs.A1[1:8,]
PDOSE = Obs.A1$Dose..mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/29677597_24.3_nsc_GI/'



##---------------------------4 nm full dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_4nm.csv")  
# SHORT TERM (WITHIN 24h)
Obs.A1 <- Obs.A1[1:7,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/4nm_nsc_GI/'
##---------------------------100 nm short dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_100nm.csv")  
# SHORT TERM (WITHIN 24h)
Obs.A1 <- Obs.A1[1:3,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","Clung")] #100nm only three organs
folder = "plots/100nm_short_nsc_GI/"

##----------------------done------------




##--------------------------Study 3: 6.2nm------------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/29677597_PEG_coated_AuNPs_6.2nm.csv") 
Obs.A1 <- Obs.A1[1:8,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/29677597_6.2_nsc_GI/'

##------------------GO: Study1: 20 nm high dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/GO/21162527_125I-NGS-PEG 10-30nm.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/GO/2116257_GI/'



##---------------------------100 nm full dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_100nm.csv")  
Obs.A1 <- Obs.A1[1:7,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","Clung")] #100nm only three organs
folder = 'plots/100nm_nsc_GI/'

##--------------------------Study 3: 42.5nm------------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/29677597_PEG coated Au NPs 42.5nm.csv") 
Obs.A1 <- Obs.A1[1:8,]
PDOSE = Obs.A1$Dose..mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/29677597_42.5_nsc_GI/'
##----------------------bad------------------
##------------------TiO2: Study1: 20 nm high dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/TiO2/1-TiO2.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/TiO2/1_GI/'


#---------------------1.Input Dataset-------------------



#---------------------2. Build mrgsolve-based PBPK Model-------
mod <- mcode ("mouse_PBPK", mousePBPK.code)

set.seed(5)

pred.mouse <- function(pars) {

    ## Get out of log domain
    pars %<>% lapply(exp) # todo: important to have because we cannot have nagetive kinetic value

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
    tsamp=tgrid(0,max(Obs.A1$Time),1)     ## Simulation time 24*7 hours (180 days)

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


pred.mouse.oral <- function(pars) {
  
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
  tsamp=tgrid(0,max(Obs.A1$Time),0.05)     ## Simulation time 24*7 hours (180 days)
  
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

# --------------2.1 initial parameters--------
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

#--------------------2. initial parameter sensitivity analysis------------

Sens.init <- sensFun(func = pred.mouse,parms = params.init)
Sens.init_abs <- Sens.init
#sensitivity value was calculated with the time range, therefore the value is 
# increase if the modeled time range is longer

# Identify the index of the 'var' column
var_index <- which(names(Sens.init) == "var")

# Exclude the 'var' column and apply the absolute function to make all values positive
Sens.init_abs[,-var_index] <- abs(Sens.init[,-var_index])

df_Sens.init_abs=summary(Sens.init_abs)
ranked_df_Sens_abs.init <-df_Sens.init_abs[order(df_Sens.init_abs$Mean, decreasing = TRUE), ]

select_par = rownames(ranked_df_Sens_abs.init [ranked_df_Sens_abs.init $Mean>0.5,])
select_par

# Print selected parameters
print(select_par)

params2fit = params.init[select_par]
graphics.off()
png(paste0(paste0(folder,"init_sens/"),"ranked_df_Sens_abs_plot_init.png"),width=3500,height=2800,res = 300)
plot(ranked_df_Sens_abs.init)
dev.off()

#-------------FIRST FITTING WITH EVERY parameters-----------------
MCcost<-function (pars, obs){
  out<- pred.mouse(pars)
  cost<- modCost(model=out,obs=obs,weight='std',x="Time")
  return(cost)
}


Fit.Result.A1<- modFit(f=MCcost, p=params.init, obs=Obs.A1, method ="Nelder-Mead", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"
summary(Fit.Result.A1)                           ## Summary of fit


res.A1=MCcost(Fit.Result.A1$par, obs=Obs.A1)$residuals$res     ## Check the residual for each time points
sum(res.A1^2)                                    ## Total residuals 

# Calculated the model output with fitted parameters

Fitted_output.A1 = pred.mouse(par=Fit.Result.A1$par)

# Create individual plots
plot_liver.a1   <- create_plot(Fitted_output.A1, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a1  <- create_plot(Fitted_output.A1, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a1    <- create_plot(Fitted_output.A1, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a1  <- create_plot(Fitted_output.A1, y = "CS", organ = "Spleen", title = "Spleen")

# Arrange the plots together
combined_mod_fit_plot.a1 <- grid.arrange(plot_liver.a1, plot_kidney.a1, plot_lung.a1, plot_spleen.a1, 
                                         ncol = 2, nrow = 2)

#--------plot fitted with observed----------


#-------------------------first fitted parameters value sensitivity---------------------------
Sens <- sensFun(func = pred.mouse,parms = Fit.Result.A1$par)

max_val <- max(Sens[3:ncol(Sens)]) # first two columns are x and var name
max_val_scientific <- format(max_val, scientific = TRUE)
min_val <- min(Sens[3:ncol(Sens)])
min_val_scientific <- format(min_val, scientific = TRUE)

for (i in 3:length(Sens)) {
  
  # Plot using ggplot
  p1.m <- ggplot(Sens, aes(x = Sens[,1], y = Sens[, i], color = var)) +
    geom_line(size = 1.5) + 
    ylim(min_val,max_val)+
    labs(x = "Time(h)", y = "normalised sensitivities of model output to parameter",title=names(params.init)[i-2]) +
    theme_minimal()+
    theme(
      axis.line = element_line(linewidth = 1),  # Increase axis width
      axis.text = element_text(size = 20),
      axis.title = element_text(size=20),
      plot.title = element_text(size = 20,  hjust = 0.1, vjust = -1),  # Adjust title position
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 16),  # Increase legend text size
      legend.title = element_blank()   # Increase legend title size
    )
  
  # Save or print combined_plot as desired
  ggsave(paste0("combined_init_sensitivity_plot_", names(params.init)[i-2], ".png"), 
         path = paste0(folder,"init_sens_fitted/"), p1.m, width = 10, height = 6,dpi = 300)
  # Clear the plot list for the next iteration
  plot_list <- list()
}


Sens_abs <- Sens
#sensitivity value was calculated with the time range, therefore the value is 
# increase if the modeled time range is longer

# Identify the index of the 'var' column
var_index <- which(names(Sens) == "var")

# Exclude the 'var' column and apply the absolute function to make all values positive
Sens_abs[,-var_index] <- abs(Sens[,-var_index])

df_Sens_abs=summary(Sens_abs)


# Rank the data frame based on the 'Score' column
ranked_df_Sens_abs <-df_Sens_abs[order(df_Sens_abs$Mean, decreasing = TRUE), ]

graphics.off()
png(paste0(paste0(folder,"init_sens_fitted/"),"ranked_df_Sens_abs_plot_fitted.png"),width=3500,height=2800,res = 300)
plot(ranked_df_Sens_abs)
dev.off()


select_par = rownames(ranked_df_Sens_abs[ranked_df_Sens_abs$Mean>0.5,])

params2fit = params.init[select_par]

Fit.Result.A2<- modFit(f=MCcost, p=params2fit, obs=Obs.A1, method ="Nelder-Mead", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"
summary(Fit.Result.A1)                           ## Summary of fit

# why sensitivity analysis give that dlc_lung is important?
res.A2=MCcost(Fit.Result.A2$par, obs=Obs.A1)$residuals$res     ## Check the residual for each time points
sum(res.A2^2)                                    ## Total residuals 

# Calculated the model output with fitted parameters

Fitted_output.A2 = pred.mouse(par=Fit.Result.A2$par)


# Create individual plots
plot_liver.a2   <- create_plot(Fitted_output.A2, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a2  <- create_plot(Fitted_output.A2, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a2    <- create_plot(Fitted_output.A2, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a2  <- create_plot(Fitted_output.A2, y = "CS", organ = "Spleen", title = "Spleen")

# Arrange the plots together
combined_mod_fit_plot.a2 <- grid.arrange(plot_liver.a2, plot_kidney.a2, plot_lung.a2, plot_spleen.a2, 
                                         ncol = 2, nrow = 2)

ggsave(paste0("sens_fitted.png"), 
       path = paste0(folder,"init_sens_fitted/"), combined_mod_fit_plot.a2, width = 10, height = 6,dpi = 300)

ggsave(paste0("all_fitted.png"), 
       path = paste0(folder,"init_sens_fitted/"), combined_mod_fit_plot.a1, width = 10, height = 6,dpi = 300)

#-----------------check the parameter sensitivity -------




# Initialize a list to store the plots
plot_list <- list()
#----------Draw the sensitivity line--------

parameters__ssens = params.init
R = pred.mouse(params.init)
folds <- seq(from = 10e-2, to = 10e18, length.out = 10)
folds = c(10e5)
lsa <- data.frame(Time = numeric(),
                  Organ = character(),
                  Fold = numeric(),
                  Parameter = character(),
                  predicted = numeric(),
                  initial_value = numeric())

for (organ in c("Liver", "Kidney", "Spleen", "Lung")) {
  cat("organ:", organ, "\n")  # Debug output
  LSA_data <- data.frame()
  for (i in 24:length(parameters__ssens)) {
    for (fold in folds) {
      pars.mouse.double <- log(c(exp(parameters__ssens[i]) * fold, exp(parameters__ssens[-i]))) 
      Rnew.double <- pred.mouse(pars.mouse.double)
      
      #pars.mouse.half <- log(c(exp(parameters__ssens[i]) / 100, exp(parameters__ssens[-i]))) 
      #Rnew.half <- pred.mouse(pars.mouse.half)
      
      if (organ == "Liver") {
        LSA_mouse_double <- as.numeric(unlist(Rnew.double$CL))
        #LSA_mouse_half <- as.numeric(unlist(Rnew.half$CL))
        initial_value = R$CL
      } else if (organ == "Kidney") {
        LSA_mouse_double <- as.numeric(unlist(Rnew.double$CK))
        #LSA_mouse_half <- as.numeric(unlist(Rnew.half$CK))
        initial_value = R$CK
      } else if (organ == "Spleen") {
        LSA_mouse_double <- as.numeric(unlist(Rnew.double$CS))
        #LSA_mouse_half <- as.numeric(unlist(Rnew.half$CS))
        initial_value = R$CS
      } else if (organ == "Lung") {
        LSA_mouse_double <- as.numeric(unlist(Rnew.double$Clung))
        #LSA_mouse_half <- as.numeric(unlist(Rnew.half$Clung))
        initial_value = R$Clung
      }
      
      combined_data <- data.frame(Organ = organ,
                                  Time = R$Time,
                                  parameter=names(parameters__ssens)[i],
                                  #half = LSA_mouse_half,
                                  predicted = LSA_mouse_double,
                                  initial_value = initial_value,
                                  Fold = factor(fold))
      
      lsa <- rbind(lsa, combined_data)
      
    }}

}

library(ggplot2)
# Filter data for the liver
lsa_liver <- subset(lsa, Organ == "Kidney")

lsa_liver <- subset(lsa_liver, parameter == "DLC_Lung")

# Plotting LSA data
ggplot(lsa_liver, aes(x = Time, y = predicted, color = Fold)) +
  geom_line(linetype="dashed") +
  
  geom_point(data = Obs.A1, aes(x = Time, y = Clung), shape = 16, size = 1, color = "red") +
  #geom_point(aes(y = initial_value), color = "blue",size=1,alpha=0.1) +
  facet_wrap(parameter ~ Organ, scales = "free_y",ncol=5) +
  labs(title = "LSA Data",
       x = "Time",
       y = "Concentration (ng/g)") +
  theme_minimal()



library(ggplot2)

# Filter data for each organ
lsa_organs <- split(lsa, lsa$Organ)

# Function to plot and save each organ's plot
plot_and_save <- function(organ_data, filename) {
  gg <- ggplot(organ_data, aes(x = Time, y = predicted, color = Fold)) +
    geom_line() +
    geom_point(aes(y = initial_value), color = "blue", size = 0.5) +
    facet_wrap(parameter ~ Organ, scales = "free_y") +
    labs(title = "LSA Data", x = "Time", y = "Concentration (ng/g)") +
    theme_minimal()
  
  ggsave(filename, gg, width = 10, height = 6)
}

# Save plots for all organs
for (organ_name in names(lsa_organs)) {
  plot_and_save(lsa_organs[[organ_name]], paste0(folder, "LSA_", organ_name, ".png"))
}




#----------Draw the sensitivity line 2--------

parameters__ssens = params.init
R = pred.mouse(params.init)
folds <- seq(from = 10e-2, to = 10e5, length.out = 10)
folds = c(10e5)
lsa <- data.frame(Time = numeric(),
                  Organ = character(),
                  Fold = numeric(),
                  Parameter = character(),
                  predicted = numeric(),
                  initial_value = numeric())

for (organ in c("Liver", "Kidney", "Spleen", "Lung")) {
  cat("organ:", organ, "\n")  # Debug output
  LSA_data <- data.frame()
  for (i in 24:length(parameters__ssens)) {
    for (fold in folds) {
      pars.mouse.double <- c(parameters__ssens[i]* fold, (parameters__ssens[-i]))
      Rnew.double <- pred.mouse(pars.mouse.double)
      
      #pars.mouse.half <- log(c(exp(parameters__ssens[i]) / 100, exp(parameters__ssens[-i]))) 
      #Rnew.half <- pred.mouse(pars.mouse.half)
      
      if (organ == "Liver") {
        LSA_mouse_double <- as.numeric(unlist(Rnew.double$CL))
        #LSA_mouse_half <- as.numeric(unlist(Rnew.half$CL))
        initial_value = R$CL
      } else if (organ == "Kidney") {
        LSA_mouse_double <- as.numeric(unlist(Rnew.double$CK))
        #LSA_mouse_half <- as.numeric(unlist(Rnew.half$CK))
        initial_value = R$CK
      } else if (organ == "Spleen") {
        LSA_mouse_double <- as.numeric(unlist(Rnew.double$CS))
        #LSA_mouse_half <- as.numeric(unlist(Rnew.half$CS))
        initial_value = R$CS
      } else if (organ == "Lung") {
        LSA_mouse_double <- as.numeric(unlist(Rnew.double$Clung))
        #LSA_mouse_half <- as.numeric(unlist(Rnew.half$Clung))
        initial_value = R$Clung
      }
      
      combined_data <- data.frame(Organ = organ,
                                  Time = R$Time,
                                  parameter=names(parameters__ssens)[i],
                                  #half = LSA_mouse_half,
                                  predicted = LSA_mouse_double,
                                  initial_value = initial_value,
                                  Fold = factor(fold))
      
      lsa <- rbind(lsa, combined_data)
      
    }}
  
}

library(ggplot2)
# Filter data for the liver
lsa_liver <- subset(lsa, Organ == "Lung")

lsa_liver <- subset(lsa_liver, parameter == "DLC_Lung")

# Plotting LSA data
ggplot(lsa_liver, aes(x = Time, y = predicted, color = Fold)) +
  geom_line(linetype="dashed") +
  geom_line(aes(y = initial_value), color = "blue",size=1,alpha=0.1) +
  facet_wrap(parameter ~ Organ, scales = "free_y",ncol=5) +
  labs(title = "LSA Data",
       x = "Time",
       y = "Concentration (ng/g)") +
  theme_minimal()



library(ggplot2)

# Filter data for each organ
lsa_organs <- split(lsa, lsa$Organ)

# Function to plot and save each organ's plot
plot_and_save <- function(organ_data, filename) {
  gg <- ggplot(organ_data, aes(x = Time, y = predicted, color = Fold)) +
    geom_line() +
    geom_point(aes(y = initial_value), color = "blue", size = 0.5) +
    facet_wrap(parameter ~ Organ, scales = "free_y") +
    labs(title = "LSA Data", x = "Time", y = "Concentration (ng/g)") +
    theme_minimal()
  
  ggsave(filename, gg, width = 10, height = 6)
}

# Save plots for all organs
for (organ_name in names(lsa_organs)) {
  plot_and_save(lsa_organs[[organ_name]], paste0(folder, "LSA_", organ_name, ".png"))
}



#----------------------rsquared for each organ prediction-------
# Initialize an empty data frame to store R-squared values
rsquared_organ <- data.frame(Organ = character(),
                            Parameter = character(),
                            Fold = numeric(),
                            R_squared = numeric())

# Define the range of folds to consider
folds <- seq(from = 1e-2, to = 1e2, length.out = 100)
for (organ in c("Liver", "Kidney", "Spleen", "Lung")) {
  for (i in 1:length(parameters__ssens)) {
    cat("i:", i, "\n")  # Debug output
    # Extract parameter for sensitivity analysis
    parameter <- names(parameters__ssens)[i]
    
    # Generate predictions for each fold of the parameter
    for (fold in folds) {
      # Calculate the new parameter values
      pars.mouse.fold <- log(c(exp(parameters__ssens[i]) * fold, exp(parameters__ssens[-i]))) 
      Rnew.fold <- pred.mouse(pars.mouse.fold)
      
      # Filter Rnew.fold to include only the rows where Time matches the Obs dataset
      Rnew.fold <- Rnew.fold[Rnew.fold$Time %in% Obs.A1$Time, ]
      # Extract observed values based on the organ
      observed <- switch(organ,
                         "Liver" = Obs.A1$CL,
                         "Kidney" = Obs.A1$CK,
                         "Spleen" = Obs.A1$CS,
                         "Lung" = Obs.A1$Clung)
      
      # Extract predicted values based on the organ
      predicted <- switch(organ,
                          "Liver" = as.numeric(unlist(Rnew.fold$CL)),
                          "Kidney" = as.numeric(unlist(Rnew.fold$CK)),
                          "Spleen" = as.numeric(unlist(Rnew.fold$CS)),
                          "Lung" = as.numeric(unlist(Rnew.fold$Clung)))
      
      # Calculate R-squared
      r_squared <- cor(observed, predicted)^2
      
      # Append the results to the data frame
      rsquared_organ <- rbind(rsquared_organ, data.frame(Organ = organ,
                                                       Parameter = parameter,
                                                       Fold = fold,
                                                       R_squared = r_squared))
      
    }
  }
}

# Print or save the data frame with R-squared values for further analysis
print(rsquared_organ)
library(ggplot2)

# Plotting R-squared values against fold for each parameter
p = ggplot(rsquared_organ, aes(x = Fold, y = R_squared, color = Organ)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.6, ymax = 1), fill = "gray80", alpha = 0.1) +
  geom_point(size=1) +
  geom_line(size=1) +
  facet_wrap(~Parameter, scales = "free") +
  labs(title = "R-squared Values vs. Fold for Each Parameter",
       x = "Fold",
       y = "R-squared",
       color = "Organ")+
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border around the plot
  )
# Save the plot as a file
ggsave(paste0(folder,"rsquared_organ_plot.png"), plot = p, width = 20, height = 10, units = "in")
#--------------------done---


#-----------------for the total rsquared----
# Initialize an empty data frame to store R-squared values
rsquared_data <- data.frame(
                            Parameter = character(),
                            Fold = numeric(),
                            R_squared = numeric())



# Define the range of folds to consider
folds <- seq(from = 1e-5, to = 1e5, length.out = 100)

for (i in 1:length(parameters__ssens)) {
    cat("i:", i, "\n")  # Debug output
    # Extract parameter for sensitivity analysis
    parameter <- names(parameters__ssens)[i]
    
    # Generate predictions for each fold of the parameter
    for (fold in folds) {
      # Calculate the new parameter values
      pars.mouse.fold <- log(c(exp(parameters__ssens[i]) * fold, exp(parameters__ssens[-i]))) 
      Rnew.fold <- pred.mouse(pars.mouse.fold)
      
      # Filter Rnew.fold to include only the rows where Time matches the Obs dataset
      Rnew.fold <- Rnew.fold[Rnew.fold$Time %in% Obs.A1$Time, ]
      
      # Extract observed values based on the organ (assuming you have a single organ's data)
      observed <- Obs.A1
      
      # Extract predicted values from Rnew.fold (adjust based on your prediction column name)
      predicted <- as.numeric(unlist(Rnew.fold[-1]))
      
      observed <- as.numeric(unlist(Obs.A1[-1]))
      
      
      # Calculate R-squared
      r_squared <- cor(observed, predicted)^2
      
      # Append the results to the data frame
      rsquared_data <- rbind(rsquared_data, data.frame(Parameter = parameter,
                                                       Fold = fold,
                                                       R_squared = r_squared))
      
    }
}

# Plotting R-squared values against fold for each parameter

ggplot(rsquared_data, aes(x = Fold, y = R_squared)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.6, ymax = 1), fill = "gray80", alpha = 0.3) +
  geom_point(size=1) +
  geom_line() +
  facet_wrap(~Parameter, scales = "free", ncol = 6) +
  labs(title = "R-squared Values vs. Fold for Each Parameter",
       x = "Fold",
       y = "R-squared") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border around the plot
  )


#---------problem is dlc_lung------

select_par = "DLC_Lung"

params2fit = params.init[select_par]

Fit.Result.A3<- modFit(f=MCcost, p=params2fit, obs=Obs.A1, method ="Nelder-Mead", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"
Fitted_output.A3 = pred.mouse(par=Fit.Result.A3$par)

Fitted_output.A0 = pred.mouse(par=params.init)

# Create individual plots
plot_liver.a0   <- create_plot(Fitted_output.A0, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a0  <- create_plot(Fitted_output.A0, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a0    <- create_plot(Fitted_output.A0, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a0  <- create_plot(Fitted_output.A0, y = "CS", organ = "Spleen", title = "Spleen")
# Arrange the plots together
combined_mod_fit_plot.a0 <- grid.arrange(plot_liver.a0, plot_kidney.a0, plot_lung.a0, plot_spleen.a0, 
                                         ncol = 2, nrow = 2)




# Create individual plots
plot_liver.a3   <- create_plot(Fitted_output.A3, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a3  <- create_plot(Fitted_output.A3, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a3    <- create_plot(Fitted_output.A3, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a3  <- create_plot(Fitted_output.A3, y = "CS", organ = "Spleen", title = "Spleen")


# Arrange the plots together
combined_mod_fit_plot.a3 <- grid.arrange(plot_liver.a3, plot_kidney.a3, plot_lung.a3, plot_spleen.a3, 
                                         ncol = 2, nrow = 2)

#------------------------plot with the changing parameters-----------------

library(ggplot2)

# Function to generate fold colors with varying alpha values
generate_fold_colors <- function(base_color, num_folds) {
  # Define the base colors (purple and yellow)
  base_colors <- c("purple", "yellow")
  # Generate a color palette from purple to yellow
  color_palette <- colorRampPalette(base_colors)(num_folds)
  # Generate a series of colors with varying alpha values
  sapply(seq(1, 1, length.out = num_folds), function(alpha) {
    sapply(color_palette, function(color) {
      alpha(color, alpha)
    })
  })
}

# Define the range of folds to consider
folds <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 50, 100, 500, 1000)

# Initialize an empty data frame to store predicted values
predicted_data <- data.frame(Organ = character(),
                             Parameter = character(),
                             Fold = numeric(),
                             Time = numeric(),
                             Predicted = numeric())

# Loop through organs and parameters
for (organ in c("Liver", "Kidney", "Spleen", "Lung")) {
  for (i in 1:length(parameters__ssens)) {
    # Extract parameter for sensitivity analysis
    parameter <- names(parameters__ssens)[i]
    
    # Generate predictions for each fold of the parameter
    for (fold in folds) {
      # Calculate the new parameter values
      pars.mouse.fold <- log(c(exp(parameters__ssens[i]) * fold, exp(parameters__ssens[-i]))) 
      Rnew.fold <- pred.mouse(pars.mouse.fold)
      
      # Extract predicted values based on the organ
      predicted_values <- switch(organ,
                                 "Liver" = as.numeric(unlist(Rnew.fold$CL)),
                                 "Kidney" = as.numeric(unlist(Rnew.fold$CK)),
                                 "Spleen" = as.numeric(unlist(Rnew.fold$CS)),
                                 "Lung" = as.numeric(unlist(Rnew.fold$Clung)))
      
      # Create a data frame for the predicted values
      predicted_df <- data.frame(Organ = organ,
                                 Parameter = parameter,
                                 Fold = fold,
                                 Time = Rnew.fold$Time,
                                 Predicted = predicted_values)
      
      # Append the results to the data frame
      predicted_data <- rbind(predicted_data, predicted_df)
    }
  }
}

# Plotting function
plot_organ <- function(organ_name, predicted_data_subset, fold_colors) {
  # Plotting
  ggplot(predicted_data_subset, aes(x = Time, y = Predicted, color = as.factor(Fold))) +
    geom_line() +
    scale_color_manual(values = fold_colors) +
    facet_wrap(Parameter ~ ., scales = "free_y", ncol = 6) +
    geom_point(data = Obs.A1, aes(x = Time, y = CL), shape = 16, size = 1, color = "red") +
    labs(title = paste("Predicted Values for Organ", organ_name),
         x = "Time",
         y = "Predicted Value",
         color = "Fold") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border around the plot
    )
}

# Generate fold colors
num_folds <- length(folds)
fold_colors <- generate_fold_colors("black", num_folds)

# Plot for Liver
predicted_data_subset <- predicted_data %>% filter(Fold != 1) %>% filter(Organ == "Liver")
#ggsave(plot_organ("Liver", predicted_data_subset, fold_colors), "Liver_plot.png")

# Plot for Kidney
predicted_data_subset <- predicted_data %>% filter(Fold != 1) %>% filter(Organ == "Kidney")
#ggsave(plot_organ("Kidney", predicted_data_subset, fold_colors), "Kidney_plot.png")

# Plot for Spleen
predicted_data_subset <- predicted_data %>% filter(Fold != 1) %>% filter(Organ == "Spleen")
#ggsave(plot_organ("Spleen", predicted_data_subset, fold_colors), "Spleen_plot.png")


# Plot for Lung
predicted_data_subset <- predicted_data %>% filter(Fold != 1) %>% filter(Organ == "Lung")
#ggsave(plot_organ("Lung", predicted_data_subset, fold_colors), "Lung_plot.png")


#-------------circle sensitivity plot-------
params.cali = Fit.Result.A1$par
## Create the matrix for normalized sensitivity coefficient data, 4 columns each represent one organ
NSC_mouse   = matrix(nrow=length(Fit.Result.A1$par),ncol=4)
NSC_mouse_end = NSC_mouse

percentage = 0.01
R = pred.mouse(Fit.Result.A1$par)
tend              <- max(Obs.A1$Time)
for (i in 1:length(Fit.Result.A1$par)) {
  NSC_mouse_T   = matrix(nrow=length(params.cali),
                         ncol=length(R$Time))
  for (j in 2:5){
    # Each cycle, generate a new value of parameter i (e.g., 10.0a), 
    # and delete parameter i, so that you can proceed to the next parameter i+1
    params.cali.new      <- log(c(exp(params.cali[i])*(1 + percentage),exp(params.cali[-i]))) 
    Rnew                <- pred.mouse(params.cali.new)
    delta.P.mouse       <- exp(params.cali[i])/(exp(params.cali[i])*percentage) # is exp here needed?
    
    ## Estimated the AUC
    
    # get the NSC for the whole time range
    delta.AUC.CLt.mouse.T = (Rnew %>% select (names(R)[j]) 
                             - R %>% select (names(R)[j]))
    
    NSC_mouse_T [i,] <-as.numeric(unlist(delta.AUC.CLt.mouse.T/ 
                                           (R %>% select (j))) 
                                  *delta.P.mouse)
    
    
    #get the NSC median value for the whole time range
    #NSC_mouse   [i, 1]      <- as.numeric((delta.AUC.CLt.mouse/Mouse.AUC.CLt.ori) * delta.P.mouse)
    NSC_mouse   [i, j-1]      <- NSC_mouse_T [i,which(Rnew$Time == 1)]
    NSC_mouse_end [i, j-1]    <- NSC_mouse_T [i,which(Rnew$Time == tend)]
  }}



colnames (NSC_mouse)      = c("NSC_AUC_CKt","NSC_AUC_CLt","NSC_AUC_CSt","NSC_AUC_CLut") 
colnames (NSC_mouse_end)  =  c("NSC_AUC_CKt","NSC_AUC_CLt","NSC_AUC_CSt","NSC_AUC_CLut") 

rownames(NSC_mouse)       = names(Fit.Result.A1$par)
rownames(NSC_mouse_T)     = names(Fit.Result.A1$par)
rownames(NSC_mouse_end)   = names(Fit.Result.A1$par)

