# ---------------------------------------------------------
# Author: Jimeng Wu
# Email: jimeng.wu@empa.ch
# Date: 2025-02-18
# Description: initial calibration of the PBPK model for each NP case, the initial parameter was set the same, 
# then using the modFit function to fit the observation (concentration-time data points) to the model, getting the fitted parameters 
# ---------------------------------------------------------




## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' package
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
library(Metrics)
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

#--------------------- Build mrgsolve-based PBPK Model-------
mod <- mcode ("mouse_PBPK", mousePBPK.code)

set.seed(5)

if (np.name == "FeO: Study2_41nm_4mg/kg") {
  tstep <- 0.5/60 # 2.5/60 in the time step
} else if (np.name %in% c("GO: Study2_243nm_1mg/kg",
                          "GO: Study2_914nm_1mg/kg_all","GO: Study2_914nm_1mg/kg_w/o_CS")) {
  # Assign the appropriate value for time step
  tstep <- 1/60   # Fill in the appropriate value here 2/60 & 5/60 in the time points
} else {
  tstep <- min(1, min(Obs.df$Time))
}


#-----------iv------
pred.mouse.iv <- function(pars) {
  
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
  tsamp=tgrid(0,max(Obs.df$Time),tstep)     ## Simulation time 24*7 hours (180 days)
  
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

#-----------oral----
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
  tsamp=tgrid(0,max(Obs.df$Time),tstep)     ## Simulation time 24*7 hours (180 days)
  
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

if (pathway == "intraperitoneal injection") {
  pred.mouse <- pred.mouse.oral
} else if (pathway == "intravenous injection") {
  pred.mouse <- pred.mouse.iv
}


MCcost<-function (pars, obs){
  out<- pred.mouse(pars)
  cost<- modCost(model=out,obs=obs,weight='mean',x="Time")
  return(cost)
}

#-------------- 1. choose the sensitive parameters for fitting -------------
#--------------draw the sensitive line one by one--------
Sens.init <- sensFun(func = pred.mouse,parms = params.init)

create_sensitivity_plots(Sens.init, names(params.init), 
                         paste0(folder,"init_sens/"))

# Identify the index of the 'var' column

var_index <- which(names(Sens.init) == "var")

# Exclude the 'var' column and apply the absolute function to make all values positive
Sens.init_abs <- Sens.init
Sens.init_abs[,-var_index] <- abs(Sens.init[,-var_index])

df_Sens.init_abs=summary(Sens.init_abs)

graphics.off()
png(paste0(paste0(folder,"init_sens/"),"ranked_df_Sens_abs_plot_init.png"),width=3500,height=2800,res = 300)
plot(df_Sens.init_abs[order(df_Sens.init_abs$Mean, decreasing = FALSE), ])
dev.off()

# Rank the parameters based on the sensitivity values

ranked_df_Sens_abs.init <-df_Sens.init_abs[order(df_Sens.init_abs$Mean, decreasing = TRUE), ]
select_par = rownames(ranked_df_Sens_abs.init [ranked_df_Sens_abs.init $Mean > 0.5,])
select_par



params2fit = params.init[select_par]

Fit.Result.A1<- modFit(f=MCcost, p=params2fit, obs=Obs.df, method ="Port", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"

res.A1=MCcost(Fit.Result.A1$par, obs=Obs.df)$residuals$res     ## Check the residual for each time points
sum(res.A1^2)      

Fitted_output.A1 = pred.mouse(par=Fit.Result.A1$par)

# Create individual plots
plot_liver.a1   <- create_plot(Fitted_output.A1, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a1  <- create_plot(Fitted_output.A1, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a1    <- create_plot(Fitted_output.A1, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a1  <- create_plot(Fitted_output.A1, y = "CS", organ = "Spleen", title = "Spleen")


# Arrange the plots together
combined_mod_fit_plot.a1 <- grid.arrange(plot_spleen.a1, 
                                         ncol = 2, nrow = 2,
                                         top = "Fitted with sensitive parameters")

combined_mod_fit_plot.a1 <- grid.arrange(plot_liver.a1, plot_lung.a1, plot_spleen.a1, 
                                         ncol = 2, nrow = 2,
                                         top = "Fitted with sensitive parameters")

combined_mod_fit_plot.a1 <- grid.arrange(plot_liver.a1, plot_kidney.a1, 
                                         plot_lung.a1, 
                                         ncol = 2, nrow = 2,
                                         top = "Fitted with sensitive parameters")

combined_mod_fit_plot.a1 <- grid.arrange(plot_liver.a1, plot_kidney.a1, 
                                         plot_lung.a1, plot_spleen.a1, 
                                         ncol = 2, nrow = 2)

# Save the combined plot
#ggsave(paste0(folder, "mod_fit/mod_fit_combined_v1_p.png"), combined_mod_fit_plot.a1,
#       width = 14, height = 10)

params.mod.a1 = params.init
# Iterate over the names of the new list
for (name in names(Fit.Result.A1$par)) {
  # Replace the element in the original list with the corresponding element from the new list
  params.mod.a1[[name]] <- Fit.Result.A1$par[[name]]
}

#write.csv(params.mod.a1, file = paste0(folder,'mod_fit/params_fitted_v1_p.csv'))


#----------2. use all parameters as fitting----------
Fit.Result.A0<- modFit(f=MCcost, p=params.init, obs=Obs.df, method ="Port", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"
res.A0=MCcost(Fit.Result.A0$par, obs=Obs.df)$residuals$res     ## Check the residual for each time points
sum(res.A0^2)      

Fitted_output.A0 = pred.mouse(par=Fit.Result.A0$par)

# Create individual plots
plot_liver.a0   <- create_plot(Fitted_output.A0, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a0  <- create_plot(Fitted_output.A0, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a0    <- create_plot(Fitted_output.A0, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a0  <- create_plot(Fitted_output.A0, y = "CS", organ = "Spleen", title = "Spleen")

# Arrange the plots together
combined_mod_fit_plot.a0 <- grid.arrange(plot_liver.a0, plot_lung.a0, plot_spleen.a0, 
                                         ncol = 2, nrow = 2,
                                         top = "Fitted with every parameters")
combined_mod_fit_plot.a0 <- grid.arrange(plot_liver.a0, plot_kidney.a0, plot_lung.a0, 
                                         ncol = 2, nrow = 2,
                                         top = "Fitted with every parameters")
combined_mod_fit_plot.a0 <- grid.arrange(plot_liver.a0, plot_kidney.a0, plot_lung.a0, plot_spleen.a0, 
                                         ncol = 2, nrow = 2,
                                         top = "Fitted with every parameters")


# Save the combined plot
ggsave(paste0(folder, "mod_fit/mod_fit_combined_v1_0_p.png"), combined_mod_fit_plot.a0, width = 14, height = 10)


params.mod.a0 = params.init
# Iterate over the names of the new list
for (name in names(Fit.Result.A0$par)) {
  # Replace the element in the original list with the corresponding element from the new list
  params.mod.a0[[name]] <- Fit.Result.A0$par[[name]]
}

write.csv(params.mod.a0, file = paste0(folder,'mod_fit/params_fitted_v1_0_p.csv'))



#----------------------RMSE for each organ prediction-------
# Initialize an empty data frame to store R-squared values
rmse_organ <- data.frame(Organ = character(),
                         Parameter = character(),
                         Fold = numeric(),
                         R_squared = numeric(),
                         RMSE = numeric())
rmse_data_tot <- data.frame(
  Parameter = character(),
  Fold = numeric(),
  R_squared = numeric(),
  RMSE = numeric())

# Define the range of folds to consider
fold1 <- seq(from = 1e-3, to = 1, length.out = 100)
fold2 <- seq(from = 1, to = 1e3, length.out = 100)

folds = c(fold1, fold2)



rmse_result = calculate_rmse(names(params.init),params.init,folds, pred.mouse, 
                             Obs.df, rmse_data_tot, rmse_organ)

rmse_data_tot = rmse_result$rmse_data_tot
rmse_organ = rmse_result$rmse_organ

# Print or save the data frame with R-squared values for further analysis



#------------plotting rmse-----------------------
# Plotting R-squared values against fold for each parameter

# Plotting R-squared values against fold for each parameter and each organ,
# four organ results were separated into 4 curves
rmse_plot = ggplot(rmse_organ, aes(x = Fold, y = RMSE, color = Organ)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "gray80", alpha = 0.1) +
  geom_point(size=1) +
  geom_line(size=1) +
  scale_x_log10() +
  facet_wrap(~Parameter, scales = "free",ncol=6) +
  labs(
       x = "Fold",
       y = "RMSE",
       color = "Organ")+
  #coord_cartesian(ylim = c(0, 1)) +
  theme_minimal()+
  theme(
   axis.title = element_text(size=18),
   axis.text.x = element_text(angle = 90),
   axis.text = element_text(size=14),
   strip.text = element_text(size=14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border around the plot
  )
rmse_plot
# Save the plot as a file
#ggsave(paste0(folder,"init_sens_v2/","rmse_organ_plot.png"), plot = rmse_plot, width = 20, height = 15, units = "in")



# Plotting R-squared values of the total organ amount against fold for each parameter 

rmse_tot = ggplot(rmse_data_tot, aes(x = Fold, y = RMSE)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.6, ymax = 1), fill = "gray80", alpha = 0.3) +
  geom_point(size=1) +
  geom_line() +
  scale_x_log10() +
  facet_wrap(~Parameter, scales = "free", ncol = 8) +
  labs(title = "Total RMSE Values vs. Fold for Each Parameter",
       x = "Fold",
       y = "RMSE",size=19) +
  theme_minimal() +
  coord_cartesian(ylim = c(min(rmse_data_tot$RMSE), max(rmse_data_tot$RMSE))) +
  theme(
    axis.text.x = element_text(angle=90,size=12), strip.text = element_text(size = 14),
    axis.text.y = element_text(size = 12), 
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border around the plot
  )


# Save the plot as a file
#ggsave(paste0(folder,"init_sens_v2/","rmse_tot_plot.png"), plot = rmse_tot, width = 20, height = 10, units = "in")


