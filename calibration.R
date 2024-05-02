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
library(Metrics)
source("helper_functions.R")
source("Mouse_PBPK.R")

#------------------TiO2: Study1: 20 nm high dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/TiO2/1-TiO2.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/TiO2/1_GI/'


#------------------GO: Study1: 20 nm high dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/GO/21162527_125I-NGS-PEG 10-30nm.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/GO/2116257_GI/'



#---------------------------100 nm full dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_100nm.csv")  
Obs.A1 <- Obs.A1[1:7,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","Clung")] #100nm only three organs
folder = 'plots/100nm_nsc_GI/'

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

#--------------draw the sensitive line one by one--------
Sens.init <- sensFun(func = pred.mouse,parms = params.init)

create_sensitivity_plots(Sens.init, names(params.init), paste0(folder,"init_sens/"))

# Identify the index of the 'var' column

var_index <- which(names(Sens.init) == "var")

# Exclude the 'var' column and apply the absolute function to make all values positive
Sens.init_abs <- Sens.init
Sens.init_abs[,-var_index] <- abs(Sens.init[,-var_index])

df_Sens.init_abs=summary(Sens.init_abs)
ranked_df_Sens_abs.init <-df_Sens.init_abs[order(df_Sens.init_abs$Mean, decreasing = TRUE), ]

graphics.off()
png(paste0(paste0(folder,"init_sens/"),"ranked_df_Sens_abs_plot_init.png"),width=3500,height=2800,res = 300)
plot(ranked_df_Sens_abs.init)
dev.off()

select_par = rownames(ranked_df_Sens_abs.init [ranked_df_Sens_abs.init $Mean>0.5,])
select_par

# Print selected parameters
print(select_par)

params2fit = params.init[select_par]

Fit.Result.A1<- modFit(f=MCcost, p=params2fit, obs=Obs.A1, method ="Nelder-Mead", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"
res.A1=MCcost(Fit.Result.A1$par, obs=Obs.A1)$residuals$res     ## Check the residual for each time points
sum(res.A1^2)      

Fitted_output.A1 = pred.mouse(par=Fit.Result.A1$par)

# Create individual plots
plot_liver.a1   <- create_plot(Fitted_output.A1, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a1  <- create_plot(Fitted_output.A1, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a1    <- create_plot(Fitted_output.A1, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a1  <- create_plot(Fitted_output.A1, y = "CS", organ = "Spleen", title = "Spleen")

# Arrange the plots together
combined_mod_fit_plot.a1 <- grid.arrange(plot_liver.a1, plot_lung.a1, plot_spleen.a1, 
                                         ncol = 2, nrow = 2)

combined_mod_fit_plot.a1 <- grid.arrange(plot_liver.a1, plot_kidney.a1, plot_lung.a1, plot_spleen.a1, 
                                         ncol = 2, nrow = 2)

# Save the combined plot
ggsave(paste0(folder, "mod_fit/mod_fit_combined_v1.png"), combined_mod_fit_plot.a1, width = 14, height = 10)

params.mod.a1 = params.init
# Iterate over the names of the new list
for (name in names(Fit.Result.A1$par)) {
  # Replace the element in the original list with the corresponding element from the new list
  params.mod.a1[[name]] <- Fit.Result.A1$par[[name]]
}

write.csv(params.mod.a1, file = paste0(folder,'mod_fit/params_fitted_v1.csv'))


#----------use all parameters as fitting----------
Fit.Result.A0<- modFit(f=MCcost, p=params.init, obs=Obs.A1, method ="Nelder-Mead", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"
res.A0=MCcost(Fit.Result.A0$par, obs=Obs.A1)$residuals$res     ## Check the residual for each time points
sum(res.A0^2)      

Fitted_output.A0 = pred.mouse(par=Fit.Result.A0$par)

# Create individual plots
plot_liver.a0   <- create_plot(Fitted_output.A0, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a0  <- create_plot(Fitted_output.A0, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a0    <- create_plot(Fitted_output.A0, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a0  <- create_plot(Fitted_output.A0, y = "CS", organ = "Spleen", title = "Spleen")

# Arrange the plots together
combined_mod_fit_plot.a0 <- grid.arrange(plot_liver.a0, plot_lung.a0, plot_spleen.a0, 
                                         ncol = 2, nrow = 2)

combined_mod_fit_plot.a0 <- grid.arrange(plot_liver.a0, plot_kidney.a0, plot_lung.a0, plot_spleen.a0, 
                                         ncol = 2, nrow = 2)


# Save the combined plot
ggsave(paste0(folder, "mod_fit/mod_fit_combined_v1_0.png"), combined_mod_fit_plot.a0, width = 14, height = 10)


params.mod.a0 = params.init
# Iterate over the names of the new list
for (name in names(Fit.Result.A0$par)) {
  # Replace the element in the original list with the corresponding element from the new list
  params.mod.a0[[name]] <- Fit.Result.A0$par[[name]]
}

write.csv(params.mod.a0, file = paste0(folder,'mod_fit/params_fitted_v1_0.csv'))



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
                             Obs.A1, rmse_data_tot, rmse_organ,tstep)

rmse_data_tot = rmse_result$rmse_data_tot
rmse_organ = rmse_result$rmse_organ

# Print or save the data frame with R-squared values for further analysis

library(ggplot2)

#------------plotting rmse-----------------------
# Plotting R-squared values against fold for each parameter

# Plotting R-squared values against fold for each parameter
rmse_plot = ggplot(rmse_organ, aes(x = Fold, y = RMSE, color = Organ)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.6, ymax = 1), fill = "gray80", alpha = 0.1) +
  geom_point(size=1) +
  geom_line(size=1) +
  scale_x_log10() +
  facet_wrap(~Parameter, scales = "free",ncol=8) +
  labs(title = "RMSE Values vs. Fold for Each Parameter",
       x = "Fold",
       y = "RMSE",
       color = "Organ")+
  #coord_cartesian(ylim = c(0, 1)) +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border around the plot
  )
# Save the plot as a file
ggsave(paste0(folder,"init_sens_v2/","rmse_organ_plot.png"), plot = rmse_plot, width = 20, height = 10, units = "in")



# Plotting R-squared values against fold for each parameter

rmse_tot = ggplot(rmse_data_tot, aes(x = Fold, y = RMSE)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.6, ymax = 1), fill = "gray80", alpha = 0.3) +
  geom_point(size=1) +
  geom_line() +
  scale_y_log10() +
  facet_wrap(~Parameter, scales = "free", ncol = 8) +
  labs(title = "Total RMSE Values vs. Fold for Each Parameter",
       x = "Fold",
       y = "RMSE") +
  theme_minimal() +
  coord_cartesian(ylim = c(min(rmse_data_tot$RMSE), max(rmse_data_tot$RMSE))) +
  theme(
    axis.text.x = element_text(angle=90),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border around the plot
  )

# Save the plot as a file
ggsave(paste0(folder,"init_sens_v2/","rmse_tot_plot.png"), plot = rmse_tot, width = 20, height = 10, units = "in")

rmse_data_tot

#-----------stats of rmse------------------
result <- rmse_data_tot %>%
  group_by(Parameter) %>%
  reframe(min_value = min(RMSE),
            first_value = RMSE[Fold == 1],
          median_value = median(RMSE),
            corresponding_Fold = Fold[which.min(RMSE)])
result = distinct(result)
filtered_result = result %>%
  #filter(corresponding_Fold != 0.001)%>%filter(corresponding_Fold != 1000)%>%
  filter(min_value/median_value < 0.95)
filtered_result

non_covered_result = filtered_result%>%filter(corresponding_Fold == 1000)
# TODO: write a while loop

folds <- seq(from = 1, to = 1e6, length.out = 1000)


rmse_result2 = calculate_rmse(non_covered_result$Parameter,params.init,folds, pred.mouse, Obs.A1, rmse_data_tot, rmse_organ)

rmse_data_tot2 = rmse_result2$rmse_data_tot
rmse_organ2 = rmse_result2$rmse_organ

result2 <- rmse_data_tot2 %>%
  group_by(Parameter) %>%
  reframe(min_value = min(RMSE),
          first_value = RMSE[Fold == 1],
          median_value = median(RMSE),
          corresponding_Fold = Fold[which.min(RMSE)])
result2 = distinct(result2)
filtered_result = result2 %>%
  #filter(corresponding_Fold != 0.001)%>%filter(corresponding_Fold != 1000)%>%
  filter(min_value/median_value < 0.96)

filtered_result


#--------use the parameter found in rmse + sensitive parameters found before----
select_par_v2 = c(select_par, "Kurine")

params2fit_v2 = params.init[select_par_v2]

Fit.Result.A2<- modFit(f=MCcost, p=params2fit, obs=Obs.A1, method ="Port", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"
res.A2=MCcost(Fit.Result.A2$par, obs=Obs.A1)$residuals$res     ## Check the residual for each time points
sum(res.A2^2)      

Fitted_output.A2 = pred.mouse(par=Fit.Result.A2$par)

# Create individual plots
plot_liver.a2   <- create_plot(Fitted_output.A2, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a2  <- create_plot(Fitted_output.A2, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a2    <- create_plot(Fitted_output.A2, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a2  <- create_plot(Fitted_output.A2, y = "CS", organ = "Spleen", title = "Spleen")

# Arrange the plots together
combined_mod_fit_plot.a2 <- grid.arrange(plot_liver.a2, plot_lung.a2, plot_spleen.a2, 
                                         ncol = 2, nrow = 2)

combined_mod_fit_plot.a2 <- grid.arrange(plot_liver.a2, plot_kidney.a2, plot_lung.a2, plot_spleen.a2, 
                                         ncol = 2, nrow = 2)

# Save the combined plot
ggsave(paste0(folder, "mod_fit/mod_fit_combined_v2.png"), combined_mod_fit_plot.a2, width = 14, height = 10)

params.mod = params.init
# Iterate over the names of the new list
for (name in names(Fit.Result.A2$par)) {
  # Replace the element in the original list with the corresponding element from the new list
  params.mod[[name]] <- Fit.Result.A2$par[[name]]
}

write.csv(params.mod, file = paste0(folder,'mod_fit/params_fitted_v2.csv'))


#-----------------------old: replace the initial value------
params.init2 = params.init
for (param in unique(filtered_result$Parameter)){
  cat("parameter:",param,"'\n")
  
  # Find the corresponding row in filtered_result for the current parameter
  
  corresponding_row <- filtered_result[filtered_result$Parameter == param, ]
  
  # Calculate the updated value for the parameter in params.init2
  updated_value <- log(exp(params.init[param]) * corresponding_row$corresponding_Fold)
  
  # Update the value in params.init2
  params.init2[param] <- updated_value
}

params.init2




#----------old: Draw the sensitivity line--------


# Initialize a list to store the plots
plot_list <- list()

parameters__ssens = params.init
R = pred.mouse(params.init)
folds <- seq(from = 10e-2, to = 10e2, length.out = 10)
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
  for (i in 1:length(parameters__ssens)) {
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

# Filter data for each organ
lsa_organs <- split(lsa, lsa$Organ)

# Function to plot and save each organ's plot
plot_and_save <- function(organ_data, organ_name) {
  gg <- ggplot(organ_data, aes(x = Time, y = predicted, color = Fold)) +
    geom_line() +
    geom_point(aes(y = initial_value), color = "blue", size = 0.5) +
    facet_wrap(parameter ~ Organ, scales = "free_y",ncol=5) +
    labs(title = paste0("LSA Data_",organ_name), x = "Time", y = "Concentration (ng/g)") +
    theme_minimal()
  
  ggsave( paste0(folder,"init_sens/",organ_name, ".png"), gg, width = 10, height = 6)
}

# Save plots for all organs
for (organ_name in names(lsa_organs)) {
  plot_and_save(lsa_organs[[organ_name]],organ_name)
}


# Filter data for each organ
lsa_paras <- split(lsa, lsa$parameter)

# Save plots for all organs
for (parameter_name in names(lsa_paras)) {
  plot_and_save(lsa_paras[[parameter_name]],parameter_name)
}



library(ggplot2)
# Filter data for the liver
lsa_liver <- subset(lsa, Organ == "Kidney")

lsa_liver <- subset(lsa, parameter == "K_max_Liver")

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




#--------------use the parameter value found in rmse-----

Sens.init2 <- sensFun(func = pred.mouse,parms = params.init2)
Sens.init_abs2 <- Sens.init2
#sensitivity value was calculated with the time range, therefore the value is 
# increase if the modeled time range is longer

# Identify the index of the 'var' column
var_index <- which(names(Sens.init2) == "var")

# Exclude the 'var' column and apply the absolute function to make all values positive
Sens.init_abs2[,-var_index] <- abs(Sens.init2[,-var_index])

df_Sens.init_abs2=summary(Sens.init_abs2)
ranked_df_Sens_abs.init2 <-df_Sens.init_abs2[order(df_Sens.init_abs2$Mean, decreasing = TRUE), ]

graphics.off()
png(paste0(paste0(folder,"init_sens_v2/"),"ranked_df_Sens_abs_plot_init.png"),width=3500,height=2800,res = 300)
plot(ranked_df_Sens_abs.init2)
dev.off()

select_par2 = rownames(ranked_df_Sens_abs.init2 [ranked_df_Sens_abs.init2 $Mean>0.5,])
select_par2 = c(select_par2, "K_max_Lung")

# Print selected parameters

params2fit2 = params.init2[select_par2]
#paramsnonfit = params.init2[setdiff(names(params.init2), select_par)]

Fit.Result.A2<- modFit(f=MCcost, p=params.init2, obs=Obs.A1, method ="Nelder-Mead", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"
res.A2=MCcost(Fit.Result.A2$par, obs=Obs.A1)$residuals$res     ## Check the residual for each time points
sum(res.A2^2)      

Fitted_output.A2 = pred.mouse(par=Fit.Result.A2$par)

# Create individual plots
plot_liver.a2   <- create_plot(Fitted_output.A2, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a2  <- create_plot(Fitted_output.A2, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a2    <- create_plot(Fitted_output.A2, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a2  <- create_plot(Fitted_output.A2, y = "CS", organ = "Spleen", title = "Spleen")

# Arrange the plots together
combined_mod_fit_plot.a2 <- grid.arrange(plot_liver.a2, plot_kidney.a2, plot_lung.a2, plot_spleen.a2, 
                                         ncol = 2, nrow = 2)

params2fit_v2 = Fit.Result.A2$par

# Iterate over the names of the new list
for (name in names(Fit.Result.A1$par)) {
  # Replace the element in the original list with the corresponding element from the new list
  params.mod[[name]] <- Fit.Result.A1$par[[name]]
}

write.csv(params.mod, file = paste0(folder,'mod_fit_v2/params_fitted.csv'))


