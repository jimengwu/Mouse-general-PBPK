library(ggplot2)     # Needed for plot
library(patchwork)
Obs.A1 <- read.csv(file ="dataset/tk/mouse/TiO2/1-TiO2.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/TiO2/1_GI/'

source("helper_functions.R")
source("Mouse_PBPK.R")
pred.mouse <- function(pars,tstep=1) {
  
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
#---------------------------13 nm full dataset--------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_13nm.csv")  
# LONG TERM (WITHIN 24h)
Obs.A1 <- Obs.A1[1:7,]
PDOSE = 0.85
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/13nm_nsc_GI/'


##------------------GO: Study1: 20 nm high dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/GO/21162527_125I-NGS-PEG 10-30nm.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/GO/2116257_GI/'

##------------------TiO2: Study1: 20 nm high dose-------------------------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/TiO2/1-TiO2.csv") 
Obs.A1 <- Obs.A1[1:5,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/TiO2/1_GI/'

#----------------------ZnO-----------
Obs.A1 <- read.csv(file ="dataset/tk/mouse/ZnO/1-ZnO.csv") 
Obs.A1 <- Obs.A1[1:3,]
PDOSE = Obs.A1$Dose.mg.kg.[1]
Obs.A1 <- Obs.A1[c("Time","CL","CS","CK","Clung")]
folder = 'plots/ZnO/1_GI/'

mod <- mcode ("mouse_PBPK", mousePBPK.code)

set.seed(5)


Mouse.MCMC        <- readRDS(file = paste0(folder,"MCMC/mouse.MCMC.rds"))[[1]] # only first
theta.MCMC <- readRDS(paste0(folder,'MCMC/theta.rds'))
#Mouse.MCMC <- MCMC[[1]]
tstep = min(0.5,min(Obs.A1$Time))

Newtime.r   = pred.mouse(theta.MCMC,tstep)$Time  


# this is the new time variable, now it has been changed to sample per day.
nrwo.r = length (Newtime.r)
outputlength = nrow(Mouse.MCMC$pars)

MC.mouse.a.CL    = matrix(nrow = nrwo.r, ncol = outputlength/100)
MC.mouse.a.CK    = matrix(nrow = nrwo.r, ncol = outputlength/100)
MC.mouse.a.Clung = matrix(nrow = nrwo.r, ncol = outputlength/100)
MC.mouse.a.CS    = matrix(nrow = nrwo.r, ncol = outputlength/100)

for(i in 1:(outputlength/100)){
  
  j = i *100  # sample parameter set once every ten sets, 
  #so you will have 5000 sets from 50000 total sets
  
  pars.mouse             = Mouse.MCMC$pars    [j,]     
  MCdata               = pred.mouse (pars.mouse,tstep)
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

# Combine the datasets
combined_data <- rbind(
  cbind(Organs = "CL", MC.mouse.CL.plot),
  cbind(Organs = "CK", MC.mouse.CK.plot),
  cbind(Organs = "Clung", MC.mouse.Clung.plot),
  cbind(Organs = "CS", MC.mouse.CS.plot)
)

# Function to calculate ratio for each column
calculate_ratio <- function(data, obs_data) {
  sapply(data[, 3:ncol(data)], function(col) col / obs_data)
}

filtered_data <- combined_data[combined_data$Time %in% Obs.A1$Time, ]


# Calculate ratio for CL
ratio_CL <- calculate_ratio(subset(filtered_data, Organs == "CL"), as.numeric(as.character(Obs.A1$CL)))

# Calculate ratio for CK
ratio_CK <- calculate_ratio(subset(filtered_data, Organs == "CK"), as.numeric(as.character(Obs.A1$CK)))

# Calculate ratio for CS
ratio_CS <- calculate_ratio(subset(filtered_data, Organs == "CS"), as.numeric(as.character(Obs.A1$CS)))

# Calculate ratio for Clung
ratio_Clung <- calculate_ratio(subset(filtered_data, Organs == "Clung"), as.numeric(as.character(Obs.A1$Clung)))

# Define upper and lower limits for highlighting
ratio_upper_limit <- 2
ratio_lower_limit <- 0.5

# Create a data frame for plotting
plot_data2 <- data.frame(Time = Obs.A1$Time, 
                        Ratio = ratio_Clung[,"median_est"],
                        Pred = subset(filtered_data, Organs == "Clung")[,"median_est"],
                        Obs = Obs.A1$Clung,
                        Organ = "Lung")
plot_data3 <- data.frame(Time = Obs.A1$Time, 
                        Ratio = ratio_CK[,"median_est"],
                        Pred = subset(filtered_data, Organs == "CK")[,"median_est"],
                        Obs = Obs.A1$CK,
                        Organ = "Kidney")
plot_data4 <- data.frame(Time = Obs.A1$Time, 
                         Ratio = ratio_CS[,"median_est"],
                         Pred = subset(filtered_data, Organs == "CS")[,"median_est"],
                         Obs = Obs.A1$CS,
                         Organ = "Spleen")
plot_data5 <- data.frame(Time = Obs.A1$Time, 
                         Ratio = ratio_CL[,"median_est"],
                         Pred = subset(filtered_data, Organs == "CL")[,"median_est"],
                         Obs = Obs.A1$CL,
                         Organ = "Liver")

# Combine all plot_data data frames
combined_plot_data <- rbind(plot_data2, plot_data3, plot_data4, plot_data5)


# Plotting
p_combined_ratio_highlighted = ggplot(combined_plot_data, aes(x = Pred, y = Ratio, shape = Organ)) +
  geom_point() +
  labs(title = "13nm Au",
       x = "Predicted value for Organ (ug/g)",
       y = "Ratio of Prediction/Observation (-)") +
  theme_minimal() +
  geom_hline(yintercept = c(0.5, 2), linetype = "dashed", color = "red") +  # Add dashed lines
  scale_shape_manual(values = c(1, 2, 3, 4))  +  # Define different shapes for each organ
  scale_y_continuous(trans = "log10", limits = c(0.1, 10)) + 
  scale_x_continuous(trans = "log10") +
  theme(panel.background = element_rect(fill = "transparent"),  # Set background to transparent
                panel.grid.major = element_blank(),  # Remove major grid lines
                panel.grid.minor = element_blank(),  # Remove minor grid lines
                axis.text = element_text(size = 12),  # Increase text size for axis labels
                axis.title = element_text(size = 14, face = "bold"),  # Increase text size and bold axis titles
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Increase text size and bold plot title
                legend.title = element_text(size = 12, face = "bold"),  # Increase text size and bold legend title
                legend.text = element_text(size = 10) # Increase text size for legend labels
        )  # Increase thickness of axis lines)  # Remove minor grid lines

# Save the combined plot
ggsave(paste0(folder, "prediction_acc_combined.png"), p_combined_ratio_highlighted, width = 8, height = 6)

#------------------------------parameter distribution-----
library(ggplot2)
library(dplyr)
library(tidyr)

# List of folders where post_mean_sd_df files are located


# List of folders where CSV files are located
folders <- c(
  #"/Users/mmm/work/Mouse-PBPK/plots/13nm_nsc_GI/mc_sens/",
  "/Users/mmm/work/Mouse-PBPK/plots/31501470_2_nsc_GI/mc_sens",
  "/Users/mmm/work/Mouse-PBPK/plots/31501470_nsc_GI/mc_sens"
)  # Add more folders as needed


folders <- c("/Users/mmm/work/Mouse-PBPK/plots/4nm_nsc_GI/mc_sens/",
             "/Users/mmm/work/Mouse-PBPK/plots/13nm_nsc_GI/mc_sens/",
             "/Users/mmm/work/Mouse-PBPK/plots/TiO2/1_GI/mc_sens/",
             "/Users/mmm/work/Mouse-PBPK/plots/GO/2116257_GI/mc_sens/"
             #"/Users/mmm/work/Mouse-PBPK/plots/ZnO/1_GI/mc_sens"
             
)  # Add more folders as needed

# Read post_mean_sd_df files from each folder into a list o data frames
data_list <- lapply(folders, function(folder) {
  read.csv(file.path(folder, "post_paras_stats.csv"))
})

# Combine data frames into a single data frame
combined_post_mean_sd_df <- bind_rows(data_list, .id = "Folder")

# modify the fold id

combined_post_mean_sd_df <- combined_post_mean_sd_df %>%
  mutate(Folder = case_when(
    Folder == 1 ~ "4nm Au",
    Folder == 2 ~ "13nm Au",
    Folder == 3 ~ "385nm TiO2",
    Folder == 4 ~ "20nm GO",
    TRUE ~ as.character(Folder)
  ))


# Convert Folder column to factor
combined_post_mean_sd_df$Folder <- as.factor(combined_post_mean_sd_df$Folder)

# Create a list to store individual plots
plots <- list()

# Loop through each parameter to create individual box plots
for (param in unique(combined_post_mean_sd_df$Par)[1:29]) {
  # Filter data for the current parameter
  param_data <- combined_post_mean_sd_df[combined_post_mean_sd_df$Par == param, ]
  
  # Create the box plot for the current parameter
  plot <- ggplot(param_data, aes(x = Folder, y = Mean, fill = Folder)) +
    #geom_boxplot() +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) + # Add error bars for SD
  
    labs(title = paste("Box Plot for Parameter:", param)) +
    theme_minimal() +
    #geom_text(aes(label = paste("SD:", sprintf("%.3f", SD))), color = "blue", size = 3) +
    geom_text(aes(label = paste("Mean:", sprintf("%.3f", Mean))), color = "red", size = 6)+
    theme(axis.text = element_text(size = 14))
  
  # Add the plot to the list
  plots[[param]] <- plot
}

# Combine all the individual plots into one plot using patchwork
combined_plot <- wrap_plots(plots, ncol = 8)

# Print the combined plot
print(combined_plot)

ggsave(paste0("/Users/mmm/work/Mouse-PBPK/plots/paras","/post_density_distribution_plot_Au_TiO2_GO.png"), combined_plot, 
       width = 110, height = 50, units = "cm")

