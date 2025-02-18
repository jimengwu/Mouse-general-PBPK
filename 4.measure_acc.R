
# ---------------------------------------------------------
# Author: Jimeng Wu
# Email: jimeng.wu@empa.ch
# Date: 2025-02-18
# Description: this function was used to calculate the prediction accuracy of the MCMC model,
# the pred-to-obs ratio was calculated and plotted was plotted, also the observation vs prediction plot was generated, 
# the adjusted R-squared and NRMSE was annotated on the plot
# ---------------------------------------------------------


library(ggplot2)     # Needed for plot
library(patchwork)
library(tidyverse)
library(Metrics)
library(dplyr)

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



# Function to calculate ratio for each column
calculate_ratio <- function(data, obs_data) {
  sapply(data[, 3:ncol(data)], function(col) col / obs_data)
}


# Define upper and lower limits for highlighting
ratio_upper_limit <- 2
ratio_lower_limit <- 0.5

# Function to calculate adjusted R-squared and NRMSE
calc_metrics <- function(data) {
  print(data)
  model <- lm(obs ~ pred, data = data)
  adj_r_squared <- summary(model)$adj.r.squared
  nrmse_value <- rmse(data$obs, data$pred) / mean(data$obs)
  
  return(data.frame(adj_r_squared = adj_r_squared, NRMSE = nrmse_value))
}


# np.name = "GO: Study1_20nm_20mg/kg"
#----------------1. calculate and plot the pred-to-obs ratio for every study----------------
nrmse_df <- data.frame(Iteration = character(), NRMSE = numeric(), stringsAsFactors = FALSE)
adjR_df <- data.frame(Iteration = character(), adjR = numeric(), stringsAsFactors = FALSE)
for (np.name in ls_np_name) {
  df_long = gen_obs_pred_data(np.name)
  df_long <- df_long %>% filter(obs != 0)
  
  folder = read_observation_data(np.name)$folder
  adj.R = summary(lm(obs~pred,data=df_long))$adj.r.squared

  NRMSE = rmse(df_long$obs,df_long$pred)/mean(df_long$obs)
  
  nrmse_df <- rbind(nrmse_df, data.frame(Iteration = np.name, NRMSE_value = NRMSE))
  
  adjR_df  <- rbind(adjR_df, data.frame(Iteration = np.name, adjR_value = adj.R))
  
  obs_pred_plot = ggplot(df_long, aes(x = obs, y = pred, shape = Variable)) +
    geom_point(size = 3) +
    geom_abline(slope = 1, linetype = "dashed") +  # Add dashed line with slope 1
    scale_shape_manual(values = c(1,2,3,4)) +  # Different shapes for each variable
    labs(
         x = "Observed NP Concentration in Organ (ng/g)",
         y = "Predicted NP Concentration in Organ (ng/g)") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "transparent"),  # Set background to transparent
          panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.ticks = element_line(),
          axis.ticks.length = unit(0.2,"cm"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          legend.position = c(0.8, 0.2),  # Adjust legend position (x, y)
          legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
          legend.key = element_blank(),  # Remove legend key
          legend.spacing.x = unit(0, "cm")  # Adjust horizontal spacing of legend
    )+
    scale_x_log10(limits = c(1e-1, 1e6))+  # Set x-axis to log scale
    scale_y_log10(limits = c(1e-1, 1e6)) +   # Set y-axis to log scale
    annotate("text", x = 0.1, y = 200000, 
             label = paste0("Adj.R^2 = ", round(adj.R, digits = 2),
                           "\nNRMSE = ", round(NRMSE, digits = 2)),
             size = 4.5, color = "black", hjust = 0, vjust = 0, fontface = "bold")
  
  # Save the combined plot
  #ggsave(paste0(folder, "obs_pred_combined.png"), obs_pred_plot, width = 8, height = 6)
  


  df_long$Ratio <- df_long$pred / df_long$obs
  
  
  
  # Plotting
  p_combined_ratio_highlighted = 
    ggplot(df_long, aes(x = pred, y = Ratio, shape = Variable)) +
    geom_point(size = 3) +
    labs(
      x = "Predicted value for Organ (ng/g)",
      y = "Ratio of Prediction-to-Observation") +
    theme_minimal() +
    geom_hline(yintercept = c(0.5, 2), linetype = "dashed", color = "red") +  # Add dashed lines
    scale_shape_manual(values = c(1,2,3,4))  +  # Define different shapes for each organ
    scale_y_continuous(trans = "log10", limits = c(0.01, 100)) + 
    scale_x_continuous(trans = "log10") +
    theme(panel.background = element_rect(fill = "transparent"),  # Set background to transparent
          panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.ticks = element_line(),
          axis.ticks.length = unit(0.2,"cm"),
          axis.text = element_text(size = 12),  # Increase text size for axis labels
          axis.title = element_text(size = 14, face = "bold"),  # Increase text size and bold axis titles
          legend.text = element_text(size = 10), # Increase text size for legend label
          legend.title = element_blank(),
          legend.position = c(0.8, 0.2),  # Adjust legend position (x, y)
          legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
          legend.key = element_blank(),  # Remove legend key
          legend.spacing.x = unit(0, "cm")  # Adjust horizontal spacing of legend
    )  # Increase thickness of axis lines)  # Remove minor grid lines
  
  
  # Save the combined plot
  #ggsave(paste0(folder, "prediction_acc_combined.png"), p_combined_ratio_highlighted, width = 8, height = 6)
  
  #write.csv(combined_ratio_data, file= paste0(folder,"pred_obs_ratio.csv"),row.names = FALSE)
  
  rm(df_obs_pred)
  
}





#-----------------2. calculate all together for total accuracy---------------

ls_np_name = c("Au: Study1_12nm_0.85mg/kg","Au: Study1_23nm_0.85mg/kg","Au: Study1_100nm_0.85mg/kg",
               "Au: Study2_34.6nm_3mg/kg","Au: Study2_55.5nm_3mg/kg","Au: Study2_77.1nm_3mg/kg",
               "Au: Study2_82.6nm_3mg/kg","Au: Study3_27.6nm_4.26mg/kg","Au: Study3_27.6nm_0.85mg/kg",
               "Si: Study1_20nm_10mg/kg","Si: Study1_80nm_10mg/kg","GO: Study1_20nm_20mg/kg",
               "GO: Study2_243nm_1mg/kg", "GO: Study2_914nm_1mg/kg_w/o_CS",
               "TiO2: Study1_385nm_10mg/kg","TiO2: Study2_220nm_60mg/kg",
               "FeO: Study1_29nm_5mg/kg","FeO: Study2_41nm_4mg/kg")



organ_study_df <- data.frame(Iteration = character(), 
                       adj_r_squared_Liver = numeric(),
                       adj_r_squared_Lung=numeric(),
                       adj_r_squared_Spleen=numeric(),
                       adj_r_squared_Kidney=numeric(),
                       NRMSE_Liver= numeric(),
                       NRMSE_Lung = numeric(),
                       NRMSE_Spleen=numeric(),
                       NRMSE_Kidney=numeric())

df_long_list <- list()

np.name ="Au: Study3_27.6nm_0.85mg/kg"



for (np.name in ls_np_name) {
  
  df_long = gen_obs_pred_data(np.name)
  df_long <- df_long %>% filter(obs != 0)
  
  
  adj.R = summary(lm(obs~pred,data=df_long))$adj.r.squared
  NRMSE = rmse(df_long$obs,df_long$pred)/mean(df_long$obs)
  

  # Assuming df_long has columns: obs, pred, and variable (to group by)
  
  # Group by the variable column and apply the calc_metrics function, variable columns
  # represents organ
  metrics_by_variable <- df_long %>%
    group_by(Variable) %>%
    do(calc_metrics(.))
  
  
  # Reshape metrics_by_variable into a single row dataset
  metrics_single_row <- metrics_by_variable %>%
    pivot_wider(names_from = Variable, values_from = c(adj_r_squared, NRMSE), names_sep = "_")
  metrics_single_row$Iteration = np.name
  
  organ_study_df <- bind_rows(organ_study_df, metrics_single_row)
  
  # Rename columns for clarity 
  colnames(metrics_single_row) <- gsub("adj_r_squared_", "adjR_", colnames(metrics_single_row))

  df_long$Ratio <- df_long$pred / df_long$obs
  df_long$ID <- np.name
  df_long_list <- append(df_long_list, list(df_long))
  

  rm(df_obs_pred)
  
}


folder_tot = "plots/"
#write.csv(nrmse_df,file=paste0(folder_tot,"nrmse_tot.csv"),row.names = FALSE)
#write.csv(adjR_df, file= paste0(folder_tot,"adjR_tot.csv"),row.names = FALSE)
#write.csv(organ_study_df, file= paste0(folder_tot,"acc_study_organ.csv"),row.names = FALSE)


# for all the study
df_tot = do.call(rbind, df_long_list)
adj.R_tot = summary(lm(obs~pred,data=df_tot))$adj.r.squared
NRMSE_tot = rmse(df_tot$obs,df_tot$pred)/mean(df_tot$obs)
adj.R_tot
NRMSE_tot
library(ggtext)
obs_pred_tot_plot <- ggplot(df_tot, aes(x = obs, y = pred, shape = Variable,color = ifelse(Time < 24, "Time < 24h", "Time >= 24h"))) +
  geom_point(size = 3) +
  geom_abline(slope = 1, linetype = "dashed") +  # Add dashed line with slope 1
  scale_shape_manual(values = c(1,2,3,4)) +  # Different shapes for each variable
  labs(
    x = "Observed NP Concentration in Organ (ng/g)",
    y = "Predicted NP Concentration in Organ (ng/g)") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent"),  # Set background to transparent
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.ticks = element_line(),
        axis.ticks.length = unit(0.2,"cm"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.2),  # Adjust legend position (x, y)
        legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
        legend.key = element_blank(),  # Remove legend key
        legend.spacing.x = unit(0, "cm")  # Adjust horizontal spacing of legend
  )+
  scale_x_log10(limits = c(1e-1, 1e6))+  # Set x-axis to log scale
  scale_y_log10(limits = c(1e-1, 1e6)) +   # Set y-axis to log scale
  annotate("text", x = 1, y = 150000, 
           label = bquote(atop("adj - R"^2 ~ ":" ~ .(round(adj.R_tot, 2)), 
                               " NRMSE" ~ ":" ~ .(round(NRMSE_tot, 2)))),
           size = 4.5,color = "black", hjust = 0, vjust = 0)+
  coord_fixed(ratio = 1)
obs_pred_tot_plot
#ggsave(paste0(folder_tot,"obs_pred_combined_time.png"), obs_pred_tot_plot, width = 6, height = 6) 

obs_pred_tot_plot <- ggplot(df_tot, aes(x = obs, y = pred, shape = Variable)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, linetype = "dashed") +  # Add dashed line with slope 1
  scale_shape_manual(values = c(1,2,3,4)) +  # Different shapes for each variable
  labs(
    x = "Observed NP Concentration in Organ (ng/g)",
    y = "Predicted NP Concentration in Organ (ng/g)") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent"),  # Set background to transparent
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.ticks = element_line(),
        axis.ticks.length = unit(0.2,"cm"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.2),  # Adjust legend position (x, y)
        legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
        legend.key = element_blank(),  # Remove legend key
        legend.spacing.x = unit(0, "cm")  # Adjust horizontal spacing of legend
  )+
  scale_x_log10(limits = c(1e-1, 1e6))+  # Set x-axis to log scale
  scale_y_log10(limits = c(1e-1, 1e6)) +   # Set y-axis to log scale
  annotate("text", x = 1, y = 150000, 
           label = bquote(atop("adj - R"^2 ~ ":" ~ .(round(adj.R_tot, 2)), 
                               " NRMSE" ~ ":" ~ .(round(NRMSE_tot, 2)))),
           size = 4.5,color = "black", hjust = 0, vjust = 0)+
  coord_fixed(ratio = 1)
obs_pred_tot_plot
#ggsave(paste0(folder_tot,"obs_pred_combined.png"), obs_pred_tot_plot, width = 6, height = 6) 


# Calculate the percentage
f2_error <- (sum(df_tot$Ratio > 0.5 & df_tot$Ratio < 2) / nrow(df_tot)) * 100
f2_error
# Calculate the percentage
f3_error <- (sum(df_tot$Ratio > 1/3 & df_tot$Ratio < 3) / nrow(df_tot)) * 100
f3_error

p_ratio_tot <- ggplot(df_tot, aes(x = pred, y = Ratio, shape = Variable)) +
  geom_point(size = 3) +
  labs(
    x = "Predicted value for Organ (ng/g)",
    y = "Ratio of Prediction-to-Observation") +
  theme_minimal() +
  geom_hline(yintercept = c(0.5, 2), linetype = "dashed", color = "red") +  # Add dashed lines
  scale_shape_manual(values = c(1,2,3,4))  +  # Define different shapes for each organ
  scale_y_log10(limits = c(0.01, 100)) + 
  scale_x_log10(limits = c(1e-1,1e6)) +
  annotate("text", x = 1, y = 20, 
           label = 
             #bquote(bold(atop("2f_error : 92.1%",
              #                      "3f_error : 97.3%")))
           bquote(atop("2f_error :" ~ .(sprintf("%.2f",f2_error))~"%", 
                       "3f_error :" ~ .(sprintf("%.2f",f3_error))~"%")),
           size = 4.5, color = "black", hjust = 0, vjust = 0, fontface = "bold") +
  theme(panel.background = element_rect(fill = "transparent"),  # Set background to transparent
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.ticks = element_line(),
        axis.ticks.length = unit(0.2,"cm"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.8),  # Adjust legend position (x, y)
        legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
        legend.key = element_blank(),  # Remove legend key
        legend.spacing.x = unit(0, "cm")  # Adjust horizontal spacing of legend
  )+ 
  coord_fixed(ratio = 7/4)


p_ratio_tot

#ggsave(paste0(folder_tot,"prediction_acc_combined.png"), p_ratio_tot, dpi=900, width =6 , height = 6) 


p_ratio_tot_highlighted <- ggplot(df_tot, aes(x = pred, y = Ratio, shape = Variable,color = ifelse(Time < 24, "Time < 24h", "Time >= 24h"))) +
  geom_point(size = 3) +
  labs(
    x = "Predicted value for Organ (ng/g)",
    y = "Ratio of Prediction-to-Observation") +
  theme_minimal() +
  geom_hline(yintercept = c(0.5, 2), linetype = "dashed", color = "red") +  # Add dashed lines
  scale_shape_manual(values = c(1,2,3,4))  +  # Define different shapes for each organ
  scale_y_log10(limits = c(0.01, 100)) + 
  scale_x_log10(limits = c(1e-1,1e6)) +
  annotate("text", x = 1, y = 20, 
           label = bquote(atop("2f_error :" ~ .(round(f2_error, 1))~"%", 
                               "3f_error :" ~ .(round(f3_error, 1))~"%")),
           size = 4.5, color = "black", hjust = 0, vjust = 0, fontface = "bold") +
  theme(panel.background = element_rect(fill = "transparent"),  # Set background to transparent
      panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.ticks = element_line(),
      axis.ticks.length = unit(0.2,"cm"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      legend.text = element_text(size=10),
      legend.title = element_blank(),
      legend.position = c(0.85, 0.8),  # Adjust legend position (x, y)
      legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
      legend.key = element_blank(),  # Remove legend key
      legend.spacing.x = unit(0, "cm")  # Adjust horizontal spacing of legend
)+ 
  coord_fixed(ratio = 7/4)


p_ratio_tot_highlighted

#ggsave(paste0(folder_tot,"prediction_acc_combined_time.png"), p_ratio_tot_highlighted, dpi=900, width =6 , height = 6) 
  

# calculate the accuracy based on organs
adj.R_tot_Liver = summary(lm(obs~pred,data= subset(df_tot, Variable=='Liver')))$adj.r.squared
adj.R_tot_Kidney = summary(lm(obs~pred,data= subset(df_tot, Variable=='Kidney')))$adj.r.squared
adj.R_tot_Spleen = summary(lm(obs~pred,data= subset(df_tot, Variable=='Spleen')))$adj.r.squared
adj.R_tot_Lung = summary(lm(obs~pred,data= subset(df_tot, Variable=='Lung')))$adj.r.squared

NRMSE_tot_Liver =  rmse(subset(df_tot, Variable=='Liver')$obs,subset(df_tot, Variable=='Liver')$pred)/mean(subset(df_tot, Variable=='Liver')$obs)
NRMSE_tot_Kidney = rmse(subset(df_tot, Variable=='Kidney')$obs,subset(df_tot, Variable=='Kidney')$pred)/mean(subset(df_tot, Variable=='Kidney')$obs)
NRMSE_tot_Spleen = rmse(subset(df_tot, Variable=='Spleen')$obs,subset(df_tot, Variable=='Spleen')$pred)/mean(subset(df_tot, Variable=='Spleen')$obs)
NRMSE_tot_Lung = rmse(subset(df_tot, Variable=='Lung')$obs,subset(df_tot, Variable=='Lung')$pred)/mean(subset(df_tot, Variable=='Lung')$obs)

# Create a data frame to save the results by different organ
organ_results <- data.frame(
  Organ = c('Liver', 'Kidney', 'Spleen', 'Lung'),
  Adjusted_R_squared = c(adj.R_tot_Liver, adj.R_tot_Kidney, adj.R_tot_Spleen, adj.R_tot_Lung),
  NRMSE = c(NRMSE_tot_Liver, NRMSE_tot_Kidney, NRMSE_tot_Spleen, NRMSE_tot_Lung)
)
#write.csv(organ_results, file= paste0(folder,"organ_tot.csv"),row.names = FALSE)
