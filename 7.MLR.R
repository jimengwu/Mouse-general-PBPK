# ---------------------------------------------------------
# Author: Jimeng Wu
# Email: jimeng.wu@empa.ch
# Date: 2025-02-18
# Description: multivariable linear regression analysis
# ---------------------------------------------------------

## loading R packages
library(magrittr)   # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)      # Needed for the pipe %>% operator
library(mrgsolve)   # Needed to run the main PBPK code
library(reshape)    # melt function to reshape the table
library(ggplot2)    # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(grid)       # for plotting the figure
library(lattice)    # for plotting the figure
library(grDevices)
library(gridExtra)
library(tidyr)
library(tibble)
library(olsrr)
library(Metrics)
library("readxl")
library(dplyr)
library(stringr)
source("helper_functions.R")
source("Mouse_PBPK.R")
source("dataset_info.R")


mlr_result_folder = "/Users/wuji/work/code/Mouse-general-PBPK/plots/mlr/"
mod <- mcode ("mouse_PBPK", mousePBPK.code)

ls_np_name = c("Au: Study1_12nm_0.85mg/kg","Au: Study1_23nm_0.85mg/kg","Au: Study1_100nm_0.85mg/kg",
               "Au: Study2_34.6nm_3mg/kg","Au: Study2_55.5nm_3mg/kg","Au: Study2_77.1nm_3mg/kg",
               "Au: Study2_82.6nm_3mg/kg","Au: Study3_27.6nm_4.26mg/kg","Au: Study3_27.6nm_0.85mg/kg",
               "Si: Study1_20nm_10mg/kg","Si: Study1_80nm_10mg/kg","GO: Study1_20nm_20mg/kg",
               "GO: Study2_243nm_1mg/kg","GO: Study2_914nm_1mg/kg_w/o_CS",
               "TiO2: Study1_385nm_10mg/kg","TiO2: Study2_220nm_60mg/kg",
               "FeO: Study1_29nm_5mg/kg","FeO: Study2_41nm_4mg/kg")
experiment_id <- c(
  "Experiment 5: Au_12nm_0.85mg/kg",
  "Experiment 6: Au_23nm_0.85mg/kg",
  "Experiment 7: Au_100nm_0.85mg/kg",
  "Experiment 8: Au_34.6nm_3mg/kg",
  "Experiment 9: Au_55.5nm_3mg/kg",
  "Experiment 10: Au_77.1nm_3mg/kg",
  "Experiment 11: Au_82.6nm_3mg/kg",
  "Experiment 12: Au_27.6nm_4.26mg/kg",
  "Experiment 13: Au_27.6nm_0.85mg/kg",
  "Experiment 3: SiO2_20nm_10mg/kg",
  "Experiment 4: SiO2_80nm_10mg/kg",
  "Experiment 14: GO_20nm_20mg/kg",
  "Experiment 15: GO_243nm_1mg/kg",
  "Experiment 16: GO_914nm_1mg/kg_w/o_CS",
  "Experiment 17: TiO2_385nm_10mg/kg",
  "Experiment 18: TiO2_220nm_60mg/kg",
  "Experiment 1: Iron_oxide_29nm_5mg/kg",
  "Experiment 2: Iron_oxide_41nm_4mg/kg"
)
for (i in 1:length(ls_np_name)) {
  
  np.name = ls_np_name[i]
  # Tell about progress
  cat('Processing study', i, 'of', length(ls_np_name),'\n')
  
  Obs.df = read_observation_data(np.name)$Obs.df
  print(max(Obs.df$Time))}

# ---------- 0. calculate the AUC and Cmax for each organ, and the delivery efficiency (all the PK parameters)------

AUC_result <- NULL  # Initialize the results dataframe
for (i in 1:length(ls_np_name)) {
  np.name = ls_np_name[i]
  # Tell about progress
  cat('Processing study', i, 'of', length(ls_np_name),'\n')
  
  Obs.df = read_observation_data(np.name)$Obs.df
  folder = read_observation_data(np.name)$folder
  pathway = read_observation_data(np.name)$pathway
  PDOSE = read_observation_data(np.name)$PDOSE
  
  if (np.name == "FeO: Study2_41nm_4mg/kg") {
    tstep <- 0.5/60 #2.5/60 in the time step
  } else if (np.name %in% c("GO: Study2_243nm_1mg/kg", "GO: Study2_914nm_1mg/kg_all",
                            "GO: Study2_914nm_1mg/kg_w/o_CS")) {
    # Assign the appropriate value for tstep
    tstep <- 1/60   # Fill in the appropriate value here 2/60 & 5/60 in the time points
  } else {
    tstep <- min(1, min(Obs.df$Time))
  }
  
  ## Loading human, rat, mouse, monkey MCMC data
  Mouse.MCMC        <- readRDS(file = paste0(folder,"MCMC/mouse.MCMC.rds"))
  
  ## loading the theta names
  theta             <- readRDS(file = paste0(folder,"MCMC/theta.rds"))
  theta.names       <- names(theta)
  which_sig         <- grep("sig", theta.names)
  tend              <- max(Obs.df$Time)
  
  ## Sensitivity analysis
  
  pars.mouse  = Mouse.MCMC[[1]]$bestpar
  names(pars.mouse) <- gsub("K_max","K_uptake", names(pars.mouse))
  
  R = Pred_auc(pars.mouse,PDOSE,tstep)
  
  # test for Trapezoidal Integration, it is fine to use the AUC in code
  #a = R$outdf
  #library(pracma)
  #trapz(a$Time, a$CLt_id_g)
  #trapz(a$Time, a$CLt_id)
  
  auc_cols <- grep("^AUC", names(R$outdf), value = TRUE)
  last_auc_row <- R$outdf[nrow(R$outdf), auc_cols, drop = FALSE]
  
  if (max(R$outdf$Time) >=168) {
    auc24_row <- R$outdf[which(R$outdf$Time == 24), auc_cols, drop = FALSE]%>%
      rename_with(~ str_replace_all(., "AUC", "AUC24"))
    
    auc168_row <- R$outdf[which(R$outdf$Time == 168), auc_cols, drop = FALSE]%>%
      rename_with(~ str_replace_all(., "AUC", "AUC168"))
  }else {auc24_row <- data.frame(
    AUC24_Kt_id_g = NA,
    AUC24_Lt_id_g = NA,
    AUC24_St_id_g = NA,
    AUC24_Lut_id_g = NA,
    AUC24_blood_id_g = NA
  )
  auc168_row <- data.frame(
         AUC168_Kt_id_g = NA,
         AUC168_Lt_id_g = NA,
         AUC168_St_id_g = NA,
         AUC168_Lut_id_g = NA,
         AUC168_blood_id_g = NA
       )
       }
  
  # Filter concentrations columns that start with "C"
  c_cols <- grep("^C", names(R$outdf), value = TRUE)
  
  # Keep only the relevant columns from the max_c_rows dataframe

  max_c_values <- sapply(c_cols, function(col) {
    max_pos <- which.max(R$outdf[[col]])
    c(max = R$outdf[[col]][max_pos], time = R$outdf$Time[max_pos])
  })
  

  max_c_row <- as.data.frame(t(max_c_values))

  max_c_row <- max_c_row %>%
    rownames_to_column(var = "id")%>%
    pivot_wider(names_from = id, values_from = !id)

  
  result <- cbind(last_auc_row, auc24_row,auc168_row,max_c_row)
  result$time_last = max(R$outdf$Time)
  result$id = np.name

  #colnames(AUC_result) = colnames(result)
  AUC_result <- rbind(AUC_result, result)
  
  
}


# caclualte the delivery efficiency based on the equation AUC/t_last
AUC_result$DE_L_id_g = AUC_result$AUC_Lt_id_g/AUC_result$time_last
AUC_result$DE_K_id_g = AUC_result$AUC_Kt_id_g/AUC_result$time_last
AUC_result$DE_S_id_g = AUC_result$AUC_St_id_g/AUC_result$time_last
AUC_result$DE_Lu_id_g = AUC_result$AUC_Lut_id_g/AUC_result$time_last

# caclualte the delivery efficiency based on the equation AUC/t_last
AUC_result$DE24_L_id_g = AUC_result$AUC24_Lt_id_g/24
AUC_result$DE24_K_id_g = AUC_result$AUC24_Kt_id_g/24
AUC_result$DE24_S_id_g = AUC_result$AUC24_St_id_g/24
AUC_result$DE24_Lu_id_g = AUC_result$AUC24_Lut_id_g/24

# caclualte the delivery efficiency based on the equation AUC/t_last
AUC_result$DE168_L_id_g = AUC_result$AUC168_Lt_id_g/168
AUC_result$DE168_K_id_g = AUC_result$AUC168_Kt_id_g/168
AUC_result$DE168_S_id_g = AUC_result$AUC168_St_id_g/168
AUC_result$DE168_Lu_id_g = AUC_result$AUC168_Lut_id_g/168


# todo: calculate the delivery efficiency based on the equation max/t_max
AUC_result$slope_CLt_id_g = AUC_result$max_CLt_id_g/AUC_result$time_CLt_id_g
AUC_result$slope_CKt_id_g = AUC_result$max_CKt_id_g/AUC_result$time_CKt_id_g
AUC_result$slope_CSt_id_g = AUC_result$max_CSt_id_g/AUC_result$time_CSt_id_g
AUC_result$slope_CLungt_id_g = AUC_result$max_CLungt_id_g/AUC_result$time_CLungt_id_g


#----combining the nanoparticle property information with the AUC result, and categorize the category values-----
dataset_info <- read_excel("dataset/tk/mouse/dataset_info.xlsx")

combined_data = merge(dataset_info, AUC_result, by = c("id"))

# Fill NA values using coalesce for hydrodynamic size and size
combined_data <- combined_data %>%
  mutate(
    Hydrodynamic_Size = coalesce(Hydrodynamic_Size, Size),
    Size = coalesce(Size, Hydrodynamic_Size)
  )

#Nanoparticles with a zeta potential between -10 and +10 mV are considered approximately neutral, 

combined_data <- combined_data %>%
  mutate(
    ZP.category = cut(zeta_potential, breaks = c(-Inf, -10, 10, Inf),
                      labels = c("negative", "neutral", "positive"),
                      include.lowest = TRUE)
  )

combined_data <- combined_data %>%
  mutate(
    ZP.category = ifelse(is.na(ZP.category),"no_info", as.character(ZP.category))
  )

#-------------------combine with the PBPK model parameter that fitted from MCMC model ------
# parameters has been logiramized so that avoid the predicted negative value
pars_data = read.csv(file=paste0("/Users/wuji/work/code/Mouse-general-PBPK/plots/paras/","pars_T_tot.csv"))
merged_data <- merge(combined_data, pars_data, by.x = "id",by.y = "Folder", all.x = TRUE) # merge with original value

# Apply log to only numeric columns
numeric_columns <- sapply(pars_data, is.numeric)  # Identify numeric columns
pars_data_numeric <- pars_data[, numeric_columns]  # Subset only numeric columns
log_pars_data_numeric <- log(pars_data_numeric)  # Apply log to numeric columns

# Combine with non-numeric columns if needed
pars_data_non_numeric <- pars_data[, !numeric_columns]  # Subset non-numeric columns
log_pars_data <- cbind(pars_data_non_numeric, log_pars_data_numeric)  # Combine both parts
log_pars_data <- log_pars_data %>%rename_with(~ paste0(., "_log"))
merged_data <- merge(merged_data, log_pars_data, by.x = "id",by.y = "pars_data_non_numeric_log", all.x = TRUE) # merge with log-transformed value


# calculate the log-transformed kinetic indicators

merged_data <- merged_data %>%
  mutate(across(starts_with("DE"), ~ log10(.), .names = "log.{col}"))

merged_data <- merged_data %>%
  mutate(across(starts_with("max_"), ~ log10(.), .names = "log.{col}"))

merged_data <- merged_data %>%
  mutate(across(starts_with("slope_"), ~ log10(.), .names = "log.{col}"))

merged_data <- merged_data %>%
  mutate(across(starts_with("slope_"), ~ ((.) - mean(.))/mean(.), .names = "norm.{col}"))

# for PBPK model parameters normalization
merged_data <- merged_data %>%
  mutate(across(all_of(gsub("max", "uptake", names(params.init))), ~ (. - mean(.)) / mean(.), .names = "norm_{col}"))

# for log transformed PBPK model parameters normalization
merged_data <- merged_data %>%
  mutate(across(all_of(paste0(gsub("max", "uptake", names(params.init)), "_log")), ~ ((.) - mean(.))/mean(.), .names = "norm.{col}"))

pars_pbpk = gsub("K_max","K_uptake", names(params.init))
#write.csv(merged_data,file=paste0("/Users/mmm/work/Mouse-PBPK/plots/mlr/","PK_results.csv"),row.names = FALSE)

# ---------- 1. normality distribution of the delivery efficiency check------



# Shapiro-Wilk normality test
# test for the PBPK model parameters
for (i in pars_pbpk){
  if (shapiro.test(merged_data[[i]])$p.value > 0.05)
  {
    print(i)
    print(shapiro.test(merged_data[[i]])$p.value)}
} # all parameters are normally distributed


# test for the log transformed PBPK model parameters
pars_log_pbpk = paste0(pars_pbpk, "_log")
for (i in pars_log_pbpk){
  
  if (shapiro.test(merged_data[[i]])$p.value > 0.05)
  {
    print(i)
    print(shapiro.test(merged_data[[i]])$p.value)}
} 
# not all parameters are normally distributed, A_cap_liver, k_release_liver, k_release_lung, 
# K_release_spleen, K_uptake_kidney, k_uptake_lung, P_kidney

# test for the normalized log transformed PBPK model parameters
colnames(merged_data)[137:165]
for (i in colnames(merged_data)[137:165]){
  if (shapiro.test(merged_data[[i]])$p.value > 0.05)
  {
    print(i)
    print(shapiro.test(merged_data[[i]])$p.value)}
} 
# same as unnormalized


# log transformed delivery efficiency 

ls_DE_log = c("log.DE_L_id_g","log.DE_K_id_g","log.DE_S_id_g","log.DE_Lu_id_g")
for (i in ls_DE_log){
  if (shapiro.test(merged_data[[i]])$p.value > 0.05)
  {
    print(i)
    print(shapiro.test(merged_data[[i]])$p.value)}
} 
# all log parameters are normally distributed, 
# except for log.DE_K_id_g, but it is also 0.053


# log transformed maximum concentration 


ls_max_log = c("log.max_CLt_id_g","log.max_CKt_id_g","log.max_CSt_id_g","log.max_CLungt_id_g")
for (i in ls_max_log){
  if (shapiro.test(merged_data[[i]])$p.value > 0.05)
  {
    print(i)
    print(shapiro.test(merged_data[[i]])$p.value)}
}
# use the log form of maximum concentration

# log transformed slope
ls_slope_log = c("log.slope_CLt_id_g","log.slope_CKt_id_g","log.slope_CSt_id_g","log.slope_CLungt_id_g")
for (i in ls_slope_log){
  if (shapiro.test(merged_data[[i]])$p.value > 0.05)
  {
    print(i)
    print(shapiro.test(merged_data[[i]])$p.value)}
} 



# draw the histogram for lig transformed delivery efficiency
hist(merged_data$log.DE24_L_id_g, probability=T) 
hist(merged_data$log.DE24_K_id_g, probability=T) 
hist(merged_data$log.DE24_S_id_g, probability=T) 
hist(merged_data$log.DE24_Lu_id_g, probability=T) 


# draw the histogram for maximum concentration
hist(merged_data$log.max_CLt_id_g, probability=T) 
hist(merged_data$log.max_CKt_id_g, probability=T) 
hist(merged_data$log.max_CSt_id_g, probability=T) 
hist(merged_data$log.max_CLungt_id_g, probability=T) 



hist(merged_data$log.slope_CLt_id_g, probability=T) 
hist(merged_data$log.slope_CKt_id_g, probability=T) 
hist(merged_data$log.slope_CSt_id_g, probability=T) 
hist(merged_data$log.slope_CLungt_id_g, probability=T) 




#-------------1.1 normality distribution of the PBPK model parameters----

for (i in pars_pbpk) {

  # Open a PNG device
  png(paste0(mlr_result_folder,"plots/PBPK/histogram_", i, ".png"))
  
  # Create histogram for the current column
  hist(merged_data[[i]], probability = TRUE, main = paste("Histogram of", i), xlab = i)
  
  # Close the PNG device
  dev.off()}


#------------2. Multivariable regression model-------
#---------------- 2.1 for PK parameters, C_max and Delivery efficiency -----

stats <- function(x){
  r2 <- summary(x)$r.squared
  r2a <- summary(x)$adj.r.squared
  f <- summary(x)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  error <- summary(x)$sigma
  y <- c(r2, r2a, f, p, error)
  return(y)
}

# Initialize vectors to store results
variable_names <- c()
adj_r_squared_full <- c()
aic_full <- c()
r_squared_full <- c()
p_full <- c()

adj_r_best <- c()
aic_best <- c()
r_squared_best <- c()
p_best <- c()
equation_best<-c()
equation_full <-c()
predictor_names <- c()
PK_pred_results <- data.frame()


ls_PK = c("max_CLt_id_g","max_CKt_id_g","max_CSt_id_g","max_CLungt_id_g",
          "log.max_CLt_id_g","log.max_CKt_id_g","log.max_CSt_id_g","log.max_CLungt_id_g",
          "log.DE_L_id_g","log.DE_K_id_g","log.DE_S_id_g","log.DE_Lu_id_g",
          "log.DE24_L_id_g","log.DE24_K_id_g","log.DE24_S_id_g","log.DE24_Lu_id_g",
          "log.DE168_L_id_g","log.DE168_K_id_g","log.DE168_S_id_g","log.DE168_Lu_id_g",
          "slope_CLt_id_g","slope_CKt_id_g","slope_CSt_id_g","slope_CLungt_id_g",
          "log.slope_CLt_id_g","log.slope_CKt_id_g","log.slope_CSt_id_g","log.slope_CLungt_id_g")

for (i in ls_PK){
  print(i)
  variable_names <- c(variable_names, i)
  
  # Remove rows with NaN in the column corresponding to 'i'
  merged_data_clean <- merged_data[!is.na(merged_data[[i]]), ]
  
  formula_full <- lm(as.formula(paste(i,"~ NP_core + Hydrodynamic_Size + coating + ZP.category + shape + Dose")), data = merged_data_clean)
  fit_all <- ols_step_all_possible(formula_full,max_order=6)
  best_fit <- fit_all$result %>% arrange(desc(adjr)) %>%  arrange(aic) %>% head(1)
  formula_best <- lm(as.formula(paste(i, "~", paste(unlist(strsplit(unlist(best_fit$predictors[[1]]), " ")), collapse = "+"))),
                     data=merged_data_clean)
  
  r_squared_full <- c(r_squared_full,tail(fit_all$result)[1,]$rsquare)
  adj_r_squared_full <- c(adj_r_squared_full,tail(fit_all$result)[1,]$adjr)
  aic_full<- c(aic_full,tail(fit_all$result)[1,]$aic)
  p_full <- c(p_full,stats(formula_full)[6])
  
  r_squared_best <- c(r_squared_best,best_fit$rsquare)
  adj_r_best <- c(adj_r_best, best_fit$adjr)
  aic_best <- c(aic_best, best_fit$aic)
  p_best <- c(p_best,stats(formula_best)[6])
  
  equation_full = c(equation_full,gene_mlr_equation(formula_full))
  equation_best = c(equation_best,gene_mlr_equation(formula_best))
  
  #--------- for the obs vs pred plot --------------
  predictor_names <- c(predictor_names,best_fit$predictors)
  predicted_full <- predict(formula_full, merged_data_clean)
  predicted_best <- predict(formula_best, merged_data_clean)
  
  # DataFrame for plotting
  PK_pred_i <- data.frame(
    Variable = i,
    case = merged_data_clean$id,
    Observed = merged_data_clean[[i]],
    Predicted_Full = predicted_full,
    Predicted_Best = predicted_best,
    adj_r_full = summary(lm(merged_data_clean[[i]]~predicted_full))$adj.r.squared,
    adj_r_best = summary(lm(merged_data_clean[[i]]~predicted_best))$adj.r.squared
  )

  # Combine the plot data
  PK_pred_results <- rbind(PK_pred_results, PK_pred_i)
}



# Create a dataframe with the results
full_PK_df <- data.frame(
  Variable = variable_names,
  Adj_R_Squared_full = adj_r_squared_full,
  aic_full = aic_full,
  r_squared_full = r_squared_full,
  p_full = p_full,
  equation_full = equation_full,
  Adj_R_Squared_best = adj_r_best,
  aic_best = aic_best,
  r_squared_best = r_squared_best,
  p_best = p_best,
  predictors = predictor_names,
  equation_best = equation_best
  
)

full_PK_df$predictors <- gsub(" ", ", ", full_PK_df$predictors)     # Replace spaces with " + "
full_PK_df$predictors <- paste("Y ~ f (", full_PK_df$predictors,")", sep = "")         # Add "Y ~" to the start of each value
unique_equation_best_PK <- full_PK_df[grep("^log", full_PK_df$Variable), c("Variable", "equation_best")]
#write.csv(PK_pred_results,file=paste0(mlr_result_folder,"MLR_PK_all.csv"),row.names = FALSE)

#write.csv(full_PK_df,file=paste0(mlr_result_folder,"MLR_PK.csv"),row.names = FALSE)

#------------------ 2.1.1 plotting --------------

#--------------------- 2.1.1.1 for log max concentration -----

log.max_data <- PK_pred_results %>% filter(startsWith(Variable, "log.max"))


log.max_data <- log.max_data %>%
  mutate(Variable = recode(Variable,
                           "log.max_CLt_id_g" = "Liver",
                           "log.max_CKt_id_g" = "Kidney",
                           "log.max_CSt_id_g" = "Spleen",
                           "log.max_CLungt_id_g" = "Lung"))

log.max_plot_full <- obs_pred_plot_func_mlr(log.max_data,"Observed","Predicted_Full","Variable","adj_r_full") + xlab("Observed log Maximum Concentration") +
  ylab("Predicted log Maximum Concentration")
log.max_plot_full
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"plots/organ_log_cmax_full.png"), plot = log.max_plot_full, width = 6, height = 6)

log.max_plot_best <- obs_pred_plot_func_mlr(log.max_data,"Observed","Predicted_Best","Variable","adj_r_best")+ xlab("Observed log Maximum Concentration") +
  ylab("Predicted log Maximum Concentration")
log.max_plot_best
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"plots/organ_log_cmax_best.png"), plot = log.max_plot_best, width = 6, height = 6)
#--------------------- 2.1.1.2 for log delivery efficiency 24 -----

log.DE24_data <- PK_pred_results %>% filter(startsWith(Variable, "log.DE24"))


log.DE24_data <- log.DE24_data %>%
  mutate(Variable = recode(Variable,
                           "log.DE24_L_id_g" = "Liver",
                           "log.DE24_K_id_g" = "Kidney",
                           "log.DE24_S_id_g" = "Spleen",
                           "log.DE24_Lu_id_g" = "Lung"))

log.DE24_plot_full <- obs_pred_plot_func_mlr(log.DE24_data,"Observed","Predicted_Full","Variable","adj_r_full")+ xlab("Observed log DE24") +
  ylab("Predicted log DE24")
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"plots/organ_log_DE24_full.png"), plot = log.DE24_plot_full, width = 6, height = 6)

log.DE24_plot_best <- obs_pred_plot_func_mlr(log.DE24_data,"Observed","Predicted_Best","Variable","adj_r_best")+ xlab("Observed log DE24") +
  ylab("Predicted log DE24")
log.DE24_plot_best
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"plots/organ_log_DE24_best.png"), plot = log.DE24_plot_best, width = 6, height = 6)

#--------------------- 2.1.1.3 for log delivery efficiency 168 -----

log.DE168_data <- PK_pred_results %>% filter(startsWith(Variable, "log.DE168"))


log.DE168_data <- log.DE168_data %>%
  mutate(Variable = recode(Variable,
                           "log.DE168_L_id_g" = "Liver",
                           "log.DE168_K_id_g" = "Kidney",
                           "log.DE168_S_id_g" = "Spleen",
                           "log.DE168_Lu_id_g" = "Lung"))

log.DE168_plot_full <- obs_pred_plot_func_mlr(log.DE168_data,"Observed","Predicted_Full","Variable","adj_r_full")+ xlab("Observed log DE168") +
  ylab("Predicted log DE168")
log.DE168_plot_full
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"plots/organ_log_DE168_full.png"), plot = log.DE168_plot_full, width = 6, height = 6)

log.DE168_plot_best <- obs_pred_plot_func_mlr(log.DE168_data,"Observed","Predicted_Best","Variable","adj_r_best")+ xlab("Observed log DE168") +
  ylab("Predicted log DE168")
log.DE168_plot_best
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"plots/organ_log_DE168_best.png"), plot = log.DE168_plot_best, width = 6, height = 6)

#--------------------- 2.1.1.4 for log rate of increase -----

log_slope_data <- PK_pred_results %>% filter(startsWith(Variable, "log.slope"))

log_slope_data <- log_slope_data %>%
  mutate(Variable = recode(Variable,
                     "log.slope_CLt_id_g" = "Liver",
                     "log.slope_CKt_id_g" = "Kidney",
                     "log.slope_CSt_id_g" = "Spleen",
                     "log.slope_CLungt_id_g" = "Lung"))


log_slope_plot_full <- obs_pred_plot_func_mlr(log_slope_data,"Observed","Predicted_Full","Variable","adj_r_full")+ xlab("Observed log Average Rate of Increase") +
  ylab("Predicted log Average Rate of Increase")
log_slope_plot_full
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"plots/organ_log_slope_full.png"), plot = log_slope_plot_full, width = 6, height = 6)

log_slope_plot_best <- obs_pred_plot_func_mlr(log_slope_data,"Observed","Predicted_Best","Variable","adj_r_best")+ xlab("Observed log Average Rate of Increase") +
  ylab("Predicted log Average Rate of Increase")
log_slope_plot_best
# Save the plot to a PNG file
ggsave(filename = paste0(mlr_result_folder,"plots/organ_log_slope_best.png"), plot = log_slope_plot_best, width = 6, height = 6)




#--------------------- 2.1.1.5 plot for all the kinetic indicators -----
# Combine data for 24 and 168 hours, with a time column
log.DE24_data <- log.DE24_data %>% mutate(time = "DE24")
log.DE168_data <- log.DE168_data %>% mutate(time = "DE168")
log.max_data <- log.max_data %>% mutate(time = "maximum concentration")
log_slope_data <- log_slope_data %>% mutate(time = "average rate of increase")

combined_KI_data <- bind_rows(log.DE24_data, log.DE168_data,log.max_data,log_slope_data)

# Call the plotting function
log_KI_plot <- obs_pred_plot_func_mlr(combined_KI_data, "Observed", "Predicted_Best",
                                      "Observed Values for Kinetic Indicators",
                                      "Predicted Values for Kinetic Indicators", "time")
log_KI_plot
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"plots/MLR_log_KI_best.png"), plot = log_KI_plot, width = 6, height = 6)


#---------------- 2.2 for PBPK model parameters -----
#------------------- 2.2.1 for log transformed PBPK model parameters----

# Initialize vectors to store results
variable_names <- c()
adj_r_squared_full <- c()
aic_full <- c()
p_full <- c()

adj_r_best <- c()
aic_best <- c()
predictor_names <- c()
formula_best <- c()
p_best <- c()

all_plot_PBPK_data_log <- data.frame()

# Replace "max" with "uptake" in the names of params.init
for (i in paste0(gsub("max", "uptake", names(params.init)), "_log")){
  print(i)
  
  # the variable that has been studied
  variable_names <- c(variable_names, i)
  
  fit_full <- lm(as.formula(paste(i,"~ NP_core + Hydrodynamic_Size + coating + ZP.category + shape + Dose")), data = merged_data)
  
  fit_all <- ols_step_all_possible(fit_full) # get every fitting combination
  best_fit <- fit_all$result %>% arrange(desc(adjr)) %>%  arrange(aic) %>% head 

  adj_r_squared_full <- c(adj_r_squared_full,tail(fit_all$result)[1,]$adjr) # r squared from the linear regression fitting process
  aic_full<- c(aic_full,tail(fit_all$result)[1,]$aic)
  
  adj_r_best <- c(adj_r_best, best_fit[1,]$adjr)
  predictor_names <- c(predictor_names,best_fit[1,]$predictors)
  aic_best <- c(aic_best, best_fit[1,]$aic)
  
  # calculate the predicted value from full equation or the best equation
  predicted_full <- predict(fit_full, merged_data)
  p_full <- c(p_full,stats(fit_full)[6])
  formula_best <- lm(as.formula(paste(i, "~", paste(unlist(strsplit(unlist(best_fit$predictors[[1]]), " ")), collapse = "+"))),data=merged_data)
  predicted_best <- predict(formula_best, merged_data)
  p_best <- c(p_best,stats(formula_best)[6])

  # DataFrame for plotting, the r squared caculated here is from the linear regression between observed and predicted value
  log_PBPK_pred_i <- data.frame(
    Variable = i,
    id = merged_data$id,
    Observed = merged_data[[i]],
    Predicted_Full = predicted_full,
    Predicted_Best = predicted_best,
    adj_r_full_pred = summary(lm(merged_data[[i]]~predicted_full))$adj.r.squared, # r squared from the linear regression between observed and predicted value
    adj_r_best_pred = summary(lm(merged_data[[i]]~predicted_best))$adj.r.squared,
    equation_full = gene_mlr_equation(fit_full),
    equation_best = gene_mlr_equation(formula_best)
  )


  # Combine the plot data，which has all the cases for each variable
  all_plot_PBPK_data_log <- rbind(all_plot_PBPK_data_log, log_PBPK_pred_i)
  
}

# Create a dataframe with the results, just the summary of each variable 
full_PBPK_df_log <- data.frame(
  Variable = variable_names,
  Adj_R_Squared_full = adj_r_squared_full,
  aic_full = aic_full,
  p_full = p_full,
  Adj_R_Squared_best = adj_r_best,
  aic_best = aic_best,
  p_best = p_best,
  predictors = predictor_names
)


# calculate the R squred for each study case
case_acc_PBPK_log = all_plot_PBPK_data_log%>%
  group_by(id) %>%
  summarise(
    r_squared = 1 - sum((Observed - Predicted_Full)^2) / sum((Observed - mean(Observed))^2)
  )

equation_best_PBPK = unique(all_plot_PBPK_data_log[, c("Variable", "equation_best")])


full_PBPK_df_log$predictors <- gsub(" ", ", ", full_PBPK_df_log$predictors)     # Replace spaces with " + "
full_PBPK_df_log$predictors <- paste("Y ~ f (", full_PBPK_df_log$predictors,")", sep = "")         # Add "Y ~" to the start of each value

#-------- 2.2.2.1 plotting ---------
# every PBPK parameter in every case, obs vs pred R squared calculation

acc_info <- bquote(
  atop(
    " adj -" ~ R^2 ~ ":"~ .(round(summary(lm(all_plot_PBPK_data_log$Observed ~ all_plot_PBPK_data_log$Predicted_Full))$adj.r.squared, 3)))
)
pbpk_log_obs_pred_plot = obs_pred_plot_func_mlr(data = all_plot_PBPK_data_log,x_col="Observed",y_col="Predicted_Full",annotation_text = acc_info)
pbpk_log_obs_pred_plot
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"plots/PBPK_log_full.png"), plot = pbpk_log_obs_pred_plot, width = 6, height = 6)


# Create the label using bquote
a_best <- bquote(
  atop(
    "  adj - R"^2 * " : " * .(round(summary(lm(all_plot_PBPK_data_log$Observed ~ 
                                                   all_plot_PBPK_data_log$Predicted_Best))$adj.r.squared, 2)))
)

pbpk_log_obs_pred_plot_best = obs_pred_plot_func_mlr(all_plot_PBPK_data_log,"Observed","Predicted_Best",
                                                     x_col_text = "Observed Values for PBPK Model Parameters",
                                                     y_col_text = "Predicted Values for PBPK Model Parameters",
                                                     annotation_text = a_best)
pbpk_log_obs_pred_plot_best
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"plots/PBPK_log_best.png"), plot = pbpk_log_obs_pred_plot_best, width = 6, height = 6)


#write.csv(full_PBPK_df_log,file=paste0(mlr_result_folder,"MLR_PBPK_log.csv"),row.names = FALSE)
sorted = full_PBPK_df_log%>%arrange(Adj_R_Squared_best)
# Plot
sorted_PBPK_mlr_acc_plot = ggplot(sorted, aes(x = reorder(Variable,Adj_R_Squared_best), y = Adj_R_Squared_best)) +
  geom_col(aes(fill = Adj_R_Squared_best < 0.6)) +  # Color bars based on the condition
  scale_fill_manual(values = c("gray", "black"), guide = "none") +  # Set colors: red for < 0.6, gray otherwise
  labs(x = "Variable", y = "adj.Rsquared Value") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +  # Keep y-axis labels horizontal
  coord_flip() + # Flip the coordinates for a horizontal bar plot
theme_bw()
sorted_PBPK_mlr_acc_plot
#ggsave(filename = paste0(mlr_result_folder,"MLR_PBPK_log.png"), 
#       plot = sorted_PBPK_mlr, width = 8, height = 4)

#---------calculate the rsaqure for each case----

# Function to calculate R-squared
calculate_r_squared <- function(obs, pred) {
  ss_res <- sum((obs - pred) ^ 2)
  ss_tot <- sum((obs - mean(obs)) ^ 2)
  r_squared <- 1 - (ss_res / ss_tot)
  return(r_squared)
}

# Group by `id` and calculate R-squared for each group
PBPK_per_id <- all_plot_PBPK_data_log %>%
  group_by(id) %>%
  summarise(R_squared = calculate_r_squared(Observed, Predicted_Best))


# plot the pbpk model parameter one by one
for (i in unique(all_plot_PBPK_data_log$Variable)){
  i = "K_release_Liver_log"
  data_i = all_plot_PBPK_data_log %>% filter(Variable == i)
  plot_i_mlr_full <-data_i %>% ggplot(aes(x = Observed, y = Predicted_Full)) +
    geom_point(color = 'red') +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'black') +
    ggtitle(paste("Observed vs Predicted (Full Model) for", i)) +
    xlab("Observed Values") +
    ylab("Predicted Values") +
    theme_minimal()+
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
          legend.spacing.x = unit(0, "cm")) + # Adjust horizontal spacing of legend
    xlim(min(data_i$Predicted_Full,data_i$Observed),max(data_i$Predicted_Full,data_i$Observed)) +
    ylim(min(data_i$Predicted_Full,data_i$Observed),max(data_i$Predicted_Full,data_i$Observed))+
    annotate("text", x = max(data_i$Predicted_Full,data_i$Observed)*0.6, y = max(data_i$Predicted_Full,data_i$Observed)*0.9, 
             label =              bquote(atop("adj - R"^2 ~ ":" ~ .(round(data_i$adj_r_full[1], digits = 2)))),
             size = 4.5, color = "black", hjust = 0, vjust = 0, fontface = "bold")

  pred = summary(lm(Observed~Predicted_Full,data=data_i))
  print(i)
  print(pred$adj.r.squared)
  # Save the plot to a PNG file-
  #ggsave(filename = paste0(mlr_result_folder,"plots/PBPK_log/",i,"_mlr.png"), plot = plot_i_mlr_full, width = 6, height = 6)
  
}

#------------------- 2.2.3 histogram plot for PBPK_log-----
combined_plots = list()
# Loop through unique variables
for (i in unique(all_plot_PBPK_data_log$Variable)) {
  
  # Filter the data for the current variable
  data_i <- all_plot_PBPK_data_log %>% filter(Variable == i)
  
  pred <- summary(lm(Observed ~ Predicted_Full, data = data_i))
  print(i)
  print(pred$adj.r.squared)
  min_axis_lim = round(min(data_i$Observed, data_i$Predicted_Full)) - 1
  max_axis_lim = round(max(data_i$Observed, data_i$Predicted_Full))
  # Create the histogram plot
  plot_i_histogram <- data_i %>%
    ggplot() +
    geom_histogram(aes(x = Observed, fill = "Observed"), 
                  bins = 30, alpha = 0.5, position = "identity") +
    geom_histogram(aes(x = Predicted_Full, fill = "Predicted"), 
                  bins = 30, alpha=0.5, position = "identity") +
    ggtitle(paste(i)) +
    xlab("Values") +
    ylab("Frequency") +
    scale_fill_manual(values = c("Observed" = "red", "Predicted" = "blue")) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "transparent"),  # Set background to transparent
          panel.border = element_rect(fill=NA, color="black", size=2, linetype="solid"),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.ticks = element_line(),
          axis.ticks.length = unit(0.2, "cm"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          legend.position = c(0.8, 0.8),  # Adjust legend position (x, y)
          legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
          legend.key = element_blank(),  # Remove legend key
          legend.spacing.x = unit(0, "cm")) + # Adjust horizontal spacing of legend
    xlim(min_axis_lim,max_axis_lim) +
    ylim(0, max(
      hist(data_i$Observed, plot = FALSE, breaks = 30)$counts,
      hist(data_i$Predicted_Full, plot = FALSE, breaks = 30)$counts
    )) +
    annotate("text", x = min_axis_lim + (max_axis_lim - min_axis_lim)*0.2, 
             y = max(
               hist(data_i$Observed, plot = FALSE, breaks = 30)$counts,
               hist(data_i$Predicted_Full, plot = FALSE, breaks = 30)$counts
             ) * 0.9, 
             label = paste0("Adj.R^2 = ", round(pred$adj.r.squared, digits = 3)),
             size = 4.5, color = "black", hjust = 0, vjust = 0, fontface = "bold")
  plot_i_histogram
  # Print the variable and adjusted R-squared

  combined_plots[[i]] <- plot_i_histogram
  
  # Save the plot to a PNG file
  #ggsave(filename = paste0(mlr_result_folder, "plots/PBPK_log/", i, "_histogram.png"), 
  #       plot = plot_i_histogram, width = 6, height = 6)
}
library(patchwork)
final_plot <- wrap_plots(plotlist = combined_plots, ncol = 6)
final_plot


#------------3. frequency of predictors in the best model for PK parameters, C_max and Delivery efficiency ------
#-----------3.1 for PBPK model parameters------------
# Flatten the list of lists
flattened_list <- unlist(strsplit(unlist(full_PBPK_df_log$predictors), " "))


flattened_list <- trimws(unlist(strsplit(gsub(".*\\((.*?)\\).*", "\\1", full_PBPK_df_log$predictors), ",")))
# Specify the numbers to count
numbers_to_count <- unique(flattened_list)

# Count the occurrences of each number
counts <- sapply(numbers_to_count, function(x) sum(flattened_list == x))

# Create a data frame for ggplot2
data <- data.frame(Number = as.factor(numbers_to_count), Frequency = counts)

# Generate the bar plot using ggplot2
freq_plot_PBPK <- ggplot(data, aes(x = reorder(Number, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "gray",color = "black") +
  xlab("Predictors") +
  ylab("Frequency") +
  theme_bw() + ylim(0, 25) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 2),  # Thicker border around the plot
    axis.text.x = element_text(size = 14, face = "bold"),  # Increase x-axis text size and make it bold
    axis.text.y = element_text(size = 14, face = "bold"),  # Increase y-axis text size and make it bold
    axis.title.x = element_text(size = 16, face = "bold"), # Increase x-axis title size
    axis.title.y = element_text(size = 16, face = "bold")
  )+ guides(color = 'none') +coord_flip()

freq_plot_PBPK
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"frequency_PBPK_log.png"), plot = freq_plot_PBPK, width = 6, height = 4)

#--------3.2 for kinetic indicators -------
# Flatten the list of lists
selected_PK_df = full_PK_df%>%
  filter(Variable %in% c("log.max_CLt_id_g","log.max_CKt_id_g","log.max_CSt_id_g","log.max_CLungt_id_g",
                         "log.DE24_L_id_g","log.DE24_K_id_g","log.DE24_S_id_g","log.DE24_Lu_id_g",
                         "log.DE168_L_id_g","log.DE168_K_id_g","log.DE168_S_id_g","log.DE168_Lu_id_g",
                         "log.slope_CLt_id_g","log.slope_CKt_id_g","log.slope_CSt_id_g","log.slope_CLungt_id_g"))
flattened_list_PK <- unlist(strsplit(unlist(selected_PK_df$predictors), " "))

flattened_list_PK <- trimws(unlist(strsplit(gsub(".*\\((.*?)\\).*", "\\1", selected_PK_df$predictors), ",")))

# Specify the numbers to count
numbers_to_count_PK <- unique(flattened_list_PK)

# Count the occurrences of each number
counts_PK <- sapply(numbers_to_count_PK, function(x) sum(flattened_list_PK == x))
# Create a data frame for ggplot2
data_PK <- data.frame(Number = as.factor(numbers_to_count_PK), Frequency = counts_PK)


# Generate the bar plot using ggplot2
freq_plot_PK <- ggplot(data_PK, aes(x = reorder(Number, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "gray",color = "black") +
  xlab("Predictors") +
  ylab("Frequency") +
  theme_bw() + ylim(0, 15) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 2),  # Thicker border around the plot
    axis.text.x = element_text(size = 14, face = "bold"),  # Increase x-axis text size and make it bold
    axis.text.y = element_text(size = 14, face = "bold"),  # Increase y-axis text size and make it bold
    axis.title.x = element_text(size = 16, face = "bold"), # Increase x-axis title size
    axis.title.y = element_text(size = 16, face = "bold")
  )+ guides(color = 'none') +coord_flip()



freq_plot_PK
# Save the plot to a PNG file
#ggsave(filename = paste0(mlr_result_folder,"frequency_PK.png"), plot = freq_plot_PK, width = 6, height = 4)






#----------- 4. calculate the model PK parameter using predicted PBPK parameters ------
# test which parameter are the most sensitive one
infos = combined_data %>% select("shape","Hydrodynamic_Size","coating","ZP.category","Dose","id")

# read the parameters that generated from MCMC model implement
pars_data = read.csv("/Users/wuji/work/code/Mouse-general-PBPK/plots/paras/pars_T_tot.csv")
pars_data = pars_data %>%
  mutate(Material = case_when(
    grepl("Au", Folder) ~ "Au",
    grepl("GO", Folder) ~ "GO",
    grepl("TiO2", Folder) ~ "TiO2",
    grepl("Si", Folder) ~ "Si",
    grepl("FeO", Folder) ~ "FeO",
    TRUE ~ "Other"
  ))
pars_data = merge(pars_data,infos, by.y = "id",by.x="Folder")

# the default value of each PBPK model parameter (the mean value across all the study)
default_pars = colMeans(pars_data[, sapply(pars_data, is.numeric)], na.rm = TRUE)
default_pars = log(default_pars) # log form of the parameters


# log variable that cannot be predicted by the mlr model, just use the mlr predicted results for all parameters
cannot_predict_vars_log = NULL
can_predict_vars_log = full_PBPK_df_log[!full_PBPK_df_log$Variable %in% cannot_predict_vars_log,]$Variable

#  ----------------- 5 use the mlr model to predict the PK profiles ---------
all_pred_mlr <- data.frame()
for (i in unique(all_plot_PBPK_data_log$id)){
  #i = "GO: Study2_243nm_1mg/kg"
  #i=unique(all_plot_PBPK_data_log$id)[1]
  print(i)
  all_plot_PBPK_data_log[all_plot_PBPK_data_log$id == i,]
  np.name = i
  
  Obs.df = read_observation_data(np.name)$Obs.df
  
  folder = read_observation_data(np.name)$folder
  pathway = read_observation_data(np.name)$pathway
  PDOSE = read_observation_data(np.name)$PDOSE
  
  material <- case_when(
    grepl("Au", np.name) ~ "Au",
    grepl("GO", np.name) ~ "GO",
    grepl("TiO2", np.name) ~ "TiO2",
    grepl("Si", np.name) ~ "Si",
    grepl("FeO", np.name) ~ "FeO",
    TRUE ~ "Other"
  )
  
  if (np.name == "FeO: Study2_41nm_4mg/kg") {
    tstep <- 0.5/60 #2.5/60 in the time step
  } else if (np.name %in% c("GO: Study2_243nm_1mg/kg",
                            "GO: Study2_914nm_1mg/kg_all","GO: Study2_914nm_1mg/kg_w/o_CS")) {
    # Assign the appropriate value for tstep
    tstep <- 1/60   # Fill in the appropriate value here 2/60 & 5/60 in the time points
  } else {
    tstep <- min(1, min(Obs.df$Time))
  }
  

 
  # variable that could not be predicted, use the default value for each PBPK model parameter
  cannot_predict_vars_log <- sub("_log*", "", cannot_predict_vars_log)
  cannot_vars_value = default_pars[names(default_pars) %in%cannot_predict_vars_log] # default value for all study
  

  
  # variable that could be predicted by the mlr model
  features = all_plot_PBPK_data_log[all_plot_PBPK_data_log$id == i,] # PBPK model parameters from the mlr model
  summary(lm(features$Observed~features$Predicted_Full))$adj.r.squared
  can_vars_value = features[features$Variable %in% can_predict_vars_log,][c("Predicted_Best","Variable")]
  can_vars_value = setNames(can_vars_value$Predicted_Best, can_vars_value$Variable)
  names(can_vars_value) <- sub("_log*", "", names(can_vars_value))
  
  #merged can_vars_value & cannot_predict_vars
  merged_vars_value <- c(cannot_vars_value, can_vars_value)

  # variables from MCMC model
  real_vars_vale <- pars_data[pars_data$Folder == i,]
  real_vars_vale = real_vars_vale[, sapply(real_vars_vale, is.numeric)]
  real_vars_vale = log(unlist(real_vars_vale))


  output = pred.mouse.iv(merged_vars_value) # parameters needs to be log transformed before input to the PBPK model


  
  p.r.L <-gen_obs_point_pred_plot(output, "Liver", "CL")+theme(legend.position = "none")+ xlab(NULL)
  p.r.K <- gen_obs_point_pred_plot(output, "Kidney", "CK")+theme(legend.position = "none")+ xlab(NULL)
  p.r.lung <- gen_obs_point_pred_plot(output, "Lung", "Clung")+theme(legend.position = "none")+ xlab(NULL)
  p.r.S <- gen_obs_point_pred_plot(output, "Spleen", "CS")+theme(legend.position = "none")+ xlab(NULL)
  
  required_columns <- c("CK", "CL", "Clung", "CS")
  library(patchwork)
  # Check if all required columns are present in the dataframe
  if (all(required_columns %in% colnames(Obs.df))) {
    # Arrange the four plots together
    MLR_plot_tot = p.r.L + p.r.K +p.r.S+p.r.lung+ plot_layout(2, 2) 
    
  } else if (!("CK" %in% colnames(Obs.df))) {
    MLR_plot_tot = p.r.L +p.r.S+p.r.lung+ plot_layout(2, 2) 
    
  } else if (!("CS" %in% colnames(Obs.df))) {
    MLR_plot_tot = p.r.L + p.r.K +p.r.lung+ plot_layout(2, 2) 
    
  }
  MLR_plot_tot = MLR_plot_tot+
    labs(tag = "Time (h)") +
    theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = "bottomleft"
    )
  #ggsave(filename = paste0(mlr_result_folder,"plots/pred_log/",gsub("/", "", i),"_mlr_pred.png"), plot = MLR_plot_tot, width = 6, height = 6)
  
  # compare predict & observed
  filtered_data <- output[output$Time %in% Obs.df$Time, ]
  
  # Define organs
  df_obs_pred  = subset(Obs.df, select = "Time")
  for (organ in colnames(Obs.df)[colnames(Obs.df)!= "Time"]){
    df_obs_pred[,paste0(organ,"_obs")] = Obs.df[,organ]
    df_obs_pred[,paste0(organ,"_pred")] = filtered_data[,organ]}
  
  #df_obs_pred = tibble(df_obs_pred)

  
  #-------PLOT OBS VS PRED----
  df_long <- df_obs_pred %>%
    pivot_longer(cols = -Time, names_to = c("Variable", "Type"), names_pattern = "(.*)_(.*)", values_to = "Value") %>%
    pivot_wider(names_from = Type, values_from = Value, values_fn = list)
  # set to numeric columns
  df_long <- df_long %>% 
    mutate(obs = as.numeric(obs)) %>%
    filter(!is.na(obs))
  
  df_long <- df_long %>%
    mutate(pred = as.numeric(pred)) %>%
    filter(!is.na(pred))
  

  df_long <- df_long %>% filter(obs != 0)
  df_long <- df_long %>% arrange(Variable) # ranking by organs
  

 
  for (organ in colnames(Obs.df)[colnames(Obs.df)!= "Time"]) {
    organ_data <- df_long[df_long$Variable == organ, ]
    
    # Filter the dataframe for the current organ
    df_long[df_long$Variable == organ, "adj.R"]  = summary(lm(obs~pred,data=organ_data))$adj.r.squared
    df_long[df_long$Variable == organ, "NRMSE"]  = rmse(organ_data$obs,organ_data$pred)/mean(organ_data$obs)
    
  }
  
  
  df_long$Variable <- ifelse(df_long$Variable == "CL", "Liver", 
                             ifelse(df_long$Variable == "CK", "Kidney",
                                    ifelse(df_long$Variable == "CS", "Spleen",
                                           ifelse(df_long$Variable == "Clung", "Lung", df_long$Variable))))
  df_long["id"] = i
  df_long["NRMSE_tot"] = rmse(df_long$obs,df_long$pred)/mean(df_long$obs)
  df_long["adj.R_tot"] = summary(lm(obs~pred,data=df_long))$adj.r.squared
  df_long["Ratio"] = df_long$pred / df_long$obs
  all_pred_mlr <- bind_rows(all_pred_mlr, df_long)
  
  
  
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
    scale_x_log10(limits = c(min(min(df_long$pred),min(df_long$pred)), max(max(df_long$pred),max(df_long$pred)))) +  # Set x-axis to log scale
    scale_y_log10(limits = c(min(min(df_long$pred),min(df_long$pred)), max(max(df_long$pred),max(df_long$pred)))) +   # Set y-axis to log scale
    annotate("text", x = 2* min(min(df_long$pred),min(df_long$pred)), y = 0.65*max(max(df_long$pred),max(df_long$pred)), 
             label = paste0("Adj.R^2 = ", round(unique(df_long["adj.R_tot"]), digits = 2),
                            "\nNRMSE = ", round(unique(df_long["NRMSE_tot"]), digits = 2)),
             size = 4.5, color = "black", hjust = 0, vjust = 0, fontface = "bold")+
    coord_fixed(ratio = 1)
  
  # Save the plot to a PNG file
  #ggsave(filename = paste0(mlr_result_folder,"plots/pred_log/",gsub("/", "", i),"_mlr.png"), plot = obs_pred_plot, width = 6, height = 6)
  
  ratio_plot <- ggplot(df_long, aes(x = pred, y = Ratio, shape = Variable)) +
    geom_point(size = 3) +
    labs(
      x = "Predicted value for Organ (ng/g)",
      y = "Ratio of Prediction-to-Observation") +
    theme_minimal() +
    geom_hline(yintercept = c(0.5, 2), linetype = "dashed", color = "red") +  # Add dashed lines
    scale_shape_manual(values = c(1,2,3,4))  +  # Define different shapes for each organ
    scale_y_continuous(trans = "log10", limits = c(1e-4, 1e4)) + 
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
  #ggsave(filename = paste0(mlr_result_folder,"plots/pred_log/",gsub("/", "", i),"_ratio_n.png"), plot = ratio_plot, width = 6, height = 6)
  
  }

all_pred_mlr
#write.csv(all_pred_mlr,file=paste0(mlr_result_folder,"MLR_PBPK_pred_obs_log_all.csv"),row.names = FALSE)

unique(all_pred_mlr[c("adj.R_tot","id","Variable","adj.R")])

case_acc_mlr = unique(all_pred_mlr[c("adj.R_tot","id")])
#write.csv(case_acc_mlr,file=paste0(mlr_result_folder,"MLR_PBPK_pred_obs_case.csv"),row.names = FALSE)

# the rank plot for each experiment's accuracy value
all_pred_mlr$Ratio <- all_pred_mlr$pred / all_pred_mlr$obs

all_pred_mlr_plot = ggplot(case_acc_mlr, aes(x = reorder(id,adj.R_tot ), y = adj.R_tot)) +
  geom_col(aes(fill = adj.R_tot < 0.6)) +  # Color bars based on the condition
  scale_fill_manual(values = c("gray", "black"), guide = "none") +  # Set colors: red for < 0.6, gray otherwise
  labs(x = "Variable", y = "adj.Rsquared Value") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +  # Keep y-axis labels horizontal
  coord_flip() + # Flip the coordinates for a horizontal bar plot
  theme_bw()
all_pred_mlr_plot
#ggsave(filename = paste0(mlr_result_folder,"plots/pred_log/","MLR_PBPK_log_pred_acc.png"), 
#       plot = all_pred_mlr_plot, width = 6, height = 4)


adj.R_tot = mean(unique(all_pred_mlr[c("adj.R_tot","id")])$adj.R_tot)
NRMSE_tot = mean(unique(all_pred_mlr[c("NRMSE_tot","id")])$NRMSE_tot)

obs_pred_tot_plot <- ggplot(all_pred_mlr, aes(x = obs, y = pred, shape = Variable)) +
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
  annotate(
    "text", x = 0.3, y = 1e5, 
    label = bquote(atop("avg_adj -" ~ R^2 ~ ":"  ~ .(round(adj.R_tot, 2)),
                             "  avg_NRMSE : " ~ .(round(NRMSE_tot, 2)))),fontface="bold",
    size = 4, color = "black", hjust = 0, vjust = 0
  ) +
  scale_x_log10(limits = c(1e-1, 1e6))+  # Set x-axis to log scale
  scale_y_log10(limits = c(1e-1, 1e6)) +   # Set y-axis to log scale
  coord_fixed(ratio = 1)
obs_pred_tot_plot
#ggsave(paste0(mlr_result_folder,"plots/pred_log/OBS_PRED_ALL_MLR.png"), plot = obs_pred_tot_plot,width = 6, height = 6)


# Calculate the percentage
f2_error <- (sum(all_pred_mlr$Ratio > 0.5 & all_pred_mlr$Ratio < 2) / nrow(all_pred_mlr)) * 100
f2_error
# Calculate the percentage
f3_error <- (sum(all_pred_mlr$Ratio > 1/3 & all_pred_mlr$Ratio < 3) / nrow(all_pred_mlr)) * 100
f3_error

mlr_ratio_tot <- ggplot(all_pred_mlr, aes(x = pred, y = Ratio, shape = Variable)) +
geom_point(size = 3) +
labs(
  x = "Predicted value for Organ (ng/g)",
  y = "Ratio of Prediction-to-Observation") +
theme_minimal() +
geom_hline(yintercept = c(0.5, 2), linetype = "dashed", color = "red") +  # Add dashed lines
scale_shape_manual(values = c(1,2,3,4))  +  # Define different shapes for each organ
scale_y_continuous(trans = "log10", limits = c(1e-4, 1e4)) + 
scale_x_continuous(trans = "log10") + 
  annotate("text", x = 1, y = 100, 
           label = bquote(atop("2f_error :" ~ .(sprintf("%.1f",f2_error))~"%", 
                               "3f_error :" ~ .(sprintf("%.1f",f3_error))~"%")),
           size = 4.5, color = "black", hjust = 0, vjust = 0, fontface = "bold") +
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
mlr_ratio_tot
# Save the plot to a file
#ggsave(paste0(mlr_result_folder,"plots/pred_log/Ratio_OBS_PRED_ALL.png"), plot = mlr_ratio_tot,width = 6, height = 6)


library(ggplot2)

# Loop through each unique id
for (i in unique(all_pred_mlr$id)) {
  print(i)
  
  # Filter the dataframe for the current id
  temp <- all_pred_mlr[all_pred_mlr$id == i, ]
  #temp = df_long
  # Loop through each organ
  for (organ in c("Liver", "Kidney", "Spleen", "Lung")) {
    # Filter the dataframe for the current organ
    temp_organ <- temp[temp$Variable == organ, ]
    
    # Create the plot
    p <- ggplot(temp_organ, aes(x = Time)) +
      geom_point(aes(y = obs), color = "black") +
      geom_line(aes(y = pred), color = "blue", linetype = "dashed") +
      labs(
        x = "Time (h)",
        y = "Concentration (ng/mL)",
        title = paste(i, "in", organ)
      ) +
      theme_minimal()+
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
      )
    
    # Define the file name
    file_name <- paste0(mlr_result_folder,"plots/pred/", i, "_", organ, ".png")
    
    # Save the plot to a file
    ggsave(file_name, plot = p)
  }
}



obs_pred_tot_plot <- ggplot(all_pred_mlr, aes(x = obs, y = pred, shape = Variable)) +
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
  scale_x_log10(limits = c(min(min(all_pred_mlr$pred), min(all_pred_mlr$obs)), max(max(all_pred_mlr$pred), max(all_pred_mlr$obs))))+  # Set x-axis to log scale
  scale_y_log10(limits = c(min(min(all_pred_mlr$pred), min(all_pred_mlr$obs)), max(max(all_pred_mlr$pred), max(all_pred_mlr$obs)))) +  # Set y-axis to log scale
  annotate("text", x = 1, y = 1E5, 
           label = paste0("Adj.R^2 = ", round(mean(all_pred_mlr$adj.R_tot), digits = 2),
                          "\nNRMSE = ", round(mean(all_pred_mlr$NRMSE_tot), digits = 2)),
           size = 4.5, color = "black", hjust = 0, vjust = 0, fontface = "bold")+
  coord_fixed(ratio = 1)
obs_pred_tot_plot

# Define the file name
#file_name <- paste0(mlr_result_folder,"plots/pred/OBS_PRED_ALL_MLR.png")

# Save the plot to a file
#ggsave(file_name, plot = obs_pred_tot_plot,width = 6, height = 6)




#-----6. draw the heatmap for the nanoparticle properties----
merged_acc = merge(PBPK_per_id,case_acc_mlr,by="id")

np_properties = merge(merged_acc, dataset_info, by = "id")

# Fill NA values using coalesce for hydrodynamic size and size
np_properties <- merged_data %>%
  mutate(
    Hydrodynamic_Size = coalesce(Hydrodynamic_Size, Size),
    Size = coalesce(Size, Hydrodynamic_Size)
  )


#Nanoparticles with a zeta potential between -10 and +10 mV are considered approximately neutral, 

np_properties <- np_properties %>%
  mutate(
    ZP.category = cut(zeta_potential, breaks = c(-Inf, -10, 10, Inf),
                      labels = c("negative", "neutral", "positive"),
                      include.lowest = TRUE)
  )

np_properties <- np_properties %>%
  mutate(
    ZP.category = ifelse(is.na(ZP.category),"no_info", as.character(ZP.category))
  )


# Create a data frame to map the names to the experiment IDs
mapping_df <- data.frame(
  np_name = ls_np_name,
  experiment_id = experiment_id
)

info_acc =np_properties[, c("id","NP_core", "Hydrodynamic_Size", "coating", "ZP.category", "shape", "Dose")]
info_acc <- info_acc %>%
  rename(Coating = coating, Shape = shape)

info_acc <- info_acc %>%
  left_join(mapping_df, by = c("id" = "np_name"))


info_acc =  as.data.frame(info_acc)
row.names(info_acc) <- info_acc[, ncol(info_acc)] # Set row names to the values in the last column
info_acc <- info_acc[, -c(1, ncol(info_acc))] # Remove the last column and the first column

# Step 1: Identify the text (character) columns in the dataframe `np_properties`
text_columns <- sapply(info_acc, is.character)

# Step 2: Convert the text values to numeric categories (1, 2, 3, etc.)
info_acc[, text_columns] <- lapply(info_acc[, text_columns], function(column) as.numeric(factor(column)))

# Step 3: Verify the changes by printing the first few rows
head(info_acc)
library(pheatmap)

heatmap_NP_props_plot<- pheatmap(info_acc,
                                 scale = "column",  # Scale the columns，
                                 border="black",
                                 cutree_rows = 6,
                                 #cutree_cols = 8,
                                 treeheight_col = 50,
                                 treeheight_row = 50,
                                 cluster_cols = FALSE,
                                 clustering_distance_rows = "canberra", # measure the distance betwee two vectors, taking into account their relative magnitudes
                                 clustering_method="ward.D2", # Ward D2 uses the sum of squared differences from the centroid as the criterion to minimize when merging clusters
                                 #clustering_method = "average",  # Clustering method
                                 color = colorRampPalette(c("navy", "white", "firebrick3"))(50),  # Color palette
                                 display_numbers = FALSE,  # Display cell values
                                 fontsize = 10  # Font size
)  # Title

heatmap_NP_props_plot
#ggsave(paste0(mlr_result_folder, "heat_map_NP_props.pdf"), heatmap_NP_props_plot, 
#       width = 10, height = 6)
#_---------finished heatmap for different np properties

