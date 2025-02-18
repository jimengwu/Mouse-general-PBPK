# --------- 1. plotting functions-------

# Function to create each individual plot
create_plot <- function(plot_data, y, organ, title) {
  ggplot(data = plot_data, aes(x = Time)) +
    geom_line(aes(y = !!rlang::sym(y)), col = "firebrick", lwd = 1.5) +
    geom_point(data = Obs.df, aes(y = !!rlang::sym(y)), size = 2, col = "black") +
    labs(title = title, y = paste("Concentration in", organ, "(ng/g)"), x = "Time (h)") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 20,hjust=0.5),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_blank(),  # Remove panel background
      panel.grid = element_blank(),  # Remove grid lines
    )
}


obs_pred_plot_func_mlr <- function(data, x_col, y_col, x_col_text = "Observed Values",y_col_text = "Predicted Values",
                                   id_col= NULL, text_col=NULL,annotation_text=NULL){
  
  # Set the axis limits
  max_axis_lim = max(max(data[x_col]),max(data[y_col]))
  print(max_axis_lim)
  min_axis_lim = min(min(data[x_col]),min(data[y_col]))
  print(min_axis_lim)
  # Ensure column names are treated as symbols
  x_col_sym <- rlang::sym(x_col)
  y_col_sym <- rlang::sym(y_col)
  
  if(!is.null(id_col)) {
    id_col_sym <- rlang::sym(id_col)
  }  else {id_col_sym = 1}
  if(!is.null(text_col)) {
    text_col_sym <- rlang::sym(text_col)
  }
  
  # Extract distinct values for annotation
  #distinct_data <- data %>% distinct(across(c(!!id_col_sym, !!text_col_sym)))
  
  # Construct annotation text for each organ for the total dataset
  #annotation_text <- paste0(
  #  "Adj.R^2:\n", distinct_data[[id_col]][1], " ", round(distinct_data[[text_col]][1], 2),
  #  ", ", distinct_data[[id_col]][2], " ", round(distinct_data[[text_col]][2], 2),
  #  "\n", distinct_data[[id_col]][3], " ", round(distinct_data[[text_col]][3], 2),
  #  ", ", distinct_data[[id_col]][4], " ", round(distinct_data[[text_col]][4], 2)
  #)
  #annotation_text <- paste("adj.R^2:", round(mean(data[[text_col]]),2)) 
  if (is.null(annotation_text)){
    #annotation_text <- paste("adj.R^2:",round(summary(lm(data[[x_col]]~data[[y_col]]))$adj.r.squared,3))
    annotation_text <- bquote(atop("adj -" ~ R^2 ~ ":"~.(round(summary(lm(data[[x_col]]~data[[y_col]]))$adj.r.squared,2))))
  }
  
  
  # Create the plot
  max_plot <- data %>%
    ggplot(aes(x = !!x_col_sym, y = !!y_col_sym, shape = if (!is.null(id_col)) !!id_col_sym else NULL)) +
    geom_point(size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'red') +
    scale_shape_manual(values = c(1, 2, 3, 4)) + 
    xlab(x_col_text) +
    ylab(y_col_text) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 2, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.position = c(0.8, 0.2),  # Adjust legend position
      legend.background = element_rect(fill = "transparent"),
      legend.key = element_blank(),
      legend.spacing.x = unit(0, "cm")
    ) +
    xlim(min_axis_lim, max_axis_lim) +
    ylim(min_axis_lim, max_axis_lim) +
    annotate(
      "text", x = -Inf, y = Inf, 
      label = annotation_text,fontface="bold",
      size = 4, color = "black", hjust = -0.3, vjust = 3
    ) + 
    coord_fixed(ratio = 1)
  
  return(max_plot)
}

obs_pred_plot_func <- function(data,x_col,y_col,id_col,text_col){
  
  axis_lim = max(abs(max(data[x_col])),max(abs(data[y_col])))
  
  print(unique(data[text_col])[1])
  print((data %>% distinct(across(c(text_col,id_col))))[id_col])
  max_plot <-data %>% ggplot(aes(x = x_col, y = y_col,shape=id_col)) +
    geom_point(size=3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'red') +
    scale_shape_manual(values = c(1,2,3,4)) + 
    xlab("Observed Values") +
    ylab("Predicted Values") +
    theme_minimal()+
    theme(panel.background = element_rect(fill = "transparent"),  # Set background to transparent
          panel.border = element_rect(fill=NA,color="black", linewidth=2, linetype="solid"),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.ticks = element_line(),
          axis.ticks.length = unit(0.2,"cm"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          legend.position.inside = c(0.8, 0.2),  # Adjust legend position (x, y)
          legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
          legend.key = element_blank(),  # Remove legend key
          legend.spacing.x = unit(0, "cm")) + # Adjust horizontal spacing of legend
    xlim(-axis_lim,axis_lim) +
    ylim(-axis_lim,axis_lim)+
    annotate("text", x = -axis_lim*0.8, y = axis_lim*0.7, 
             label = paste0("Adj.R^2:\n",(data %>% distinct(across(c(text_col,id_col))))[id_col][1,]," ", round((data %>% distinct(across(c(text_col,id_col))))[text_col][1,], digits = 2),
                            ", ",(data %>% distinct(across(c(text_col,id_col))))[id_col][2,], " ", round((data %>% distinct(across(c(text_col,id_col))))[text_col][2,], digits = 2),
                            "\n",(data %>% distinct(across(c(text_col,id_col))))[id_col][3,], " ", round((data %>% distinct(across(c(text_col,id_col))))[text_col][3,], digits = 2),
                            ", ", (data %>% distinct(across(c(text_col,id_col))))[id_col][4,], " ", round((data %>% distinct(across(c(text_col,id_col))))[text_col][4,], digits = 2)
             ),
             size = 4, color = "black", hjust = 0, vjust = 0)
  return (max_plot)
}



gen_MCMC_plot <- function(data, title, ylabel) {
  
  ggplot(data) + 
    geom_ribbon(aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est, color = "95% CI"), 
                fill = "lightblue", alpha = 0.3) +
    geom_ribbon(aes(x = Time, ymin = ci_q1, ymax = ci_q3, color = "50% CI"), 
                fill = "mediumblue", alpha = 0.3) +
    geom_line(aes(x = Time, y = median_est, color = "Median"), size = 1) +
    geom_point(data = Obs.df, aes(x = Time, y = !!sym(ylabel), color = "Observed"), 
               shape = 1, fill = "white", size = 2, stroke = 2) +
    scale_color_manual(values = c("95% CI" = alpha("lightblue", 0.3), 
                                  "50% CI" = alpha("mediumblue", 0.3), 
                                  "Median" = "blue", 
                                  "Observed" = "red"),
                       labels = c("95% CI", "50% CI", "Median", "Observed")) + 
    labs(y = NULL, x = "Time (h)") +
    theme_minimal() +
    theme(
      axis.ticks = element_line(),
      axis.ticks.length = unit(0.2,"cm"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5))+
    xlim(0, round(max(data$Time),1)) +
    
    scale_y_continuous(labels = ~ sprintf(fmt = "%0.1e", .))+
    # Add the title inside the plot
    annotate("text", x = 0.9 * max(data$Time), y = 0.9 * max(data$ci_upper_est, na.rm = TRUE), 
             label = title, 
             size = 4, fontface = "bold", color = "black", hjust = 0.5)
}



gen_obs_point_pred_plot <- function(data, title, ylabel) {
  ggplot(data) +
    geom_line(aes(x = Time, y = !!rlang::sym(ylabel), color = "Predicted"), linewidth = 1) +
    geom_point(data = Obs.df, aes(x = Time, y = !!sym(ylabel), color = "Observed"), 
               shape = 1, fill = "white", size = 2, stroke = 2) +
    scale_color_manual(values = c( 
      "Predicted" = "blue", 
      "Observed" = "red"),
      labels = c( "Predicted", "Observed")) + 
    labs(y = paste(title, "(ng/g)") , x = "Time (h)") +
    theme_minimal() +
    theme(
      axis.ticks = element_line(),
      axis.ticks.length = unit(0.2,"cm"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5))+
    xlim(0, round(max(data$Time),1))
}


#------------ 2. statistical functions -------

stats <- function(x){
  r2 <- summary(x)$r.squared
  r2a <- summary(x)$adj.r.squared
  f <- summary(x)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  error <- summary(x)$sigma
  y <- c(r2, r2a, f, p, error)
  return(y)
}

gene_mlr_equation <- function (lm_model){
  coefficients <- coef(lm_model)
  variable_names <- names(coefficients)
  
  # Initialize equation string with the intercept
  equation <- paste("Y =", round(coefficients[1], 4))  # intercept
  
  # Loop through each coefficient, append to equation with corresponding variable name
  for (i in 2:length(coefficients)) {
    coef_value <- coefficients[i]
    var_name <- variable_names[i]
    
    # Handle NA coefficients (if any)
    if (!is.na(coef_value)) {
      sign <- ifelse(coef_value > 0, "+", "-")
      equation <- paste(equation, sign, abs(round(coef_value, 4)), "*", var_name)
    }
  }
  return(equation)
}

MCcost<-function (pars, obs){
  out<- pred.mouse(pars)
  cost<- modCost(model=out,obs=obs,weight='mean',x="Time")
  return(cost)
}

ls_to_formula <- function(ls){
  return(paste(unlist(strsplit(unlist(ls), " ")), collapse = "+"))
}

Pred_auc <- function (pars,PDOSE,tstep){

  F_W_Liver  = 0.055  
  F_W_Brain  = 0.017 
  F_W_Lung   = 0.007  
  F_W_Kidney = 0.017
  F_W_Spleen = 0.005 
  F_W_Plasma = 0.029  
  
  BW           = 0.02                              ## kg, body weight
  tinterval    = 1                                 ## hr, Time interval for input
  TDoses       = 1                                 ## Dose times, only one dose
  
  ## Get out of log domain
  # pars is the entire parameter set; pars [-which_sig] means to keep 
  # parameters without "sig" only, and then do exp transformation, then reassign to pars
  pars <- lapply(pars [-which_sig],exp) 
  
  ## Repeat dose exposure scenario: 
  DOSE    = PDOSE*BW            ## mg; amount of oral dose
  
  ex         <- ev(ID=1, amt= DOSE, ii=tinterval, 
                         addl=TDoses-1, cmt="MBV", replicate = FALSE)
  
  
  ## set up the exposure time
  ## Simulated for 24*365 hours after dosing, but only obtained data at 24 h
  tsamp     = tgrid(0,max(Obs.df$Time),tstep)          
  
  
  ## Get a prediction
  # The code can produce time-dependent NSC values, but at time = 0, 
  # NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  out <- 
    mod %>%
    param(pars) %>%
    update(atol = 1E-80,maxsteps = 5000000)%>%
    mrgsim_d(data = ex, tgrid = tsamp)%>%
    filter(time!=0) 
  
  outdf = cbind.data.frame (Time       = out$time, 
                            AUC_Kt_id_g     = out$AUC_Kt/(DOSE*10^6)*100,#from ng/g to %ID/g
                            AUC_Lt_id_g     = out$AUC_Lt/(DOSE*10^6)*100,
                            AUC_St_id_g     = out$AUC_St/(DOSE*10^6)*100,
                            AUC_Lut_id_g    = out$AUC_Lut/(DOSE*10^6)*100,
                            AUC_blood_id_g = out$AUC_blood/(DOSE*10^6)*100,
                            CLt_id_g        = out$Liver_t/(DOSE*10^6)*100,
                            CKt_id_g        = out$Kidney_t/(DOSE*10^6)*100,
                            CSt_id_g        = out$Spleen_t/(DOSE*10^6)*100,
                            CLungt_id_g     = out$Lung_t/(DOSE*10^6)*100,
                            CBlood_id_g      = out$Plasma/(DOSE*10^6)*100
                            ) 
  return (list("outdf"  = outdf))
  
}

unit_conversion <-function(data,DOSE,input_unit = "ng/g"){
  if (input_unit == "ng/g") {
    data = data/(DOSE*10^6)*100
  }
  return (data)
}

calculate_rmse <- function(ls_parameters, params_init, folds, pred_mouse, Obs_A1, 
                           rmse_data_tot, rmse_organ) {
  for (parameter in ls_parameters) {
    cat("i:", parameter, "\n")  # Debug output
    
    # Generate predictions for each fold of the parameter
    for (fold in folds) {
      # Calculate the new parameter values
      pars_mouse_fold <- log(c(exp(params_init[parameter]) * fold, 
                               exp(params_init[-which(names(params_init) == parameter)]))) 
      R_new_fold <- pred_mouse(pars_mouse_fold)
      
      # Filter R_new_fold to include only the rows where Time matches the Obs dataset
      R_new_fold <- R_new_fold[R_new_fold$Time %in% Obs_A1$Time, ]
      
      # Extract observed and predicted values, removing NaNs
      observed_tot <- as.numeric(unlist(Obs_A1[-1]))
      predicted_tot <- as.numeric(unlist(R_new_fold[-1]))
      non_nan_indices <- which(!is.na(observed_tot))
      
      predicted_tot <- predicted_tot[non_nan_indices]
      observed_tot <- observed_tot[non_nan_indices]
      
      # Calculate RMSE
      rmse_tot <- rmse(observed_tot, predicted_tot)
      
      # Append the results to the data frame
      rmse_data_tot <- rbind(rmse_data_tot, data.frame(Parameter = parameter,
                                                       Fold = fold,
                                                       RMSE = rmse_tot))
      
      for (organ in c("Liver", "Kidney", "Spleen", "Lung")) {
        # Extract observed values based on the organ
        observed <- switch(organ,
                           "Liver" = Obs_A1$CL,
                           "Kidney" = Obs_A1$CK,
                           "Spleen" = Obs_A1$CS,
                           "Lung" = Obs_A1$Clung)
        
        # Extract predicted values based on the organ
        predicted <- switch(organ,
                            "Liver" = as.numeric(unlist(R_new_fold$CL)),
                            "Kidney" = as.numeric(unlist(R_new_fold$CK)),
                            "Spleen" = as.numeric(unlist(R_new_fold$CS)),
                            "Lung" = as.numeric(unlist(R_new_fold$Clung)))
        
        non_nan_indices <- which(!is.na(observed))
        predicted <- predicted[non_nan_indices]
        observed <- observed[non_nan_indices]
        
        # Calculate RMSE
        rmse <- rmse(observed, predicted)
        
        # Append the results to the data frame
        rmse_organ <- rbind(rmse_organ, data.frame(Organ = organ,
                                                   Parameter = parameter,
                                                   Fold = fold,
                                                   RMSE = rmse))
      }
    }
  }
  
  return(list(rmse_data_tot = rmse_data_tot, rmse_organ = rmse_organ))
}


# ------------- 3. PBPK model related functions --------
#-----------iv------
pred.mouse.iv <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp) # important to have because we cannot have nagetive kinetic value
  
  ## Define the exposure scenario
  
  BW           = 0.02                              ## kg, body weight
  tinterval    = 1                                 ## hr, Time interval
  TDoses       = 1                                 ## Dose times, only one dose

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

#----------parameters----------
params.init <- log(c(
  K_release_Liver = 0.001,  # h-1
  K_max_Liver = 20,         # h-1
  A_cap_liver = 1000,
  K_release_GI = 0.001,     # h-1
  K_max_GI = 0.075,         # h-1
  K_GI_b = 4e-5,
  K_release_Spleen = 0.001, # h-1
  K_max_Spleen = 40,        # h-1
  K_release_Kidney = 0.0004, # h-1
  K_max_Kidney = 0.075,
  K_release_Lung = 0.003,   # h-1
  K_max_Lung = 0.075,
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

gen_obs_pred_data <- function(np.name)  {
  print(np.name)
  Obs.df = read_observation_data(np.name)$Obs.df
  folder = read_observation_data(np.name)$folder
  pathway = read_observation_data(np.name)$pathway
  PDOSE = read_observation_data(np.name)$PDOSE
  
  # Define organs
  organs <- colnames(Obs.df)[colnames(Obs.df)!= "Time"]
  
  #-------------measure the prediction accuracy----------
  MC_plots = readRDS(paste0(folder,"MCMC/combined_MCMC_mean_sd_data.rds"))
  
  # Combine the four organs data into one dataset
  combined_data <- rbind(
    cbind(Organs = "CL", MC_plots$CL),
    cbind(Organs = "CK", MC_plots$CK),
    cbind(Organs = "Clung", MC_plots$Clung),
    cbind(Organs = "CS", MC_plots$CS)
  )
  
  
  filtered_data <- combined_data[combined_data$Time %in% Obs.df$Time, ]
  
  # adjusted R squared
  df_obs_pred  = subset(Obs.df, select = "Time")
  for (organ in organs){
    df_obs_pred[,paste0(organ,"_obs")] = Obs.df[,organ]
    df_obs_pred[,paste0(organ,"_pred")] = subset(filtered_data, Organs == organ)[,3]}
  
  df_obs_pred = tibble(df_obs_pred)
  
  #-------generate observed data vs prediction data----
  df_long <- df_obs_pred %>%
    pivot_longer(cols = -Time, names_to = c("Variable", "Type"), names_pattern = "(.*)_(.*)", values_to = "Value") %>%
    pivot_wider(names_from = Type, values_from = Value, values_fn = list)
  
  df_long <- df_long %>%
    mutate(obs = as.numeric(obs)) %>%
    filter(!is.na(obs))
  
  df_long <- df_long %>%
    mutate(pred = as.numeric(pred)) %>%
    filter(!is.na(pred))
  
  df_long$Variable <- ifelse(df_long$Variable == "CL", "Liver", 
                             ifelse(df_long$Variable == "CK", "Kidney",
                                    ifelse(df_long$Variable == "CS", "Spleen",
                                           ifelse(df_long$Variable == "Clung", "Lung", df_long$Variable))))
  return (df_long)
}
