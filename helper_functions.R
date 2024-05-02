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


#todo : the tstep inside function pred.mouse
MCcost<-function (pars, obs){
  out<- pred.mouse(pars)
  cost<- modCost(model=out,obs=obs,weight='mean',x="Time")
  return(cost)
}


create_sensitivity_plots <- function(sensitivity_data, parameter_names, folder_path) {
  max_val <- max(sensitivity_data[3:ncol(sensitivity_data)])
  min_val <- min(sensitivity_data[3:ncol(sensitivity_data)])
  
  for (i in 3:length(sensitivity_data)) {
    # Plot sensitivity
    p1.m <- ggplot(sensitivity_data, aes(x = sensitivity_data[,1], y = sensitivity_data[, i], color = var)) +
      geom_line(size = 1.5) + 
      ylim(min_val, max_val) +
      labs(x = "Time (h)", y = "Normalized sensitivities of model output to parameter", 
           title = parameter_names[i - 2]) +
      theme_minimal() +
      theme(
        axis.line = element_line(linewidth = 1),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, hjust = 0.1, vjust = -1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_blank()
      )
    
    # Save plot
    ggsave(paste0("combined_init_sensitivity_plot_", parameter_names[i - 2], ".png"), 
           path = folder_path, plot = p1.m, width = 10, height = 6, dpi = 300)
  }
}

calculate_rmse <- function(ls_parameters, params_init, folds, pred_mouse, Obs_A1, 
                           rmse_data_tot, rmse_organ,tstep = 1) {
  for (parameter in ls_parameters) {
    cat("i:", parameter, "\n")  # Debug output
    
    # Generate predictions for each fold of the parameter
    for (fold in folds) {
      # Calculate the new parameter values
      pars_mouse_fold <- log(c(exp(params_init[parameter]) * fold, 
                               exp(params_init[-which(names(params_init) == parameter)]))) 
      R_new_fold <- pred_mouse(pars_mouse_fold,tstep)
      
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