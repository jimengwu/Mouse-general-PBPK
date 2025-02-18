


# ---------------------------------------------------------
# Author: Jimeng Wu
# Email: jimeng.wu@empa.ch
# Date: 2025-02-18
# Description: this function was used to generate the MCMC analysis for the mouse PBPK model
# prior parameter distribution plot, posterior parameter distribution plot, 
# MCMC results distribution plot, trace plot and probability density function plot,
# also the R value and the 95% CI of the posterior parameters were calculated and saved as a CSV file
# ---------------------------------------------------------



library(bayesplot)   # Package for MCMC traceplot
library(ggplot2)
library(FME) 
library(dplyr)  
library(gridExtra)
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
source("helper_functions.R")
source("Mouse_PBPK.R")
source("dataset_info.R")
ls_np_name = c("Au: Study1_12nm_0.85mg/kg","Au: Study1_23nm_0.85mg/kg","Au: Study1_100nm_0.85mg/kg",
               "Au: Study2_34.6nm_3mg/kg","Au: Study2_55.5nm_3mg/kg","Au: Study2_77.1nm_3mg/kg",
               "Au: Study2_82.6nm_3mg/kg","Au: Study3_27.6nm_4.26mg/kg","Au: Study3_27.6nm_0.85mg/kg",
               "Si: Study1_20nm_10mg/kg","Si: Study1_80nm_10mg/kg","GO: Study1_20nm_20mg/kg",
               "GO: Study2_243nm_1mg/kg", "GO: Study2_914nm_1mg/kg_w/o_CS",
               "TiO2: Study1_385nm_10mg/kg","TiO2: Study2_220nm_60mg/kg",
               "FeO: Study1_29nm_5mg/kg","FeO: Study2_41nm_4mg/kg")


mod <- mcode ("mouse_PBPK", mousePBPK.code)

for (np.name in ls_np_name) {

  print(np.name)
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
  
  
  if (pathway == "intraperitoneal injection") {
    pred.mouse <- pred.mouse.oral
  } else if (pathway == "intravenous injection") {
    pred.mouse <- pred.mouse.iv
  }
  
  
  #--------------------------------------plots-------------------------
  #--------------1.1 plotting for prior parameter setting--------
  theta.MCMC = readRDS(paste0(folder,'MCMC/theta.rds'))
  
  
  combined_prior_density_df = readRDS(paste0(folder,"mc_sens/pri_paras_density.rds"))
  combined_prior_density_df$Par <- gsub("K_max","K_uptake", combined_prior_density_df$Par)
  pri_mean_sd_df <- read.csv(paste0(folder,"mc_sens/pri_paras_stats.csv"))
  
  # Plot density curves with facets for each parameter
  pri_den_plot<- ggplot(combined_prior_density_df, aes(x = x, y = density, fill = "Density")) +
    geom_line(color="black") +
    geom_area(fill="#9EC7C5", alpha = 0.7)+
    #geom_area(fill = "gray70") +  # Fill the area under the curve with gray color
    facet_wrap(~ Par, scales = "free", ncol = 6) +
    labs(x="value") +
    theme_minimal()+
    theme(panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "transparent"),
          axis.text.x = element_text(angle = 30, hjust = 1)) +  # Set background to transparent
    geom_text(data = pri_mean_sd_df, aes(label = paste0("Mean:", sprintf("%.0e", Mean), 
                                                        "\nSD:", sprintf("%.0e", SD))), 
              x = Inf, y = Inf, hjust = 1, vjust = 1.5, size = 3,color="black") +  # Add mean and sd as text annotations
    scale_fill_manual(values = "gray")

  # Save the plot
  #ggsave(paste0(folder,"mc_sens/pri_density_distribution_plot.png"), pri_den_plot,
  #       width = 30, height = 18, units = "cm")
  
  

  #--------------------1.2 Densities of posterior parameter uncertainty distributions of the population mean (μ). ----------
  
  
  Mouse.MCMC_chain1        <- readRDS(file = paste0(folder,"MCMC/mouse.MCMC.rds"))[[1]] # use only first chain
  parar_all_iters = exp(Mouse.MCMC_chain1$pars) # needs to exp form of parameters
  
  ## loading the theta names
  theta.names       <- names(readRDS(file = paste0(folder,"MCMC/theta.rds")))
  
  ## Mouse posterior distributions
  M.Mouse  <- parar_all_iters %>% apply(2,mean)
  SD.Mouse <- parar_all_iters %>% apply(2,sd)
  
  
  #------ plot from the mean and sd value generated from chain for posterior parameter plotting---------
  # Initialize lists to store samples for each parameter
  post_param_samples <- list()
  
  # Number of points to use in density estimation
  n_points <- 10000  # Adjust as needed for smoother or rougher curves
  
  
  for (param_index in seq_along(M.Mouse)) {
    # Draw samples from the normal distribution for the current parameter
    post_param_samples[[param_index]] <- rnorm(n_points, 
                                               mean = M.Mouse[[param_index]], 
                                               sd = SD.Mouse[[param_index]])
  }
  
  
  # Initialize an empty list to store density estimation dataframes for each parameter
  post_density_data <- list()
  # Generate density estimation and store it for each parameter
  for (param_index in seq_along(M.Mouse)) {
    #print(param_index)
    # Generate density estimation for the current parameter
    post_density_estimation <- density(post_param_samples[[param_index]], n = n_points)
  
  # Create a dataframe for the density estimation data
    post_density_df <- data.frame(x = post_density_estimation$x, 
                                  density = post_density_estimation$y, 
                                  Par = names(M.Mouse)[param_index])
  
  # Store the density dataframe
    post_density_data[[param_index]] <-post_density_df
  }
  #-----------plot from the mean and sd value generated from chain for posterior parameter plotting---------
  
  # Combine all density dataframes into one dataframe
  combined_post_density_df <- do.call(rbind, post_density_data)
  
  combined_post_density_df$Par <- gsub("K_max","K_uptake", combined_post_density_df$Par)
  
  # from posterior parameters to generate the posterior distribution's mean and sd value for text in the plot
  

  # Calculate mean and standard deviation for each parameter
  post_mean_sd_df <- data.frame(Par = names(M.Mouse),
                                Mean = unlist(M.Mouse),
                                SD = unlist(SD.Mouse))
  
  
  # Plot density curves with facets for each parameter the nrow and ncol might need to be changed
  
  
  post_den_plot <- ggplot(combined_post_density_df %>%
                            filter(!grepl("sig", Par)), aes(x = x, y = density)) +
    geom_area(fill="#D17E5E", alpha = 0.7)+
    geom_line(color="black") +
    facet_wrap(~ Par, scales = "free",ncol = 6) +
    labs(x = "value") +
    theme_minimal() +
    theme(panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "transparent"),
          axis.text.x = element_text(angle = 30, hjust = 1)) + #"#f7f7f7"
    geom_text(data = post_mean_sd_df %>%
                filter(!grepl("sig", Par)), 
              aes(label = paste0("Mean:", sprintf("%.2e", Mean), "\nSD:", sprintf("%.1e", SD))), 
              x = Inf, y = Inf, hjust = 1, vjust = 1.5, size = 3, color = "black")  # Add mean and sd as text annotations

  prior_post_den_plot <- ggplot() +
    geom_area(data = combined_prior_density_df, aes(x = x, y = density, fill = "Prior"), color = "black", alpha = 0.7) +
    geom_area(data = combined_post_density_df %>% filter(!grepl("sig", Par)), 
              aes(x = x, y = density, fill = "Posterior"), color = "black", alpha = 0.7) +
    facet_wrap(~ Par, scales = "free", nrow = 7, ncol = 6) +
    labs(x="value") +
    theme_minimal() +
    theme(panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "transparent"),
          axis.text.x = element_text(angle = 30, hjust = 1)) + #"#f7f7f7"
    scale_fill_manual(name = "Distribution", values = c("Prior" = "#9EC7C5", "Posterior" = "#D17E5E")) + # Add legend for distribution
    guides(fill = guide_legend(title = "Distribution")) # Add legend title

  
  # Save the plot
  #ggsave(paste0(folder,"mc_sens/post_density_distribution_plot.png"), post_den_plot, 
  #       width = 30, height = 18, units = "cm")
  
  #ggsave(paste0(folder,"mc_sens/prior_post_density_distribution_plot.png"), prior_post_den_plot, 
  #       width = 30, height = 18, units = "cm")
  
  
  # Save to CSV
  #write.csv(post_mean_sd_df, file = paste0(folder,"mc_sens/post_paras_stats.csv"), row.names = FALSE)
  
  #-------- Save the posterior parameters (95% CI)----------
  MC.mouse.1 = as.mcmc (Mouse.MCMC_chain1$pars)
  quan.mouse = exp(summary(MC.mouse.1)$quantiles)  
  
  #write.csv(quan.mouse,file=paste0(folder,"mc_sens/mouse.summary_pos.csv"))
  
  #-----------------------------2. get range plots results -------------
  print("starting prediction results plotting...")
  Mouse.MCMC_chain1        <- readRDS(file = paste0(folder,"MCMC/mouse.MCMC.rds"))[[1]] # only first chain
  
  which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG
  
  Newtime.r   = pred.mouse(theta.MCMC[-which_sig])$Time  
  
  # this is the new time variable, now it has been changed to sample per day.
  nrwo.r = length (Newtime.r)
  outputlength = nrow(Mouse.MCMC_chain1$pars)
  
  # Define the matrix names
  matrix_names <- c("CL", "CK", "Clung", "CS")
  
  # Initialize lists to store matrices and plots
  MC_matrices <- list()
  MC_plots <- list()
  
  # Loop over matrix names
  for (name in matrix_names) {
    # Create the matrix
    MC_matrices[[name]] <- matrix(nrow = nrwo.r, ncol = outputlength / 1000)
    
    # Loop over columns
    for (i in 1:(outputlength / 1000)) {
      j <- i * 1000
      pars.mouse <- Mouse.MCMC_chain1$pars[j, ]
      MCdata <- pred.mouse(pars.mouse)
      MC_matrices[[name]][, i] <- MCdata[[name]]
    }
    
    # Calculate mean and standard deviation
    M <- apply(MC_matrices[[name]], 1, mean)
    SD <- apply(MC_matrices[[name]], 1, sd)
    
    # Create data for plot
    MC_plot <- cbind(
      Time = Newtime.r,
      as.data.frame(t(apply(MC_matrices[[name]], 1, function(y_est) c(
        median_est = median(y_est, na.rm = TRUE),
        ci_q1 = quantile(y_est, probs = 0.25, names = FALSE, na.rm = TRUE),
        ci_q3 = quantile(y_est, probs = 0.75, names = FALSE, na.rm = TRUE),
        ci_10 = quantile(y_est, probs = 0.1, names = FALSE, na.rm = TRUE),
        ci_90 = quantile(y_est, probs = 0.9, names = FALSE, na.rm = TRUE),
        ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE, na.rm = TRUE),
        ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE, na.rm = TRUE)
      ))))
    )
    
    # Store plot in list
    MC_plots[[name]] <- MC_plot
  }
  
  #saveRDS(MC_plots, paste0(folder,"MCMC/combined_MCMC_mean_sd_data.rds"))
  MC_plots = readRDS(paste0(folder,"MCMC/combined_MCMC_mean_sd_data.rds"))
  
  
  # Create plots for each parameter
  p.r.L <- gen_MCMC_plot(MC_plots$CL, "Liver", "CL")
  p.r.K <- gen_MCMC_plot(MC_plots$CK, "Kidney", "CK")
  p.r.lung <- gen_MCMC_plot(MC_plots$Clung, "Lung", "Clung")
  p.r.S <- gen_MCMC_plot(MC_plots$CS, "Spleen", "CS")
  
  # Modify the first two plots to remove the x-axis
  p.r.L_no_x <- p.r.L + theme(axis.title.x = element_blank())
  p.r.K_no_x <- p.r.K + theme(axis.title.x = element_blank())
  p.r.Lung_no_x <- p.r.lung + theme(axis.title.x = element_blank())
  
  #p.r.S
  
  #ggsave(paste0(folder,"MCMC/MCMC_fit_spleen.png"),width = 10, height = 8)
  required_columns <- c("CK", "CL", "Clung", "CS")
  
  # Check if all required columns are present in the dataframe
  if (all(required_columns %in% colnames(Obs.df))) {
    # Arrange the four plots together


    # Create the combined plot
    combined_MCMC_plot <- grid.arrange(
      arrangeGrob(grobs = list(p.r.L_no_x, p.r.K_no_x, p.r.lung, p.r.S),
        ncol = 2, nrow = 2),
      left = textGrob("Concentration in Organ (ng/g)",rot = 90, 
        gp = gpar(fontsize = 14, fontface = "bold")))
    
  } else if (!("CK" %in% colnames(Obs.df))) {
    # Arrange the three plots together without CK

    # Create the combined plot
    combined_MCMC_plot <- grid.arrange(
      arrangeGrob(grobs = list(p.r.L_no_x, p.r.Lung_no_x, p.r.S),
                  ncol = 2, nrow = 2),
      left = textGrob("Concentration in Organ (ng/g)",rot = 90, 
                      gp = gpar(fontsize = 14, fontface = "bold")))
    
  } else if (!("CS" %in% colnames(Obs.df))) {
    # Arrange the three plots together without CS

    # Create the combined plot
    combined_MCMC_plot <- grid.arrange(
      arrangeGrob(grobs = list(p.r.L_no_x, p.r.K_no_x, p.r.lung),
                  ncol = 2, nrow = 2),
      left = textGrob("Concentration in Organ (ng/g)",rot = 90, 
                      gp = gpar(fontsize = 14, fontface = "bold")))
  }
  
  # Save the combined plot
  #ggsave(paste0(folder, "MCMC/mcmc_fit_combined.png"), combined_MCMC_plot, width = 14, height = 10)
  
  
  
  #-----------------------------
  
  
  
  #---------------------3. trace plot and probability density function plot----------
  ## Trace plot using bayes plot
  ## Convergences plot
  print("starting convergence results plotting...")
  color_scheme_set("blue")

  combinedchains        <- readRDS(file = paste0(folder,'MCMC/mouse.MCMC.comb.rds')) # use only first chain
  R <-gelman.diag(combinedchains) # Gel man convergence diagnosis
  
  #heidel.diag (combinedchains)          # covergence diagnosis/Heidelberger and Welch's convergence diagnostic
  #gelman.plot (combinedchains)          # gelman plot
  
  
  trace_plot <- mcmc_trace (
    combinedchains,
    pars = names(theta.MCMC[1:length(theta.MCMC [-which_sig])]),
    size = 0.5,
    facet_args = list(nrow = 5)) +
    #ggplot2::scale_color_brewer() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))  # Adjust angle of x-axis labels
  #ggsave(paste0(folder,"mc_sens/trace_plot.png"),trace_plot, width = 10, height = 6)
  
  # todo probabilistic distribution of this difference with Fig3 
  plot_prob_chains <- mcmc_dens_overlay(
    combinedchains,
    pars = names(theta.MCMC[1:length(theta.MCMC [-which_sig])]),
    facet_args = list(nrow=5))  +
    #ggplot2::scale_color_brewer() +
    theme(axis.text.x = element_text(angle = 30))
 
  #ggsave(paste0(folder,"mc_sens/prob_chains.png"),plot_prob_chains, width = 10, height = 6)
  
  #write.csv(R[1],paste0(folder,"mc_sens/r.csv"))
  #write.csv(R[2],paste0(folder,"mc_sens/r_summary.csv"))
}


R_df <- data.frame(Iteration = character(), R_value = numeric(),
                   mean_Point_est = numeric(), 
                   mean_Upper_ci = numeric(), stringsAsFactors = FALSE)


for (np.name in ls_np_name) {
  print(np.name)
  folder = read_observation_data(np.name)$folder
  combinedchains        <- readRDS(file = paste0(folder,'MCMC/mouse.MCMC.comb.rds')) # use only first chain
  R <-gelman.diag(combinedchains) # Gel man convergence diagnosis
  R_df <- rbind(R_df, data.frame(Iteration = np.name, NRMSE_value = R[2],
                                 mean_Point_est = colMeans(R[1]$psrf, na.rm = TRUE)[1],
                                 mean_Upper_ci = colMeans(R[1]$psrf, na.rm = TRUE)[2]))
}
folder_tot = "plots/"
#write.csv(R_df, file= paste0(folder_tot,"R_hat_tot.csv"),row.names = FALSE)

