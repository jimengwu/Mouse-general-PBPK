
# ---------------------------------------------------------
# Author: Jimeng Wu
# Email: jimeng.wu@empa.ch
# Date: 2025-02-18
# Description:this script was generated  to compare two settings of parameter, 
# one generated from one case and another generated from another, however,
# the two cases was got from the same NP setting, but different time range. 
# longer time range could cover the short time range, but vice versa is not the case.

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

mod <- mcode ("mouse_PBPK", mousePBPK.code)



print(np.name)
np.name = "Au: Study3_27.6nm_0.85mg/kg" 
np.name_other = "Au: Study1_23nm_0.85mg/kg" 


np.name_other = "Au: Study3_27.6nm_0.85mg/kg" 
np.name = "Au: Study1_23nm_0.85mg/kg" 



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


cal_out <- function(np.name) {
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
  
  theta.MCMC = read.csv(file = paste0(folder,"mc_sens/post_paras_stats.csv"))
  theta.MCMC = setNames(theta.MCMC$Mean, theta.MCMC$Par)
  names(theta.MCMC) <- gsub("K_max","K_uptake", names(theta.MCMC))
  which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG
  
  results = pred.mouse(log(theta.MCMC[-which_sig]))
  

  return(results)
  
}


results_other_pars = cal_out(np.name_other)
results_own_pars = cal_out(np.name)

create_double_plot <- function(plot_data_1,plot_data_2, y, organ, title,legend= FALSE) {
  ggplot(data = plot_data_1, aes(x = Time)) +
    geom_line(aes(y = !!rlang::sym(y), color = "parameters from its own"),lwd = 1) +
    geom_line(data = plot_data_2, aes(y = !!rlang::sym(y),color = "parameters from another"), lwd = 1) +  # Second line
    geom_point(data = Obs.df, aes(y = !!rlang::sym(y)), size = 2, col = "black") +
    labs(title = title, x = "Time (h)") +
    scale_color_manual(values = c("parameters from its own" = "#B3001B", 
                                  "parameters from another" = "#255C99",
                                  "observed" = "black")) +
    guides(color = guide_legend(title = NULL)) +  # Remove legend box title
    theme_minimal() +
    theme(
      #legend.position = if (legend) c(0.45,0.9) else "none",
      
      legend.position = if (legend) c(0.6,0.6) else "none", # for S1
      legend.text = element_text(size=14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_blank(),  # Remove y-axis title
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 18,hjust=0.5),
      panel.background = element_blank(),  # Remove panel background
      panel.grid = element_blank(),  # Remove grid lines
      axis.line = element_line(size = 0.5, color = "black"),  # Add axis lines
      axis.ticks.length = unit(0.2, "cm"),  # Customize tick length
      axis.ticks.x = element_line(size = 0.5, color = "black"),  # X-axis ticks
      axis.ticks.y = element_line(size = 0.5, color = "black"),  # Y-axis ticks
      panel.border = element_rect(color = "black", fill = NA, size = 1.5) # Remove the panel border
    )

}


plot_liver.a1   <- create_double_plot(results_own_pars,results_other_pars, y = "CL", organ = "Liver", title = "Liver")
plot_kidney.a1  <- create_double_plot(results_own_pars,results_other_pars, y = "CK", organ = "Kidney", title = "Kidney")
plot_lung.a1    <- create_double_plot(results_own_pars,results_other_pars, y = "Clung", organ = "Lung", title = "Lung")
plot_spleen.a1  <- create_double_plot(results_own_pars,results_other_pars, y = "CS", organ = "Spleen", title = "Spleen",legend = TRUE)

library(grid)
y_axis_title <- textGrob("Concentration in Organs (ng/g)", rot = 90, gp = gpar(fontsize = 20))
plot.a1 <- grid.arrange(
  arrangeGrob(y_axis_title, ncol = 1, nrow = 1),
  arrangeGrob(plot_liver.a1, plot_kidney.a1, plot_lung.a1, plot_spleen.a1, ncol = 4),
  ncol = 2,
  widths = c(1, 20)  # Adjust width ratio as needed
)


plot.a1
#ggsave(paste0(folder, "double:",strtrim(np.name, 20),"_",strtrim(np.name_other, 20),".png"), plot.a1,
#       width = 20, height = 5)
