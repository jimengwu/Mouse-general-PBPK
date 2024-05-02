########################################## Fig. 6 ################################ 
# circular plot of sensitivity analysis                                           #
###################################################################################
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

folder = "/Users/mmm/work/Mouse-PBPK/plots/100nm_short_nsc/"
Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_100nm_short.csv")  
## Loading pbpk model 
mod         <- mcode ("mouse_PBPK", mousePBPK_100nm.code)

folder = 'plots/13nm_short_nsc/'
#Obs.A1 <- read.csv(file ="dataset/tk/mouse/R_input_mouse_study1_13nm_short.csv") 
mod         <- mcode ("mouse_PBPK", mousePBPK.code)


## Input PBPK model
#micePBPK.code     <- readRDS (file = "micePBPK.RDS")

folder = 'plots/13nm_nsc_GI/'
## Loading human, rat, mouse, monkey MCMC data
Mouse.MCMC        <- readRDS(file = paste0(folder,"MCMC/mouse.MCMC.rds"))

## loading the theta names
theta             <- readRDS(file = paste0(folder,"MCMC/theta.rds"))
theta.names       <- names(theta)
which_sig         <- grep("sig", theta.names)
tend              <- max(Obs.A1$Time)
## Sensitivity analysis
Pred <- function (pars.mouse){
  
  #which_sig <- grep("sig", theta.names)
  ## Get out of log domain
  # pars.mouse is the entire parameter set; pars.mouse [-which_sig] means to keep 
  # parameters with "sig" only, and then do exp transformation, then reassign to pars.mouse
  pars.mouse <- lapply(pars.mouse [-which_sig],exp) 
  
  ## Repeat dose exposure scenario: 
  
  BW.mouse          = 0.02                               ## mouse body weight

  
  tinterval         = 0.1                                  ## Time interval
  TDoses.mouse      = 1                                   ## The number of dosing in the mouse

  
  PDOSEoral.mouse   = 0.85                                   ## mg/kg; BW Oral dose, Change et al., 2012

  DOSEoral.mouse    = PDOSEoral.mouse*BW.mouse            ## mg; amount of oral dose

  ex.mouse         <- ev(ID=1, amt= DOSEoral.mouse, ii=tinterval, 
                         addl=TDoses.mouse-1, cmt="MBV", replicate = FALSE)

  
  ## set up the exposure time
  ## Simulated for 24*365 hours after dosing, but only obtained data at 24 h
  tsamp.mouse     = tgrid(0,max(Obs.A1$Time),tinterval)          

  
  ## Get a prediction
  # The code can produce time-dependent NSC values, but at time = 0, 
  # NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  out.mouse <- 
    mod %>%
    param(pars.mouse) %>%
    update(atol = 1E-80,maxsteps = 5000000)%>%
    mrgsim_d(data = ex.mouse, tgrid = tsamp.mouse)%>%
    #Req(AUC_Liver_t)
    filter(time!=0) 
  
  outdf.mouse = cbind.data.frame (Time       = out.mouse$time, 
                                  AUC_Kt     = out.mouse$AUC_Kt,
                                  AUC_Lt     = out.mouse$AUC_Lt,
                                  AUC_St     = out.mouse$AUC_St,
                                  AUC_Lut    = out.mouse$AUC_Lut,
                                  CLt        = out.mouse$Liver_t,
                                  CKt        = out.mouse$Kidney_t) 
  return (list("outdf.mouse"  = outdf.mouse))
  
}

pars.mouse  = Mouse.MCMC[[1]]$bestpar


#pars.mouse = theta.MCMC
R = Pred(pars.mouse)

## Create the matrix for normalized sensitivity coefficient data, 4 columns each represent one organ
NSC_mouse   = matrix(nrow=length(pars.mouse[-which_sig]),ncol=4)
NSC_mouse_end = NSC_mouse

percentage = 0.01

for (i in 1:length(pars.mouse[-which_sig])) {
  NSC_mouse_T   = matrix(nrow=length(pars.mouse[-which_sig]),ncol=length(R$outdf.mouse$Time))
  for (j in 2:5){
        # Each cycle, generate a new value of parameter i (e.g., 10.0a), 
        # and delete parameter i, so that you can proceed to the next parameter i+1
        pars.mouse.new      <- log(c(exp(pars.mouse[i])*(1 + percentage),exp(pars.mouse[-i]))) 
        Rnew                <- Pred(pars.mouse.new)
        delta.P.mouse       <- exp(pars.mouse[i])/(exp(pars.mouse[i])*percentage) # is exp here needed?
      
        ## Estimated the AUC
        
        # get the NSC for the whole time range
        delta.AUC.CLt.mouse.T = (Rnew$outdf.mouse %>% select (names(R$outdf.mouse)[j]) 
                                - R$outdf.mouse %>% select (names(R$outdf.mouse)[j]))
        
        NSC_mouse_T [i,] <-as.numeric(unlist(delta.AUC.CLt.mouse.T/ 
                                               (R$outdf.mouse %>% select (j))) 
                                      *delta.P.mouse)
    
        
        #get the NSC median value for the whole time range
        #NSC_mouse   [i, 1]      <- as.numeric((delta.AUC.CLt.mouse/Mouse.AUC.CLt.ori) * delta.P.mouse)
        NSC_mouse   [i, j-1]      <- NSC_mouse_T [i,which(Rnew$outdf.mouse$Time == 1)]
        NSC_mouse_end [i, j-1]    <- NSC_mouse_T [i,which(Rnew$outdf.mouse$Time == tend)]
}}



colnames (NSC_mouse)      = c("NSC_AUC_CKt","NSC_AUC_CLt","NSC_AUC_CSt","NSC_AUC_CLut") 
colnames (NSC_mouse_end)  =  c("NSC_AUC_CKt","NSC_AUC_CLt","NSC_AUC_CSt","NSC_AUC_CLut") 

rownames(NSC_mouse)       = theta.names[1:length(pars.mouse[-which_sig])]
rownames(NSC_mouse_T)     = theta.names[1:length(pars.mouse[-which_sig])]
rownames(NSC_mouse_end)   = theta.names[1:length(pars.mouse[-which_sig])]


##################################### Circle barplot function ###############################################
## plot modifed from "R graph gallery: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/ "  #
#############################################################################################################

Circle.plot <- function (melt.data){ # melt.data is an argument of Circle.plot function.
  melt.data= melt.data[order(melt.data$value, decreasing = TRUE), ]
  melt.data$id=seq(1, nrow(melt.data)) # id is the number of rows. In total, there were 68 rows.
  
  
  # Get the name and the y position of each label
  label_data=melt.data
  number_of_bar=nrow(label_data) # in total, there were 30 rows, 30 parameters * 1 species
  # substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  angle= 90 - 360 * (label_data$id - 0.5) /number_of_bar     
  label_data$hjust<-ifelse( angle < -90, 1, 0)
  label_data$angle<-ifelse(angle < -90, angle+180, angle)

  # prepare a data frame for base lines
  base_data=melt.data %>%
  #group_by(group) %>%
  summarize(start=min(id), end=max(id)) %>%
  rowwise() %>%
  mutate(title=mean(c(start, end)))
  
  # Make the plot
  #windowsFonts(Times=windowsFont("Times New Roman"))
 
  p.cir.plot <- 
    ggplot(melt.data, aes(x = as.factor(id), y = (value * 100))) +       
    geom_bar(aes(x = as.factor(id), y = abs(value * 100)), 
             stat = "identity", alpha = 0.6,fill = ifelse(melt.data$value < 0, "red", "black")) +
    ylim(-100, 200) +  # Setting the y-axis limit
    theme_minimal() +
    theme(
      legend.position = "none",
      text = element_text(family = "Times"),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1, 4), "cm")
    ) +
    coord_polar() +
    geom_text(data = label_data, aes(x = id, y = 10, label = par, hjust = hjust),
              fontface = "bold", alpha = 1, size = 10,
              angle = label_data$angle,color = "black") + 
    geom_text(data = label_data, aes(x = id, y = - 75, label = paste(round(value * 100, 1), "%"), 
              hjust = hjust), color = ifelse(label_data$value < 0, "red", "black")
              ,fontface = "bold", alpha = 1, size = 10,
              angle = label_data$angle)  +
    geom_segment(data = base_data, aes(x = start, y = -5, xend = end+1, yend = -5), 
                 colour = "black", alpha = 0.8, size = 1.0)
  

  
return (p.cir.plot)
}
library(gridExtra)
###################### Fig. 8a; AUC of plasma ####################
# Initially, NSC_mouse is a dataset with two rows. melt can reshape it to a dataset with two columns.
plot_list <- list()
for (i in 1:4){
  melt.mouse.t1         = melt(NSC_mouse[,i]) 
  melt.mouse.t1$par     = rownames(NSC_mouse) # Get the row names (parameter names), repeat four times, so that you can the parameter names for each of the four species
  p1.t1                = Circle.plot (melt.mouse.t1%>%filter(abs(value)>0.01))+
  ggtitle(paste("Sensitity to", colnames(NSC_mouse)[i])) + # Set title for each plot
    theme(plot.title = element_text(vjust = -10,hjust=0.5))
  #p1.t1
  # Save the plot to the list
  plot_list[[i]] <-p1.t1
  
  if (i %% 4 == 0) {
    # Combine every three plots into one plot
    combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 2, nrow =2))
    # Apply color scale manually to the combined plot
    
    # Save or print combined_plot as desired
    ggsave(paste0("mc_sens/", "pars_sens.png"), path = folder, combined_plot, width = 20, height = 8)
    # Clear the plot list for the next iteration
    plot_list <- list()
  }
}




melt.mouse.end         = melt(NSC_mouse_end[,1]) 
#melt.mouse.CA$group  = c("Mouse") # Add a third column of group with the values of Mouse for all rows.
melt.mouse.end$par     = rep(rownames(NSC_mouse_end),1) # Get the row names (parameter names), repeat four times, so that you can the parameter names for each of the four species

p1.e                = Circle.plot (melt.mouse.end%>%filter(abs(value)>0.001))
p1.e


plot_list <- list()
for (i in 1:4){
  melt.mouse.end         = melt(NSC_mouse_end[,i]) 
  melt.mouse.end$par     = rownames(NSC_mouse) # Get the row names (parameter names), repeat four times, so that you can the parameter names for each of the four species
  p1.e                = Circle.plot (melt.mouse.end%>%filter(abs(value)>0.01))+
                    ggtitle(paste("Sensitity to", colnames(NSC_mouse)[i])) + # Set title for each plot
                    theme(plot.title = element_text(vjust = -10,hjust=0.5),
                          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
  # Save the plot to the list
  plot_list[[i]] <-p1.e
  
  ggsave(paste0("mc_sens/", colnames(NSC_mouse)[i],"_pars_sens_end.png"), 
         path = folder, p1.e, width = 24, height = 20)
  if (i %% 4 == 0) {
    # Combine every three plots into one plot
    combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 2, nrow =2))
    # Apply color scale manually to the combined plot
    
    # Save or print combined_plot as desired
    ggsave(paste0("mc_sens/", "pars_sens_end.png"), path = folder, combined_plot, width = 24, height = 20)
    # Clear the plot list for the next iteration
    plot_list <- list()
  }
}



####### Save plot #######
#ggsave("Sens@t1.tiff",scale = 1,
#       plot = p1.t1,
#       path = paste0(folder,"post_sens/"),
#       width = 45, height = 20, units = "cm",dpi=320)

#ggsave("Sens_end.tiff",scale = 1,
#       plot = p1.m,
#       path = paste0(folder,"post_sens/"),
#       width = 45, height = 20, units = "cm",dpi=320)

#dev.off()

#-----------------------local sensitivity analysis------------------



## Create the matrix for normalized sensitivity coefficient data
LSA_mouse_double <- matrix(nrow = length(pars.mouse[-which_sig]), ncol = length(R$outdf.mouse$Time))
LSA_mouse_half   <- matrix(nrow = length(pars.mouse[-which_sig]), ncol = length(R$outdf.mouse$Time))



# Initialize a list to store the plots
plot_list <- list()


#-----------for Liver-------------------
for (i in 1:length(pars.mouse[-which_sig])) {

  # Each cycle, generate a new value of parameter i (e.g., 10.0a), 
  # and delete parameter i, so that you can proceed to the next parameter i+1
  pars.mouse.double <- log(c(exp(pars.mouse[i]) * 2, exp(pars.mouse[-i]))) 
  Rnew.double <- Pred(pars.mouse.double)
  
  pars.mouse.half <- log(c(exp(pars.mouse[i]) / 2, exp(pars.mouse[-i]))) 
  Rnew.half <- Pred(pars.mouse.half)
  
  # Extract CLt values
  LSA_mouse_double[i,] <- as.numeric(unlist(Rnew.double$outdf.mouse$CLt))
  LSA_mouse_half[i, ] <- as.numeric(unlist(Rnew.half$outdf.mouse$CLt))
  
  # Combine the data frames into a single data frame
  combined_data <- data.frame(Time = R$outdf.mouse$Time,
                              half = LSA_mouse_half[i,],
                              double = LSA_mouse_double[i,],
                              predicted = R$outdf.mouse$CLt)
  
  # Plot using ggplot
  p1.m <- ggplot(combined_data, aes(x = Time)) +
    geom_point(data = Obs.A1, aes(Time, CL), size = 1.5, color = "black") + 
    geom_line(aes(y = half), color = "#FFC107", alpha=0.5,linewidth = 1,linetype="dashed") +  # Green color for LSA_mouse_half
    geom_line(aes(y = double), color = "#FFC107",alpha=0.5, linewidth = 1,linetype="dashed") +  # Yellow color for LSA_mouse_double
    geom_line(aes(y = predicted), color = "#2196F3",alpha=0.7, linewidth = 1) +  # Blue color for R_outdf_mouse_CLt
    ylim(0, 15000) + 
    labs(x = "Time(h)", y = "Concentration in the liver (ng/g)",title=names(pars.mouse[-which_sig])[i]) +
    theme_minimal()+
    theme(
      axis.line = element_line(linewidth = 1),  # Increase axis width
      axis.text = element_text(size = 20),
      axis.title = element_text(size=20),
      plot.title = element_text(size = 20,  hjust = 0.1, vjust = -8),  # Adjust title position
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  
  # Save the plot to the list
  plot_list[[if (i %% 3 == 0) 3 else i %% 3]] <- p1.m

  # Save every three plots into one plot
  if (i %% 3 == 0) {
    # Combine every three plots into one plot
    combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 3, nrow = 1))
    # Apply color scale manually to the combined plot
    
    # Save or print combined_plot as desired
    ggsave(paste0("mc_sens/LSA_CL_", i/3, ".png"), path = folder, combined_plot, width = 20, height = 6)
    # Clear the plot list for the next iteration
    plot_list <- list()
  }
}



#-----------for Kidney-------------------
for (i in 1:length(pars.mouse[-which_sig])) {
  
  # Each cycle, generate a new value of parameter i (e.g., 10.0a), 
  # and delete parameter i, so that you can proceed to the next parameter i+1
  pars.mouse.double <- log(c(exp(pars.mouse[i]) * 2, exp(pars.mouse[-i]))) 
  Rnew.double <- Pred(pars.mouse.double)
  
  pars.mouse.half <- log(c(exp(pars.mouse[i]) / 2, exp(pars.mouse[-i]))) 
  Rnew.half <- Pred(pars.mouse.half)
  
  # Extract CKt values
  LSA_mouse_double[i,] <- as.numeric(unlist(Rnew.double$outdf.mouse$CKt))
  LSA_mouse_half[i, ] <- as.numeric(unlist(Rnew.half$outdf.mouse$CKt))
  
  # Combine the data frames into a single data frame
  combined_data <- data.frame(Time = R$outdf.mouse$Time,
                              half = LSA_mouse_half[i,],
                              double = LSA_mouse_double[i,],
                              predicted = R$outdf.mouse$CLt)
  
  # Plot using ggplot
  p1.m <- ggplot(combined_data, aes(x = Time)) +
    geom_point(data = Obs.A1, aes(Time, CK), size = 1.5, color = "black") + 
    geom_line(aes(y = half), color = "#FFC107", alpha=0.5,linewidth = 1,linetype="dashed") +  # Green color for LSA_mouse_half
    geom_line(aes(y = double), color = "#FFC107",alpha=0.5, linewidth = 1,linetype="dashed") +  # Yellow color for LSA_mouse_double
    geom_line(aes(y = predicted), color = "#2196F3",alpha=0.7, linewidth = 1) +  # Blue color for R_outdf_mouse_CLt
    ylim(0, 15000) + 
    labs(x = "Time(h)", y = "Concentration in the liver (ng/g)",title=names(pars.mouse[-which_sig])[i]) +
    theme_minimal()+
    theme(
      axis.line = element_line(linewidth = 1),  # Increase axis width
      axis.text = element_text(size = 20),
      axis.title = element_text(size=20),
      plot.title = element_text(size = 20,  hjust = 0.1, vjust = -8),  # Adjust title position
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  
  # Save the plot to the list
  plot_list[[if (i %% 3 == 0) 3 else i %% 3]] <- p1.m
  
  # Save every three plots into one plot
  if (i %% 3 == 0) {
    # Combine every three plots into one plot
    combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 3, nrow = 1))
    # Apply color scale manually to the combined plot
    
    # Save or print combined_plot as desired
    ggsave(paste0("mc_sens/LSA_CK_", i/3, ".png"), path = folder, combined_plot, width = 20, height = 6)
    # Clear the plot list for the next iteration
    plot_list <- list()
  }
}


