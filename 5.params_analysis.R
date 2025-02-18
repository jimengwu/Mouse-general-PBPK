
# ---------------------------------------------------------
# Author: Jimeng Wu
# Email: jimeng.wu@empa.ch
# Date: 2025-02-18
# Description:this function is used to analyze the posterior parameter distribution and the sensitivity of the parameters 
# ---------------------------------------------------------


#------------------------------parameter distribution-----
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(plyr)
library(reshape)
library(scales)
library(ggbreak)
library("ggExtra")
library(ggridges)
library(ggrepel)
library(readr)
library(readxl)
library(stringr)
library(RColorBrewer)
source("helper_functions.R")
source("Mouse_PBPK.R")
source("dataset_info.R")


# Create a function to generate data points for each distribution
generate_data <- function(mean, sd, n = 1000) {
  data <- rnorm(n, mean, sd)
  return(data)
}


#-------------------------0.1 read parameter distribution---------------
ls_np_name = c("Au: Study1_12nm_0.85mg/kg","Au: Study1_23nm_0.85mg/kg","Au: Study1_100nm_0.85mg/kg",
               "Au: Study2_34.6nm_3mg/kg","Au: Study2_55.5nm_3mg/kg","Au: Study2_77.1nm_3mg/kg",
               "Au: Study2_82.6nm_3mg/kg","Au: Study3_27.6nm_4.26mg/kg","Au: Study3_27.6nm_0.85mg/kg",
               "Si: Study1_20nm_10mg/kg","Si: Study1_80nm_10mg/kg","GO: Study1_20nm_20mg/kg",
               "GO: Study2_243nm_1mg/kg", "GO: Study2_914nm_1mg/kg_w/o_CS",
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
folders = list()

for (np.name in ls_np_name) {
  print(np.name)
  folder = read_observation_data(np.name)$folder
  folders <- append(folders, paste0(folder,"mc_sens/"))
}


# Read post_mean_sd_df files from each folder into a list o data frames
data_list <- lapply(folders, function(folder) {
  read.csv(file.path(folder, "post_paras_stats.csv"))
})

# Combine data frames into a single data frame
combined_post_mean_sd_df <- bind_rows(data_list, .id = "Folder") # modify the fold id


# Convert the Folder column from numeric to factor
combined_post_mean_sd_df$Folder <- as.factor(combined_post_mean_sd_df$Folder)

# Create a named vector to map folder indices to ls_np_name
folder_names <- setNames(ls_np_name, seq_along(ls_np_name))

# Update the Folder column with names from ls_np_name
combined_post_mean_sd_df <- combined_post_mean_sd_df %>%
  mutate(Folder = folder_names[as.character(Folder)])

# Convert Folder column to factor
combined_post_mean_sd_df$Folder <- as.factor(combined_post_mean_sd_df$Folder)
combined_post_mean_sd_df$Par <- gsub("K_max","K_uptake", combined_post_mean_sd_df$Par)
print(combined_post_mean_sd_df)

# get all the parameters instead of the sig values.
subset = combined_post_mean_sd_df%>%
  filter(!grepl("^sig", Par))



# rescale the mean value as the percentage changes compared to the mean value 
# (absolute value and values)
subset.m = ddply(subset, .(Par), transform, rescale = rescale(Mean),
                 deviation_percentage = (Mean - mean(Mean)) / mean(Mean) * 100) # rescale for the each parameter
subset.m = ddply(subset.m, .(Par), transform, rescale = rescale(Mean),
                 deviation_percentage_abs = abs(Mean - mean(Mean)) / mean(Mean) * 100) # rescale for the each parameter
subset.m = ddply(subset.m, .(Par), transform, rescale = rescale(Mean),
                 deviation_abs = abs(Mean - mean(Mean))) # rescale for the each parameter


val_range <- max(abs(subset.m$deviation_percentage))
subset.m_g = subset.m %>%
  mutate(Material = case_when(
    grepl("Au", Folder) ~ "Au",
    grepl("GO", Folder) ~ "GO",
    grepl("TiO2", Folder) ~ "TiO2",
    grepl("Si", Folder) ~ "Si",
    grepl("FeO", Folder) ~ "FeO",
    TRUE ~ "Other"
  ))


subset.m_g = subset.m_g %>%
  mutate(NM = parse_number(sub(".*_(.*nm).*", "\\1", Folder))
  )


bins <- c(0, 50, 100, 500,1000)
labels <- c("0 - 50nm","50 - 100nm", "100 - 500nm", "500 - 1000nm")

# Categorize the deviation_percentage_abs values
subset.m_g <- subset.m_g %>%
  mutate(Size = cut(NM, breaks = bins, 
                        labels = labels, include.lowest = TRUE))



pars_result_folder = "/Users/mmm/work/Mouse-PBPK/plots/paras/"
#write.csv(subset.m_g,file=paste0(pars_result_folder,"pars_tot.csv"),row.names = FALSE)

pars_mean_T = subset.m_g%>%
      select(Folder, Par, Mean) %>%
     pivot_wider(names_from = Par, values_from = Mean)
#write.csv(pars_mean_T,file=paste0(pars_result_folder,"pars_T_tot.csv"),row.names = FALSE)
dataset_info <- read_excel("dataset/tk/mouse/dataset_info.xlsx")


dataset_info <- dataset_info %>%
  mutate(
    Hydrodynamic_Size = coalesce(Hydrodynamic_Size, Size),
    Size = coalesce(Size, Hydrodynamic_Size)
  )

#Nanoparticles with a zeta potential between -10 and +10 mV are considered approximately neutral, 

dataset_info <- dataset_info %>%
  mutate(
    ZP.category = cut(zeta_potential, breaks = c(-Inf, -10, 10, Inf),
                      labels = c("negative", "neutral", "positive"),
                      include.lowest = TRUE)
  )

dataset_info <- dataset_info %>%
  mutate(
    ZP.category = ifelse(is.na(ZP.category),"no_info", as.character(ZP.category))
  )


# get the mean value for each parameter
pars_mean_T_group = unique(subset(subset.m_g, select = c("Folder", "Material","Size")))
pars_mean_T_group =  as.data.frame(pars_mean_T_group)
pars_mean_T_group = merge(pars_mean_T_group, dataset_info, by.x = "Folder",by.y = "id")
row.names(pars_mean_T_group) <- pars_mean_T_group[, 1] # Set row names to the values in the first column
pars_mean_T_group <- subset(pars_mean_T_group, select = c("ZP.category","coating"))




#---calculate the dispersion of parameters-----



# calculate the range, standard deviation. variance and coefficient of variation
range_df <- subset.m_g %>%
  group_by(Par) %>%
  dplyr::summarise(max_mean = max(Mean),
                   min_mean = min(Mean),range = max(Mean)- min (Mean),
                   range_norm = (max(Mean)- min (Mean))/mean(Mean),ratio = max(Mean)/min(Mean),
                   std = sd(Mean),var = var(Mean),cov = sd(Mean)/mean(Mean))
range_df <- range_df %>%
  arrange(Par)




#-----------------------2. heatmap------------------------

pars_mean_T =  as.data.frame(pars_mean_T)

row.names(pars_mean_T) <- pars_mean_T[, 1] # Set row names to the values in the first column
pars_mean_T <- pars_mean_T[, -1]

# Create a data frame to map ls_np_name to experiment_id
mapping_df <- data.frame(ls_np_name = ls_np_name, experiment_id = experiment_id)

# Replace row names in pars_mean_T by mapping the ls_np_name to experiment_id
rownames(pars_mean_T) <- mapping_df$experiment_id[match(rownames(pars_mean_T), mapping_df$ls_np_name)]


#install.packages("pheatmap")
library(pheatmap)


annoCol <- c("#FF0000", "#92DCe5","#51CB20","#731DD8")
names(annoCol) <- unique(pars_mean_T_group$Size)


my_colour = list(
  Material = c(Au = "#FFFD82", FeO = "#697A21",GO = "#050404", 
Si = "#D5B9B2", TiO2 = "#96ADC8"),Size = annoCol
)

annoCol_coating <- c("#FF69B4", "#33FF57", "#3357FF", "#8A2BE2", "#FFD733", "#33FFF3")

# Assign names based on the unique values in the 'Size' column
names(annoCol_coating) <- unique(pars_mean_T_group$coating)

my_colour = list(
  coating = annoCol_coating
)

# Create a cluster heatmap scaled based on column
heatmap_plot <- pheatmap(pars_mean_T,
         scale = "column",  # Scale the columns，
         annotation_row  = pars_mean_T_group,
         border="black",
         cutree_rows = 6,
         #cutree_cols = 8,
         treeheight_col = 50,
         treeheight_row = 50,
         cluster_cols = FALSE,
         annotation_colors = my_colour,
         clustering_distance_rows = "canberra", # measure the distance betwee two vectors, taking into account their relative magnitudes
         clustering_method="ward.D2", # Ward D2 uses the sum of squared differences from the centroid as the criterion to minimize when merging clusters
         #clustering_distance_rows = "euclidean",  # Distance measure for rows
         clustering_distance_cols = "euclidean",  # Distance measure for columns
         #clustering_method = "average",  # Clustering method
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),  # Color palette
         display_numbers = FALSE,  # Display cell values
         fontsize = 10  # Font size
         )  # Title
heatmap_plot
#ggsave(paste0("/Users/mmm/work/Mouse-PBPK/plots/paras/",
#              "heat_map_zp_coating.png"), heatmap_plot, 
#       width = 13, height = 6)

heatmap_plot_wo_anot = pheatmap(pars_mean_T,
                         scale = "column",  # Scale the columns，
                         border="black",
                         cutree_rows = 6,
                         cutree_cols = 8,
                         treeheight_col = 50,
                         treeheight_row = 50,
                         cluster_cols = FALSE,
                         annotation_colors = my_colour,
                         clustering_distance_rows = "canberra", # measure the distance betwee two vectors, taking into account their relative magnitudes
                         clustering_method="ward.D2", # Ward D2 uses the sum of squared differences from the centroid as the criterion to minimize when merging clusters
                         clustering_distance_cols = "euclidean",  # Distance measure for columns

                         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),  # Color palette
                         display_numbers = FALSE,  # Display cell values
                         fontsize = 10  # Font size
)  # Title
#ggsave(paste0("/Users/wuji/work/code/Mouse-general-PBPK/plots/paras/",
#              "heat_map_wo_anot2.pdf"), heatmap_plot_wo_anot, 
#       width = 9, height = 6,dpi=300)


#-------------- 3. mountain & histogram plot------
combined_plots <- list()
n_points =1000
ls_multimode <- list()
for (j in sort( unique(combined_post_mean_sd_df$Par)[1:29])){
  combined_data <- data.frame()
  subset_par = subset.m_g[subset.m_g$Par == j,]
  
  for (i in 1:nrow(subset_par)) {
    data <- generate_data(subset_par$Mean[i], subset_par$SD[i],n_points)
    temp_df <- data.frame(Value = data, Label = subset_par$Folder[i], Material = subset_par$Material[i])
    combined_data <- rbind(combined_data, temp_df)
    
    # Create the ridge plot
    ordered_labels <- subset_par[order(subset_par$Folder), "Folder"]
    combined_data$Label <- factor(combined_data$Label, levels = ordered_labels)}
    library(diptest)
    if(dip.test(combined_data$Value)$p.value<0.05) {
      print(c(j,"multimodal",dip.test(combined_data$Value)$p.value))
      ls_multimode <- c(ls_multimode,j)
      }
    
  density_plot <-ggplot(combined_data, aes(x = Value, fill = Label)) +
    geom_density(alpha = 0.1, show.legend = FALSE) +
    theme_bw(base_size = 10) +
    labs(
      x = "Value",
      y = "Density") +
    theme(
          axis.text.x = element_text(angle = 90))+
    scale_x_log10()
  
  # Plot histogram
  histogram_plot <-ggplot(combined_data, aes(x = Value)) +
    geom_histogram(fill = "skyblue", color = "black", alpha = 0.7,bins=50) +
    labs(title =j, y = "Frequency") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),  # Remove x-axis title
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    scale_x_log10()
  

  library(cowplot)
  # Arrange the plots vertically
  combined_plots[[j]] <- plot_grid(histogram_plot,density_plot, ncol = 1,
                                   align = "v",rel_heights = c(0.4, 0.6))
  
  
  #ggsave(paste0("/Users/mmm/work/Mouse-PBPK/plots/paras/hist_plot/",j,
  #              "_hist_plot.png"), combined_plots[[j]], 
  #       width = 4, height = 8)
  
}
unlist(ls_multimode)
# Combine all the plots into one
final_plot <- wrap_plots(plotlist = combined_plots, ncol = 6)
final_plot



ggsave(paste0("/Users/mmm/work/Mouse-PBPK/plots/paras/hist_plot/",
              "hist_final_plot.png"), final_plot, 
       width = 30, height = 15)







