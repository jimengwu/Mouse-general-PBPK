ls_np_name = c("Au: Study1_12nm_0.85mg/kg","Au: Study1_23nm_0.85mg/kg","Au: Study1_100nm_0.85mg/kg",
               "Au: Study2_34.6nm_3mg/kg","Au: Study2_55.5nm_3mg/kg","Au: Study2_77.1nm_3mg/kg",
               "Au: Study2_82.6nm_3mg/kg","Au: Study3_27.6nm_4.26mg/kg","Au: Study3_27.6nm_0.85mg/kg",
               "Si: Study1_20nm_10mg/kg","Si: Study1_80nm_10mg/kg","GO: Study1_20nm_20mg/kg",
               "GO: Study2_243nm_1mg/kg", "GO: Study2_914nm_1mg/kg_all","GO: Study2_914nm_1mg/kg_w/o_CS",
               "TiO2: Study1_385nm_10mg/kg","TiO2: Study2_220nm_60mg/kg",
               "FeO: Study1_29nm_5mg/kg","FeO: Study2_41nm_4mg/kg")

read_observation_data <- function(np_name) {
  if (np_name == "Au: Study1_12nm_0.85mg/kg") {
    Obs.df <- read.csv(file = "dataset/tk/mouse/gold/R_input_mouse_Study1_12nm.csv")
    Obs.df <- Obs.df[1:7,]
    PDOSE <- Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder <- 'plots/Au/S1_12nm_nsc_GI/'
  } 
  else if (np_name == "Au: Study1_23nm_0.85mg/kg") {
    Obs.df <- read.csv(file = "dataset/tk/mouse/gold/R_input_mouse_Study1_23nm.csv")
    Obs.df <- Obs.df[1:7,]
    PDOSE <- Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder <- 'plots/Au/S1_23nm_nsc_GI/'
  } 
  else if (np_name == "Au: Study1_100nm_0.85mg/kg") {
    Obs.df <- read.csv(file = "dataset/tk/mouse/gold/R_input_mouse_study1_100nm.csv")
    Obs.df <- Obs.df[1:7,]
    PDOSE <- Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","Clung")] #100nm only three organs, kidney only have two data points
    folder <- 'plots/Au/S1_100nm_nsc_GI/'
  } 
  else if (np_name == "Au: Study2_34.6nm_3mg/kg"){
    #--------------------------Study 2: 34.6nm------------------------------
    Obs.df <- read.csv(file ="dataset/tk/mouse/gold/Study2_PEG_coated_AuNPs_34.6nm.csv") 
    Obs.df <- Obs.df[1:6,]
    pathway <- Obs.df$Input.pathway[1]
    PDOSE = Obs.df$Dose.mg.kg.[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder = 'plots/Au/S2_34.6_nsc_GI/'
  }
  else if (np_name == "Au: Study2_55.5nm_3mg/kg"){
    #--------------------------Study 3: 55.5nm------------------------------
    Obs.df <- read.csv(file ="dataset/tk/mouse/gold/Study2_PEG_coated_AuNPs_55.5nm.csv") 
    Obs.df <- Obs.df[1:8,]
    pathway <- Obs.df$Input.pathway[1]
    PDOSE = Obs.df$Dose..mg.kg.[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder = 'plots/Au/S2_55.5_nsc_GI/'}
  else if (np_name == "Au: Study2_77.1nm_3mg/kg"){
    #--------------------------Study 3: 77.1nm------------------------------
    Obs.df <- read.csv(file ="dataset/tk/mouse/gold/Study2_PEG_coated_AuNPs_77.1nm.csv") 
    Obs.df <- Obs.df[1:8,]
    PDOSE = Obs.df$Dose..mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder = 'plots/Au/S2_77.1_nsc_GI/'}
  else if (np_name == "Au: Study2_82.6nm_3mg/kg"){
    #--------------------------Study 3: 82.6nm------------------------------
    Obs.df <- read.csv(file ="dataset/tk/mouse/gold/Study2_PEG_coated_AuNPs_82.6nm.csv") 
    Obs.df <- Obs.df[1:8,]
    PDOSE = Obs.df$Dose..mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder = 'plots/Au/S2_82.6_nsc_GI/'
    }
  else if (np_name == "Au: Study3_27.6nm_4.26mg/kg"){
    #------------------Study3: 27.6 nm high dose-------------------------
    Obs.df <- read.csv(file ="dataset/tk/mouse/gold/Study3_Gold dextran_h.csv") 
    Obs.df <- Obs.df[1:5,]
    pathway <- Obs.df$Input.pathway[1]
    PDOSE = Obs.df$Dose.mg.kg.[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder = 'plots/Au/S3_27.6nm_Au_h/'
    Obs.df$Time[1] =5/60
    
  }
  else if (np_name == "Au: Study3_27.6nm_0.85mg/kg"){
    #------------------Study3: 27.6 nm low dose-------------------------
    Obs.df <- read.csv(file ="dataset/tk/mouse/gold/Study3_Gold dextran_l.csv") 
    Obs.df <- Obs.df[1:5,]
    pathway <- Obs.df$Input.pathway[1]
    PDOSE = Obs.df$Dose.mg.kg.[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder = 'plots/Au/S3_27.6nm_Au_l/'
    Obs.df$Time[1] =5/60
    
  }
  else if (np_name == "Si: Study1_20nm_10mg/kg") {
    Obs.df <- read.csv(file = "dataset/tk/mouse/silica/Study1_125I-SiNPs 20nm.csv")
    Obs.df <- Obs.df[1:5,]
    PDOSE <- Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder <- 'plots/silica/S1_20nm/'
  } 
  else if (np_name == "Si: Study1_80nm_10mg/kg") {
    Obs.df <- read.csv(file = "dataset/tk/mouse/silica/Study1_125I-SiNPs 80nm.csv")
    Obs.df <- Obs.df[1:5,]
    PDOSE <- Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder <- 'plots/silica/S1_80nm/'
  } 
  else if (np_name == "GO: Study1_20nm_20mg/kg") {
    Obs.df <- read.csv(file = "dataset/tk/mouse/GO/Study1_125I-NGS-PEG 10-30nm.csv")
    Obs.df <- Obs.df[1:5,]
    PDOSE <- Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder <- 'plots/GO/S1_20nm/'
  } 
  else if (np_name == "GO: Study2_243nm_1mg/kg") {
    Obs.df <- read.csv(file = "dataset/tk/mouse/GO/Study2_125I-s-GO 243nm.csv")
    Obs.df <- Obs.df[1:5,]
    PDOSE <- Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder <- 'plots/GO/S2_243nm/'
    Obs.df$Time[1] = 2/60
    
    Obs.df$Time[2] = 5/60
    
    Obs.df$Time[3] = 10/60
  } 
  else if (np_name == "GO: Study2_914nm_1mg/kg_all") {
    Obs.df <- read.csv(file = "dataset/tk/mouse/GO/Study2_125I-l-GO 914nm.csv")
    Obs.df <- Obs.df[1:5,]
    PDOSE <- Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder <- 'plots/GO/S2_914nm/'
    Obs.df$Time[1] = 2/60
    
    Obs.df$Time[2] = 5/60
    
    Obs.df$Time[3] = 10/60
  } 
  else if (np_name == "GO: Study2_914nm_1mg/kg_w/o_CS") {
    Obs.df <- read.csv(file = "dataset/tk/mouse/GO/Study2_125I-l-GO 914nm.csv")
    Obs.df <- Obs.df[1:5,]
    PDOSE <- Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CK","Clung")]
    folder <- 'plots/GO/S2_914nm_wo_CS/'
    Obs.df$Time[1] = 2/60
    
    Obs.df$Time[2] = 5/60
    
    Obs.df$Time[3] = 10/60
  } 
  else if (np_name == "TiO2: Study1_385nm_10mg/kg") {
    #------------------TiO2: Study1: 385 nm-------------------------
    Obs.df <- read.csv(file ="dataset/tk/mouse/TiO2/Study1-TiO2.csv") 
    Obs.df <- Obs.df[1:5,]
    PDOSE = Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder = 'plots/TiO2/S1_385nm_TiO2/'
    
  } 
  else if (np_name == "TiO2: Study2_220nm_60mg/kg") {
    #------------------TiO2: Study1: 100 nm-------------------------
    Obs.df <- read.csv(file ="dataset/tk/mouse/TiO2/Study2-TiO2.csv") 
    Obs.df = Obs.df[1:3,]
    PDOSE = Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","Clung","CK")]
    folder = 'plots/TiO2/S2_220nm_TiO2_60/'
    Obs.df$Time[1] = 5/60
  } 
  else if (np_name == "FeO: Study1_29nm_5mg/kg") {
    #------------------IO: Study1: 29 nm-------------------------
    Obs.df <- read.csv(file ="dataset/tk/mouse/iron/Study2_EDT_IONP_29nm.csv") 
    Obs.df <- Obs.df[1:4,]
    PDOSE = Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder = 'plots/iron/S1_29nm_IO/'
    Obs.df$Time[1] = 1/60
    
    Obs.df$Time[2] = 2/60
    
    Obs.df$Time[3] = 15/60
    
  } 
  else if (np_name == "FeO: Study2_41nm_4mg/kg") {
    #------------------IO: Study1: 29 nm-------------------------
    Obs.df <- read.csv(file ="dataset/tk/mouse/iron/Study3_USPIO_41nm.csv") 
    Obs.df <- Obs.df[1:12,]
    PDOSE = Obs.df$Dose.mg.kg.[1]
    pathway <- Obs.df$Input.pathway[1]
    Obs.df <- Obs.df[c("Time","CL","CS","CK","Clung")]
    folder = 'plots/iron/S3_41nm_IO/'
    Obs.df$Time[1] = 1/60
    
    Obs.df$Time[2] = 2.5/60
    
    Obs.df$Time[3] = 5/60
    
  } 
  else {
    # Handle other cases or provide a default behavior
    Obs.df <- NULL
    folder <- NULL
  }
  
  return(list(Obs.df = Obs.df, folder = folder,pathway= pathway,PDOSE=PDOSE))
}
