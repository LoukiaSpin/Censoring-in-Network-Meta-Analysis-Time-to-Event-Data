#*******************************************************************************
#*                                                                                                              
#*                                                                                                              
#*  NMA of Survival Data with Fractional Polynomials to Censored participants              
#*            <Second order fractional polynomial (PMID: 21548941)>                                                                
#*                        (MISSING-AT-RANDOM ASSUMPTION)                                              
#*                                                                                                              
#*                                                                                                                                    
#*******************************************************************************



### Install necessary libraries ----
list.of.packages <- c("R2jags", "coda", "dplyr", "mcmcplots", "ggplot2", "ggrepel", "ggpubr")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages) 



## Load functions ----
source("./32_Model functions/data_preparation_function.R")
source("./32_Model functions/prepare_fppm_function.R")
source("./32_Model functions/model_fppm_function.R")
source("./32_Model functions/heterogeneity_param_prior_function.R")
source("./32_Model functions/missingness_param_prior_function.R")
source("./32_Model functions/leverage_plot_function.R")
source("./32_Model functions/fppm_plot_function.R")
source("./32_Model functions/plot_hr_surv_function.R")



## Replication of the example
# Dataset - Trial-arm information ----
MTCData <- read.table("./30_Analysis & Results/Data1.txt", header = TRUE)

# Dataset - Trial information ----
MTCAddData <- read.table("./30_Analysis & Results/Add_Data.txt", header = TRUE)


## Run random-effects NMA 
# Run FP NMA (exclusion of MOD) ----
res <- model_fppm(data_points = MTCData[, -2], 
                  data_trial = MTCAddData,
                  max_time = 60, 
                  model = "RE",
                  assumption = "IND-UNCORR",
                  heter_prior = list("halfnormal", 0, 1),
                  mean_misspar = 0,
                  var_misspar = 1,
                  P1 = -2, 
                  P2 = 1, 
                  D = 1,   # D = 1, beneficial outcome; D = 0, harmful outcome
                  n_chains = 2, 
                  n_iter = 8000, 
                  n_burnin = 2000, 
                  n_thin = 4)

# Run FP-PM NMA (Independent, uncorrelated) ----
res_pm <- model_fppm(MTCData, 
                     MTCAddData, 
                     max_time = 60,
                     model = "RE",
                     assumption = "IND-UNCORR",
                     heter_prior = list("halfnormal", 0, 1),
                     P1 = -2, 
                     P2 = 1, 
                     D = 1,   # D = 1, beneficial outcome; D = 0, harmful outcome
                     n_chains = 2, 
                     n_iter = 8000, 
                     n_burnin = 2000, 
                     n_thin = 4)



# Obtain results ----
drug_names <- c("A", "B", "C", "D")
fig <- fppm_plot(full = res_pm,
                 control = "C",
                 drug_names = drug_names,
                 time_title = "months"); fig



# Provisional: Plot hazard ratio and overall survival ----
tiff("./30_Analysis & Results/Figure HR & Survival.tiff", height = 20, width = 37, units = "cm", compression = "lzw", res = 300)
ggarrange(fig$`Hazard ratio`, fig$`Overall survival`, labels = c("A)", "B)"))
dev.off()


