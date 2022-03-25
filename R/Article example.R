#*******************************************************************************
#*                                                                                                              
#*                                                                                                              
#*  NMA of Survival Data with Fractional Polynomials to Censored participants              
#*            <Second order fractional polynomial (PMID: 21548941)>                                                                
#*                        (MISSING-AT-RANDOM ASSUMPTION)                                              
#*                                                                                                              
#*                                                                                                                                    
#*******************************************************************************



## Install necessary libraries ----
list.of.packages <- c("R2jags", "coda", "dplyr", "mcmcplots", "ggplot2", "ggrepel", "ggpubr")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages) 

options(mc.cores = parallel::detectCores())
set.seed(12345)

## Load functions ----
source("./R/data_preparation_function.R")
source("./R/prepare_fppm_function.R")
source("./R/model_fppm_function.R")
source("./R/heterogeneity_param_prior_function.R")
source("./R/missingness_param_prior_function.R")
source("./R/leverage_plot_function.R")
source("./R/fppm_plot_function.R")
source("./R/plot_hr_surv_function.R")



## Replication of the example
# Dataset - Trial-arm information ----
MTCData <- read.table("./data/Data1.txt", header = TRUE)

# Dataset - Trial information ----
MTCAddData <- read.table("./data/Add_Data.txt", header = TRUE)



# Run random-effects FP NMA (exclusion of MOD) ----
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
                  n_iter = 80000, 
                  n_burnin = 20000, 
                  n_thin = 4)

# Run random-effects FP-PM NMA (Independent, uncorrelated) ----
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
                     n_iter = 80000, 
                     n_burnin = 20000, 
                     n_thin = 4)



# Obtain results ----
drug_names <- c("Docetaxel", "Best Supportive Care", "Pemetrexed", "Gefitinib")
fig <- fppm_plot(full = res,
                 control = "Docetaxel",
                 drug_names = drug_names,
                 time_title = "Time (months)"); fig


