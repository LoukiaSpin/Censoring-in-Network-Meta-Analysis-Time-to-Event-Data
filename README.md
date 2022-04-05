
# Censoring in Network Meta-Analysis Time-to-Event Data

## Introduction

Jansen [[1]](https://doi.org/10.1186/1471-2288-11-61) proposed the use of the family of fractional polynomial (FP) models which offers a multidimensional treatment effect approach to conduct a network meta-analysis of time-to-event data that does not rely on distributional assumptions. The FP approach relies on aggregate data from digitized survival curves divided into multiple consecutive time intervals over the follow-up period. The observed number of event times in each interval is then modeled as a function of the number of participants alive at time, assuming censoring occurs before the deaths in the consecutive follow-up intervals. Censoring can be considered a form of the missing outcome data (MOD) problem, assuming that participants are censored at random (CAR). CAR assumes that the outcome is similarly distributed among those leaving and those remaining in the study over time, conditional on fully observed covariates.

In the case of aggregate data, statistical modeling of MOD has received considerable attention in the last years for binary and continuous outcomes,creating an elegant framework that acknowledges the uncertainty about the assumed missingness scenarios. [[2]](https://doi.org/10.1016/j.jclinepi.2018.09.002) [[3]](https://doi.org10.1002/sim.6365) This can be achieved under a model that reflects the distribution of the outcome in completers and missing participants, known as the pattern-mixture (PM) model.[[4]](https://doi.org/10.2307/2290705) Modeling MOD via the PM model offers a thorough investigation of the underlying missingness mechanisms across different studies and interventions. The PM model incorporates an informative missingness parameter that reflects belief(s) about the missingness mechanism for each arm of every study.

## Description of the repository

The repository offers the typical structure of separate folders for data, and R (code/scripts), respectively.
* The _data_ folder includes two text files: __Data1.txt__ contains the dataset in the arm-level, timepoint long format, and __Add_Data.txt__ includes the interventions for each arm of every study;
* The _R_ folder includes eight function scripts and one script to run the example (__Article example.R__) using the functions of the repository and replicate the analyses of the article.<br>

[JAGS](http://mcmc-jags.sourceforge.net/) must be installed to employ the [R2jags](https://github.com/suyusung/R2jags/issues/) package. After downloading/cloning the repo, the user can use the .Rproj file to source all code.

The next sections briefly illustrate the _functions_ of this repository _to perform network meta-analysis with fractional polynomial and pattern-mixture model_ and _illustrate the results_.

## Output 

Prerequisite R packages: [R2jags](https://CRAN.R-project.org/package=R2jags), [code](https://cran.r-project.org/web/packages/coda/index.html),
[mcmcplots](https://cran.r-project.org/web/packages/mcmcplots/index.html), [dplyr](https://CRAN.R-project.org/package=dplyr), [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html), and [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)

### Run the model

To perform a one-stage network meta-analysis with fractional polynomial while addressing aggregate censored participants via the pattern-mixture model, we have developed the function `model_fppm()` which has the following syntax: 

```r
model_fppm(data_points, 
           data_trial,
           max_time, 
           model,
           assumption,
           heter_prior,
           mean_misspar,
           var_misspar,
           P1, 
           P2, 
           D,   
           n_chains, 
           n_iter, 
           n_burnin, 
           n_thin)
```
#### Explaining the arguments

* data_points: A data-frame of the one-timepoint-per-row format containing arm-level data for each study. This argument is inherited by the `data_preparation()` function. 
* data_trial: A data-frame of the one-trial-per-row format with the intervention identifier per arm. This argument is also inherited by the `data_preparation()` function. 
* max_time: Positive integer for the maximum timepoint of interest 
* model: Character string indicating the analysis model with values `"RE"` for the random-effects and `"FE"` for the fixed-effect model. The default argument is `"RE"`.
* assumption: Character string indicating the structure of the informative missingness hazard ratio parameter with values `"HIE-COMMON"`, `"HIE-TRIAL"`, `"HIE-ARM"`, `"IDE-COMMON"`, `"IDE-TRIAL"`, `"IDE-ARM"`, or `"IND-UNCORR"`. The default argument is `"IDE-ARM"`. The abbreviations `"IDE"`, `"HIE"`, and `"IND"` stand for identical, hierarchical and independent, respectively, and `"UNCORR"` stands for uncorrelated.
* heter_prior: A list of three elements with the following order: 1) a character string indicating the distribution with (currently available) values `"halfnormal"`, or `"uniform"`; 2) two numeric values that refer to the parameters of the selected distribution.  See 'Details' in the function `heterogeneity_param_prior()`.
* mean_misspar: A numeric value or a vector of two numeric values for the mean of the normal distribution of the informative missingness hazard ratio parameter in the logarithmic scale . The default argument is 0 and corresponds to the censored-at-random assumption. See also 'Details' in the function `missingness_param_prior()`.
* var_misspar: A positive non-zero number for the variance of the normal distribution of the informative missingness hazard ratio parameter in the logarithmic scale. The default argument is 1. 
* P1: A scalar for the power transformation of the time variable in the first-order fractional polynomial model.
* P2 A scalar for the power transformation of the time variable in the second-order fractional polynomial model.
* D: A binary number for the direction of the outcome with values `1` for beneficial outcome and `0` for harmful outcome.
* n_chains: Positive integer specifying the number of chains for the MCMC sampling; an argument of the `jags()` function of the R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The default argument is 2.
* n_iter: Positive integer specifying the number of Markov chains for the MCMC sampling; an argument of the `jags()` function of the R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The default argument is 10000.
* n_burnin: Positive integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the `jags()` function of the R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The default argument is 1000.
* n_thin: Positive integer specifying the thinning rate for the MCMC sampling; an argument of the `jags()` function of the R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The default argument is 1.

For details about the arguments, read the __documentation__ of the function.

The model runs in `JAGS` and the progress of the simulation appears on the R console. The output of `model_fppm()` is used as an S3 object by the functions `plot_hr_surv()`, `leverage_plot()`, and `fppm_plot()` to be processed further and provide an end-user-ready output.

#### Using the example 

```r
res <- model_fppm(data_points = MTCData, 
                  data_trial = MTCAddData, 
                  max_time = 60,
                  model = "RE",
                  assumption = "IND-UNCORR",
                  heter_prior = list("halfnormal", 0, 1),
                  P1 = -2, 
                  P2 = 1, 
                  D = 1,   
                  n_chains = 2, 
                  n_iter = 80000, 
                  n_burnin = 20000, 
                  n_thin = 4)
```

### Plot the Hazard Ratio and Survival Probability 

To plot the hazard ratio over time for comparisons with the selected intervention of the network as well as the survival probability for all interventions at once, use the function `plot_hr_surv()` which has the following syntax: 

```r
plot_hr_surv(full, 
             type, 
             drug_names, 
             control, 
             time_title)
```
#### Explaining the arguments

* full: An object of S3 class `model_fppm()`. See 'Value' in function `model_fppm()`.
* type: Character string indicating the plot with values `"HR"` for the hazard ratio plot and `"Surv"` for the survival probability. 
* drug_names: A vector of labels with the name of the interventions in the order they appear in the argument `data_trial` of `model_fppm()`.
* control: A character to indicate the comparator intervention. It must be any name found in `drug_names`.
* time_title: A title to indicate the time measure (e.g., days, weeks, months, and so on).

#### Using the example 

```r
# The names of the interventions in the order they appear in the dataset
interv_names <- c("Docetaxel", "Best Supportive Care", "Pemetrexed", "Gefitinib")

plot_hr_surv(full = res,
             type = "HR",
             drug_names = interv_names,
             control = "Docetaxel",
             time_title = "Time in months")
```
### The leverage plot

To plots the leverage against the square root of the posterior mean of residual deviance of the trial-arms under the model of interest, we have developed the function `leverage_plot()` which has the following syntax: 

```r
leverage_plot(full, title) 
```
#### Explaining the arguments

* full: An object of S3 class `model_fppm()`. See 'Value' in function `model_fppm()`.
* title: A title to indicate the model (e.g., fractional polynomial with pattern-mixture model, or fractional polynomial with excluded missing participants).

#### Using the example

```r
leverage_plot(full = res,
              time = "Fractional polynomial with pattern-mixture model")
```

### End-user-ready results 

To plot all results (hazard ratio, survival probability, leverage plot and tabulated estimates) at once, use the function `fppm_plot()` which has the following syntax: 

```r
fppm_plot(full, 
          drug_names, 
          control, 
          time_title,
          save_xls)
```
#### Explaining the arguments 

* full: An object of S3 class `model_fppm()`. See 'Value' in function `model_fppm()`.
* drug_names: A vector of labels with the name of the interventions in the order they appear in the argument `data_trial` of `model_fppm()`.
* control: A character to indicate the comparator intervention. It must be any name found in `drug_names`.
* time_title: A title to indicate the time measure (e.g., days, weeks, months, and so on).
* save_xls: Logical to indicate whether to export the tabulated results to an _xlsx_ file (via the `write_xlsx()` function of the R-package [writexl](https://CRAN.R-project.org/package=writexl) to the working directory of the user. The default is `FALSE` (do not export).

#### Using the example

```r
fppm_plot(full = res,
          control = "Docetaxel",
          drug_names = interv_names,
          time_title = "Time in months")
```
