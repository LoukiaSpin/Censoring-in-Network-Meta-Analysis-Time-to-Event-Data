
------------------------------------------------------------------------

# Censoring-in-Network-Meta-Analysis-Time-to-Event-Data

## Description

Jansen [[1]](https://doi.org/10.1186/1471-2288-11-61) proposed the use of 
the family of fractional polynomial (FP) models which offers a 
multidimensional treatment effect approach to
conduct a network meta-analysis of time-to-event data that does not rely
on distributional assumptions. The FP approach relies on aggregate data
from digitized survival curves divided into multiple consecutive time
intervals over the follow-up period. The observed number of event times
in each interval is then modeled as a function of the number of
participants alive at time, assuming censoring occurs before the deaths
in the consecutive follow-up intervals. Censoring can be considered a
form of the missing outcome data (MOD) problem, assuming that
participants are censored at random (CAR). CAR assumes that the outcome
is similarly distributed among those leaving and those remaining in the
study over time, conditional on fully observed covariates.

In the case of aggregate data, statistical modeling of MOD has received considerable attention in the last years for binary and continuous outcomes,creating an elegant framework that acknowledges the uncertainty about the assumed missingness scenarios. [[2]](https://doi.org/10.1016/j.jclinepi.2018.09.002) [[3]](https://doi.org10.1002/sim.6365) This can be achieved under a model that reflects the distribution of the outcome in completers and missing participants, known as the pattern-mixture (PM) model.[[4]](https://doi.org/10.2307/2290705) Modeling MOD via the PM model offers a thorough investigation of the underlying missingness mechanisms across different studies and interventions. The PM model incorporates an informative missingness parameter that reflects belief(s) about the missingness mechanism for each arm of every study.


## Description of the repository

The repository offers the typical structure of separate folders for data, and R (code/scripts), respectively.
* The _data_ folder includes two text files: __Add_Data.txt__ contains the dataset in the arm-level timepoint wide format. Finally, the text files contain the posterior summaries from the Bayesian analysis (via the [R2jags](https://github.com/suyusung/R2jags/issues/) package) under all pre-defined scenarios for the missingness mechanism;
* The _R_ folder includes two analysis scripts (__A.Run empirical analysis_PMA & NMA.R__ and __Β.Determine robustness_PMA & NMA.R__). The first script sources the RData files with the datasets and performs the Bayesian analyses (one-stage random-effects PMA and NMA for continuous and binary outcomes with incorporation of the pattern-mixture model). The second script calculates the robustness index and apply our proposed decision framework for robustness of the primary analysis results for each analysis. Furthermore, it produces 1) the heatmap of robustness index for all possible comparisons of the network and 2) the panel of density plots of the summary log OR under all missingness scenarios and the Kullback-Leibler divergence. The network of Liu et al. (the motivating example in our article) has been used for that purpose. Finally, the remaining R scripts refer to the necessary functions to apply the aforementioned analysis scripts.<br>

[JAGS](http://mcmc-jags.sourceforge.net/) must be installed to employ the [R2jags](https://github.com/suyusung/R2jags/issues/) package. After downloading/cloning the repo, the user can use the .Rproj file to source all code.

The next sections briefly illustrate the functions of this repository.


