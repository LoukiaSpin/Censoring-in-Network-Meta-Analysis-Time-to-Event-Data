
------------------------------------------------------------------------

# Censoring-in-Network-Meta-Analysis-Time-to-Event-Data

## Description

Jansen [[1]](https://doi.org/10.1186/1471-2288-11-61) proposed the use of the family of fractional polynomial (FP) models which offers a multidimensional treatment effect approach to conduct a network meta-analysis of time-to-event data that does not rely on distributional assumptions. The FP approach relies on aggregate data from digitized survival curves divided into multiple consecutive time intervals over the follow-up period. The observed number of event times in each interval is then modeled as a function of the number of participants alive at time, assuming censoring occurs before the deaths in the consecutive follow-up intervals. Censoring can be considered a form of the missing outcome data (MOD) problem, assuming that participants are censored at random (CAR). CAR assumes that the outcome is similarly distributed among those leaving and those remaining in the study over time, conditional on fully observed covariates.

In the case of aggregate data, statistical modeling of MOD has received considerable attention in the last years for binary and continuous outcomes,creating an elegant framework that acknowledges the uncertainty about the assumed missingness scenarios. [[2]](https://doi.org/10.1016/j.jclinepi.2018.09.002) [[3]](https://doi.org10.1002/sim.6365) This can be achieved under a model that reflects the distribution of the outcome in completers and missing participants, known as the pattern-mixture (PM) model.[[4]](https://doi.org/10.2307/2290705) Modeling MOD via the PM model offers a thorough investigation of the underlying missingness mechanisms across different studies and interventions. The PM model incorporates an informative missingness parameter that reflects belief(s) about the missingness mechanism for each arm of every study.


## Description of the repository

The repository offers the typical structure of separate folders for data, and R (code/scripts), respectively.
* The _data_ folder includes two text files: __Data1.txt__ contains the dataset in the arm-level, timepoint long format, and __Add_Data.txt__ includes the interventions for each arm of every study;
* The _R_ folder includes eight function scripts and one script to run the example (__Article example.R__) using the functions of the repository and replicate the analyses of the article.<br>

[JAGS](http://mcmc-jags.sourceforge.net/) must be installed to employ the [R2jags](https://github.com/suyusung/R2jags/issues/) package. After downloading/cloning the repo, the user can use the .Rproj file to source all code.

The next sections briefly illustrate the functions of this repository.

## Output 

Prerequisite R packages: [R2jags](https://CRAN.R-project.org/package=R2jags), [code](https://cran.r-project.org/web/packages/coda/index.html),
[mcmcplots](https://cran.r-project.org/web/packages/mcmcplots/index.html), [dplyr](https://CRAN.R-project.org/package=dplyr), [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html), and [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)

### Data preparation 

To prepare the data in the proper format to use [R2jags](https://CRAN.R-project.org/package=R2jags), we have developed the function `data.preparation()` which has the following syntax:

```r
data.preparation(data_points, data_trial)
```

#### Explaining the arguments

* data: The input is a data-frame of a one-timepoint-per-row format containing arm-level, timepoint data for each study. This format is used in Jansen [[1]](https://doi.org/10.1186/1471-2288-11-61). The columns of `data` refer to the following elements for a continuous outcome:
__time__, the observed timepoint measured in days, months, or years, and so on;
__m__, the number of missing participant outcome data in each timepoint and arm of each study. If a study does **not** report this information for any investigated arm, insert `NA` in the corresponding row(s);
__r__, the number of observed events in each timepoint and arm of each study;
__n__, the number of participants randomised in each timepoint and arm of each study;
__a__, the investigated arm (for a two-arm study, the maximum __a__ is two, for a three-arm study, the maximum __a__ is three, and so on);
__dt__, the length of the interval;
__s__, the study identifier;
__c__, the number of completers in each timepoint and arm of each study;

* measure: A character string indicating the effect measure with values `OR`, `MD`, `SMD`, or `ROM` for the odds ratio, mean difference, standardised mean difference and ratio of means, respectively.

#### Output of the function

This function returns a list with the necessary data to run [R2jags](https://CRAN.R-project.org/package=R2jags) through the function `sensitivity.analysis.mod()`.

