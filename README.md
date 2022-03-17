
`{r setup, include=FALSE} knitr::opts_chunk$set(echo = TRUE)`

------------------------------------------------------------------------

# Censoring-in-Network-Meta-Analysis-Time-to-Event-Data

## Description

Jansen [[1]](https://doi.org/10.1186/1471-2288-11-61) proposed the use of the family of fractional polynomial (FP)
models which offers a multidimensional treatment effect approach to
conduct a network meta-analysis of time-to-event data that does not rely
on distributional assumptions. The FP approach relies on aggregate data
from digitized survival curves divided into multiple consecutive time
intervals over the follow-up period. The observed number of event times
in each interval is then modeled as a function of the number of
participants alive at time , assuming censoring occurs before the deaths
in the consecutive follow-up intervals. Censoring can be considered a
form of the missing outcome data (MOD) problem, assuming that
participants are censored at random (CAR). CAR assumes that the outcome
is similarly distributed among those leaving and those remaining in the
study over time, conditional on fully observed covariates.


