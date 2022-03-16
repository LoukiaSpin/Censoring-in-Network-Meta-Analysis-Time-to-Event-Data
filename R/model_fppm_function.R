#' Perform Bayesian network meta-analysis with fractional polynomials
#'
#' @description
#'   Performs a one-stage network meta-analysis with second-order fractional 
#'   polynomials for aggregate time-to-event data while addressing aggregate 
#'   (binary) missing participant outcome data via the pattern-mixture model.
#'
#' @param data_points A data-frame in long format with trial-arm data for each 
#'   time-point. See 'Format' below.
#' @param data_trial A data-frame of the one-trial-per-row format with arm-level
#'   data interventions studied in each arm of every trial.
#'   In essence, an intervention identifier in each arm.
#' @param model Character string indicating the analysis model with values
#'   \code{"RE"}, or \code{"FE"} for the random-effects and fixed-effect model,
#'   respectively. The default argument is \code{"RE"}.
#' @param assumption Character string indicating the structure of the
#'   informative missingness hazard ratio parameter (submitted).
#'   Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"},
#'   \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"},
#'   \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, or \code{"IND-UNCORR"}.
#'   The default argument is \code{"IDE-ARM"}. The abbreviations \code{"IDE"},
#'   \code{"HIE"}, and \code{"IND"} stand for identical, hierarchical and
#'   independent, respectively. and \code{"UNCORR"} stands for uncorrelated.
#' @param heter_prior A list of three elements with the following order:
#'   1) a character string indicating the distribution with
#'   (currently available) values \code{"halfnormal"}, or \code{"uniform"}; 
#'   2) two numeric values that refer to the parameters of the selected 
#'   distribution. For \code{"halfnormal"}, these numbers refer to zero and the 
#'   scale parameter (equal to 4 or 1 being the corresponding precision of the 
#'   scale parameter 0.5 or 1). For \code{"uniform"}, these numbers refer to the
#'   minimum and maximum value of the distribution.
#'   See 'Details' in \code{\link{heterogeneity_param_prior}}.
#' @param mean_misspar A numeric value or a vector of two numeric values for the
#'   mean of the normal distribution of the informative missingness hazard ratio 
#'   parameter in the logarithmic scale (see 'Details'). The default argument is
#'   0 and corresponds to the missing-at-random assumption.
#'   See also 'Details' in \code{\link{missingness_param_prior}}.
#' @param var_misspar A positive non-zero number for the variance of the
#'   normal distribution of the informative missingness hazard ratio parameter 
#'   in the logarithmic scale.
#'   When the \code{measure} is \code{"OR"}, \code{"MD"}, or \code{"SMD"}
#'   the default argument is 1. When the \code{measure} is \code{"ROM"}
#'   the default argument is 0.04.
#' @param P1 A scalar chosen from a set within {-2, -1, -0.5, 0, 0.5, 1, 2, 3} 
#'   for the power transformation of the time variable in the first-order 
#'   fractional polynomial model (Jansen, 2011; Royston and Altman, 1994).
#' @param P2 A scalar chosen from a set within {-2, -1, -0.5, 0, 0.5, 1, 2, 3}
#'   for the power transformation of the time variable in the second-order 
#'   fractional polynomial model (Jansen, 2011; Royston and Altman, 1994).
#' @param D A binary number for the direction of the outcome.
#'   Set \code{D = 1} for beneficial outcome and \code{D = 0} for harmful
#'   outcome.
#' @param n_chains Positive integer specifying the number of chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n_iter Positive integer specifying the number of Markov chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n_burnin Positive integer specifying the number of iterations to
#'   discard at the beginning of the MCMC sampling; an argument of the
#'   \code{\link[R2jags:jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n_thin Positive integer specifying the thinning rate for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @format The data-frame \code{data_points} includes the following columns:
#'   \item{time}{The measured time-points (they may differ across the 
#'   trials).}
#'   \item{r}{The number of observed events of the outcome at each
#'   time-point and trial-arm.}
#'   \item{m}{The number of missing participant outcome data at each
#'   time-point and trial-arm.}
#'   \item{c}{The number of completers at each time-point and 
#'   trial-arm.}
#'   \item{n}{The number of randomised participants at each time-point and
#'   trial-arm.}
#'   \item{a}{The intervention arm at each time-point and trial.}
#'   \item{s}{The trial identifier.}
#'   \item{dt}{The length of time-interval at each trial-arm.}
#'   The number of rows of \code{data_points} equals the sum of time-points per 
#'   arm across the trials.
#'
#'   In \code{data_trial}, the intervention identifier needs to be in an 
#'   ascending order across the arms of each trial. \code{model_fppm} considers 
#'   the first column in \strong{t} as being the control arm for every trial. 
#'   Thus, the interventions with a lower identifier are consistently treated as
#'   the control arm in each trial. This case is relevant in non-star-shaped 
#'   networks. The number of rows in \code{data_trial} equals the number of 
#'   collected trials.\strong{t} appears in \code{data_trial} as many times as 
#'   the maximum number of interventions compared in a trial of the dataset.
#'   In network meta-analysis without multi-arm trials, the maximum number of 
#'   arms is inherently two. In the case of network meta-analysis with multi-arm
#'   trials, the maximum number of arms exceeds two. 
#'
#' @return A list of R2jags output on the summaries of the posterior
#'   distribution, and the Gelman-Rubin convergence diagnostic
#'   (Gelman et al., 1992) of the following monitored parameters for a
#'   fixed-effect network meta-analysis:
#'   \item{EM}{The estimated summary hazard ratio in the logarithmic scale for 
#'   all possible pairwise comparisons of interventions in the network.}
#'   \item{EM_ref}{The estimated summary hazard ratio in the logarithmic scale 
#'   for all comparisons with the reference intervention in the network.}
#'   \item{log_HR}{The estimated summary hazard ratio in the logarithmic scale 
#'   for all possible pairwise comparisons of interventions in the network and 
#'   for time-points in the interval [1, max] with max being the maximum 
#'   time-point observed in the dataset.}
#'   \item{Surv}{The survival probability for all interventions in the network 
#'   and for time-points in the interval [1, max] with max being the maximum 
#'   time-point observed in the dataset.}
#'   \item{expect_Surv}{The expected survival time for all interventions in the 
#'   network.}
#'   \item{dev_o}{The deviance contribution of each trial-arm and time-point 
#'   based on the observed outcome.}
#'   \item{hat_par}{The fitted outcome at each trial-arm and time-point.}
#'   \item{phi}{The informative missingness hazard ratio parameter in the
#'   logarithmic scale.}
#'   \item{SUCRA1}{The surface under the cumulative ranking curve for each
#'   intervention under the proportional hazards model.}
#'   \item{SUCRA2}{The surface under the cumulative ranking curve for each
#'   intervention under the first-order fractional polynomial.}
#'   \item{SUCRA3}{The surface under the cumulative ranking curve for each
#'   intervention under the second-order fractional polynomial.}
#'   \item{effectiveneness1}{The ranking probability of each intervention for
#'   every rank under the proportional hazards model.}
#'   \item{effectiveneness2}{The ranking probability of each intervention for
#'   every rank under the first-order fractional polynomial.}
#'   \item{effectiveneness3}{The ranking probability of each intervention for
#'   every rank under the second-order fractional polynomial.}
#'
#'   For a random-effects network meta-analysis, the output additionally
#'   includes the following elements:
#'   \item{EM_pred}{The predicted summary hazard ratio in the logarithmic scale 
#'   for all possible pairwise comparisons of interventions in the network.}
#'   \item{pred_ref}{The predicted summary hazard ratio in the logarithmic 
#'   scale for all comparisons with the reference intervention in the network.}
#'   \item{delta}{The estimated trial-specific effect measure (according to the
#'   argument \code{measure}).}
#'   \item{tau}{The between-trial standard deviation.}
#'   \code{tau} is typically assumed to be common for all observed comparisons
#'   in the network. For a multi-arm trial, we estimate a total of \emph{T-1}
#'   \code{delta} for comparisons with the baseline intervention of the trial
#'   (found in the first column of the element \bold{t} in \code{data_trial}), 
#'   with \emph{T} being the number of interventions in the trial.
#'
#'   Furthermore, the output includes the following elements:
#'   \item{lev_o}{The leverage for the observed outcome at each trial-arm
#'   and time-point.}
#'   \item{sign_o}{The sign of the difference between observed and fitted
#'   outcome at each trial-arm and time-point.}
#'   \item{model_assessment}{A data-frame on the measures of model assessment:
#'   deviance information criterion, number of effective parameters, and total
#'   residual deviance.}
#'   \item{jagsfit}{An object of S3 class \code{\link[R2jags:jags]{jags}} with
#'   the posterior results on all monitored parameters to be used in the
#'   \code{\link{mcmc_diagnostics}} function.}
#'
#'   The \code{model_fppm} function also returns the arguments \code{model}, 
#'   \code{assumption}, and \code{D}, as specified by the user to be considered 
#'   in other functions.
#'
#' @details The model runs in \code{JAGS} and the progress of the simulation
#'   appears on the R console. The output of \code{model_fppm} is used as an S3
#'   object by other functions to be processed further and provide an 
#'   end-user-ready output.
#'
#'   The \code{\link{data_preparation}} function is called to prepare the data
#'   for the Bayesian analysis. \code{\link{data_preparation}} checks whether
#'   the element \strong{m} exists in the \code{data_points}. If this element is
#'   missing, \code{\link{data_preparation}} creates 1) a pseudo-data-frame for
#'   \strong{m} that has the zero value for the observed trial-arms, and
#'   \code{NA} for the unobserved trial-arms, and 2) the pseudo-data-frame
#'   \code{I} that is identical with the pseudo-data-frame for \code{m}.
#'   If the element \strong{m} exists in the \code{data_points} and has values 
#'   only for some trial-arms, the pseudo-data-frame for \strong{m} is identical
#'   to \strong{m} for the corresponding trial-arms, and the pseudo-data-frame
#'   \code{I} has the value one for these trial-arms. Both pseudo-data-frames
#'   aim to retain the trials without information on missing participant outcome
#'   data.
#'
#'   To perform a Bayesian network meta-analysis, the 
#'   \code{\link{prepare_fppml}} function is called which contains the WinBUGS
#'   code as written by Jansen (2011) for aggregate time-to-event data.
#'   \code{\link{prepare_fppml}} uses the consistency model (as described in
#'   Lu and Ades (2006)) to estimate all possible comparisons in the network.
#'   It also accounts for the multi-arm trials by assigning conditional
#'   univariate normal distributions on the basic parameters of these trials,
#'   namely, effect parameters between the non-baseline arms and the baseline
#'   arm of the multi-arm trial (Dias et al., 2013).
#'
#'   The code of Jansen (2011) has been extended to incorporate the
#'   pattern-mixture model to adjust the underlying outcome in each arm of
#'   every trial for aggregate missing participant outcome data (Spineli, 2019a; 
#'   Turner et al., 2015). The assumptions about the missingness parameter are 
#'   specified using the arguments \code{mean_misspar} and \code{var_misspar}. 
#'   Specifically, \code{model_fppm} considers the informative missingness 
#'   hazard ratio in the logarithmic scale (submitted).
#'
#'   When \code{assumption} is trial-specific (i.e., \code{"IDE-TRIAL"} or
#'   \code{"HIE-TRIAL"}), or independent (\code{"IND-UNCORR"}), only one numeric
#'   value can be assigned to \code{mean_misspar} because the same missingness 
#'   scenario is applied to all trials and trial-arms of the dataset, 
#'   respectively. When \code{assumption} is \code{"IDE-ARM"} or 
#'   \code{"HIE-ARM"}, a maximum of two \emph{different} or \emph{identical} 
#'   numeric values can be assigned as a vector to \code{mean_misspars}: the 
#'   first value refers to the experimental arm, and the second value refers to 
#'   the control arm of a trial. Specifically, the first value is considered for
#'   all non-reference interventions and the second value is considered for the
#'   reference intervention of the network (i.e., the intervention with
#'   identifier equal to one). This is necessary to ensure transitivity
#'   in the assumptions for the missingness parameter across the network
#'   (Spineli, 2019b).
#'
#' @author {Loukia M. Spineli}, {Chrysostomos Kalyvas}
#'
#' @seealso \code{\link{data_preparation}},
#'   \code{\link{heterogeneity_param_prior}}, \code{\link[R2jags:jags]{jags}},
#'   \code{\link{missingness_param_prior}}, \code{\link{prepare_fppm}}
#'
#' @references
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
#' making 2: a generalized linear modeling framework for pairwise and network
#' meta-analysis of randomized controlled trials. \emph{Med Decis Making}
#' 2013;\bold{33}(5):607--617. \doi{10.1177/0272989X12458724}
#'
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}:457--472.
#'
#' Jansen JP. Network meta-analysis of survival data with fractional 
#' polynomials. \emph{BMC Med Res Methodol} 2011;\bold{11}:61.
#' \doi{10.1186/1471-2288-11-61}
#'
#' Kalyvas C, Papadimitropoulou K, Malbecq W, Spineli LM. Dealing with censoring
#' in a network meta-analysis of time-to-event data. 2022 (submitted)
#'
#' Lu G, Ades AE. Assessing evidence inconsistency in mixed treatment
#' comparisons. \emph{J Am Stat Assoc} 2006;\bold{101}:447--459.
#' \doi{10.1198/016214505000001302}
#'
#' Royston P, Altman DG. Regression Using Fractional Polynomials of 
#' Continuous Covariates: Parsimonious Parametric Modelling. \emph{Appl Stat} 
#' 1994;\bold{43}(3):429--453.
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019a;\bold{19}(1):86.
#' \doi{10.1186/s12874-019-0731-y}
#'
#' Spineli LM. Modeling missing binary outcome data while preserving
#' transitivity assumption yielded more credible network meta-analysis
#' results. \emph{J Clin Epidemiol} 2019b;\bold{105}:19--26.
#'
#' Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account
#' for uncertainty due to missing binary outcome data in pairwise
#' meta-analysis. \emph{Stat Med} 2015;\bold{34}(12):2062--2080.
#' \doi{10.1002/sim.6475}
#'
model_fppm <- function(data_points, 
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
                       n_thin) {
  
  item <- data_preparation(data_points, data_trial)
  
  # Missing and default arguments
  max_time <- if (missing(max_time)) {
    stop("The argument 'max_time' needs to be defined", call. = FALSE)
  } else {
    max_time
  }
  model <- if (missing(model)) {
    "RE"
  } else if (!is.element(model, c("RE", "FE"))) {
    stop("Insert 'RE', or 'FE'", call. = FALSE)
  } else {
    model
  }
  assumption <- if (missing(assumption)) {
    message("The 'IDE-ARM' has been used as the default.")
    "IDE-ARM"
  } else if (min(item$I) == 0 & max(item$I) == 1 & assumption != "IND-UNCORR") {
    aa <- "Missing participant outcome data have been collected partially."
    bb <- "Insert 'IND-UNCORR'."
    stop(paste(aa, bb), call. = FALSE)
  } else {
    assumption
  } 
  heterog_prior <- heterogeneity_param_prior(model, heter_prior)
  mean_misspar <- missingness_param_prior(assumption, mean_misspar)
  var_misspar <- if (missing(var_misspar)) {
    1
  } else {
    var_misspar
  }
  P1 <- if (missing(P1)) {
    stop("The argument 'P1' needs to be defined", call. = FALSE)
  } else {
    P1
  }
  P2 <- if (missing(P2)) {
    stop("The argument 'P2' needs to be defined", call. = FALSE)
  } else {
    P2
  }
  D <- if (missing(D)) {
    stop("The argument 'D' needs to be defined", call. = FALSE)
  } else {
    D
  }
  n_chains <- ifelse(missing(n_chains), 2, n_chains)
  n_iter <- ifelse(missing(n_iter), 10000, n_iter)
  n_burnin <- ifelse(missing(n_burnin), 1000, n_burnin)
  n_thin <- ifelse(missing(n_thin), 1, n_thin)

  # Data in list format for R2jags
  data_jag <- list("dpoints" = item$dpoints,
                   "n" = item$n,
                   "r" = item$r, 
                   "m" = item$mod, 
                   "c" = item$c, 
                   "a" = item$a, 
                   "s" = item$s, 
                   "time" = item$time,
                   "I" = item$I,
                   "dt" = item$dt, 
                   "P1" = P1, 
                   "P2" = P2, 
                   "t" = item$t, 
                   "ns" = item$ns,
                   "nt" = item$nt, 
                   "ref" = 1, 
                   "D" = D, 
                   "na" = item$na, 
                   "maxt" = max_time,
                   "meand_phi" = mean_misspar,
                   "precd_phi" = 1 / var_misspar)
  
  data_jag <- if (model == "RE") {
    append(data_jag, list("heter_prior" = heterog_prior))
  } else {
    data_jag
  }
  
  param_jags <- c("Beta",
                  "Surv",
                  "expect_Surv", 
                  "EM", 
                  "EM_ref",
                  "log_HR",
                  "SUCRA1", 
                  "SUCRA2", 
                  "SUCRA3", 
                  "rhat_o", 
                  "totresdev_o",
                  "dev_o")
  
  param_jags <- if (is.element(assumption,
                               c("HIE-COMMON", "HIE-TRIAL", "HIE-ARM"))) {
    append(param_jags, "mean_phi")
  } else {
    append(param_jags, "phi")
  }
  
  param_jags <- if (model == "FE") {
    param_jags
  } else {
    append(param_jags, c("EM_pred", "pred_ref", "tau", "delta"))
  }
  
  ## Run model 
  jagsfit <- jags(data = data_jag, 
                  parameters.to.save = param_jags, 
                  model.file = textConnection(
                    prepare_fppm(model,
                                 assumption)
                    ),
                  n.chains = n_chains, 
                  n.iter = n_iter, 
                  n.burnin = n_burnin, 
                  n.thin = n_thin, 
                  DIC = F, 
                  digits = 4, 
                  working.directory = getwd())
  
  # Turn R2jags object into a data-frame
  get_results <- as.data.frame(t(jagsfit$BUGSoutput$summary))
  
  # Log hazard per arm and order
  Beta <- t(get_results %>% dplyr::select(starts_with("Beta[")))
  
  # Within-trial effects size
  delta <- t(get_results %>% dplyr::select(starts_with("delta[") & 
                                             !ends_with(",1]")))
  
  # Effect size of all unique pairwise comparisons for each order
  EM <- t(get_results %>% dplyr::select(starts_with("EM[")))
  
  # Effect size of comparisons with the reference for each order
  EM_ref <- t(get_results %>% dplyr::select(starts_with("EM_ref[")))
  
  # Predictive effects of all unique pairwise comparisons for each order
  EM_pred <- t(get_results %>% dplyr::select(starts_with("EM_pred[")))
  
  # Predictive effects of all comparisons with the reference for each order
  pred_ref <- t(get_results %>% dplyr::select(starts_with("pred_ref[")))
  
  # Effect size of all unique pairwise comparisons for a range of time-points
  log_HR <- t(get_results %>% dplyr::select(starts_with("log_HR[")))

  # Estimated missingness parameter
  phi <- if (min(item$I) == 1) {
    t(get_results %>% dplyr::select(starts_with("phi") |
                                    starts_with("mean.phi") |
                                    starts_with("mean.phi[") |
                                    starts_with("phi[")))
  } else if (min(item$I) == 0 & max(item$I) == 1) {
    t(get_results %>% dplyr::select(starts_with("phi") |
                                      starts_with("mean.phi") |
                                      starts_with("mean.phi[") |
                                      starts_with("phi["))) * 
      matrix(rep(item$I , 9), nrow = item$dpoints, ncol = 9)
  } else if (max(item$I) == 0) {
    NULL
  }
  
  # Between-trial standard deviation
  tau <- t(get_results %>% dplyr::select(starts_with("tau")))
  
  # Total residual deviance
  dev <- t(get_results %>% dplyr::select(starts_with("totresdev")))
  
  # Trial-arm deviance contribution for observed outcome
  devi_o <- t(get_results %>% dplyr::select(starts_with("dev_o")))
  
  # Fitted/predicted outcome
  rhat_o <- t(get_results %>% dplyr::select(starts_with("rhat_o[")))
  
  # Ranking probability of each intervention for every rank (PH model)
  effectiveness1 <- t(get_results %>% dplyr::select(
    starts_with("effectiveness1")))
  
  # Ranking probability of each intervention for every rank (first order)
  effectiveness2 <- t(get_results %>% dplyr::select(
    starts_with("effectiveness2")))
  
  # Ranking probability of each intervention for every rank (second order)
  effectiveness3 <- t(get_results %>% dplyr::select(
    starts_with("effectiveness3")))
  
  # SUrface under the Cumulative RAnking curve values (PH model)
  SUCRA1 <- t(get_results %>% dplyr::select(starts_with("SUCRA1[")))
  
  # SUrface under the Cumulative RAnking curve values (first order)
  SUCRA2 <- t(get_results %>% dplyr::select(starts_with("SUCRA2[")))
  
  # SUrface under the Cumulative RAnking curve values (second order)
  SUCRA3 <- t(get_results %>% dplyr::select(starts_with("SUCRA3[")))
  
  # Cumulative survival probability per intervention and across the timepoints
  Surv <- t(get_results %>% dplyr::select(starts_with("Surv[")))
  
  # Expected survival time per intervention and across the timepoints
  expect_Surv <- t(get_results %>% dplyr::select(starts_with("expect_Surv[")))
  
  ## Calculate the deviance at posterior mean of fitted values
  # Correction for zero events in trial-arm. Remove when z equals 0
  r0 <- subset(ifelse(item$r == 0, item$r + 0.01, 
                      ifelse(item$r == item$c, item$r - 0.01, item$r)), 
               item$c > 0)
  
  # Remove when z equals 0
  c0 <- subset(item$c, item$c > 0)
  rhat_o_new <- subset(rhat_o[, 1], item$c > 0)
  dev_o <- subset(devi_o[, 1], item$c > 0)
  
  # Deviance at the posterior mean of the fitted response
  dev_post_o <- 2*(r0*(log(r0) - log(as.vector(rhat_o_new))) + 
                     (c0 - r0)*(log(c0 - r0) - log(c0 - as.vector(rhat_o_new))))
  
  # Sign of the difference between observed and fitted response
  sign_o <- sign(r0 - as.vector(rhat_o_new))
  
  # Obtain the leverage for observed outcomes
  lev_o <- as.vector(dev_o) - dev_post_o
  
  # Number of effective parameters
  pD <- dev[, 1] - sum(dev_post_o)
  
  # Deviance information criterion
  DIC <- pD + dev[, 1]
  
  # Measures of model assessment: DIC, pD, and total residual deviance
  model_assessment <- data.frame(DIC, pD, dev[, 1])
  colnames(model_assessment) <- c("DIC", "pD", "deviance")
  
  # Return a list of results
  results <- list(Beta = Beta, 
                  EM = EM, 
                  EM_ref = EM_ref,
                  log_HR = log_HR,
                  phi = phi,
                  SUCRA1 = SUCRA1, 
                  SUCRA2 = SUCRA2, 
                  SUCRA3 = SUCRA3, 
                  effectiveness1 = effectiveness1, 
                  effectiveness2 = effectiveness3, 
                  effectiveness3 = effectiveness3, 
                  Surv = Surv, 
                  expect_Surv = expect_Surv,
                  model_assessment = model_assessment,
                  dev_o = dev_o,
                  sign_o = sign_o,
                  lev_o = lev_o,
                  D = D,
                  max_time = max_time,
                  model = model,
                  assumption = assumption,
                  jagsfit = jagsfit)
  
  if (model == "FE") {
    results
  } else {
    results <- append(results, list(delta = delta,
                                    EM_pred = EM_pred,
                                    pred_ref = pred_ref,
                                    tau = tau))
  }

  return(results)
}
