#' End-user-ready results for network meta-analysis with fractional polynomials 
#'
#' @description \code{fppm_plot} hosts a toolkit of functions that facilitates
#'   the visualisation and tabulation of results from the Bayesian network 
#'   meta-analysis model with fractional polynomials via the 
#'   \code{\link{model_fppm}} function.
#'
#' @param full An object of S3 class \code{\link{model_fppm}}. See 'Value' in
#'   \code{\link{model_fppm}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data_trial} of
#'   \code{\link{model_fppm}}.
#' @param control A character indicating the intervention of interest
#'   \emph{exactly} as defined in \code{drug_names}.
#' @param time_title A x axis label to indicate how the scale of time-points 
#'   (e.g., years, months).
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   to an 'xlsx' file (via the \code{\link[writexl:write_xlsx]{write_xlsx}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) to the working
#'   directory of the user. The default is \code{FALSE} (do not export).
#'
#' @return \code{fppm_plot} returns the following list of elements:
#'   \item{table_effect_size}{A data-frame with the posterior median, posterior 
#'   standard deviation, and 95\% credible interval of the summary hazard ratio 
#'   (in the logarithmic scale) for each possible pairwise comparisons in the 
#'   network and each order of the fractional polynomial model, namely, 
#'   proportional hazards model, first- and second-order (Jansen, 2011; 
#'   Kalyvas et al., 2022)).}
#'   \item{table_model_assessment}{The DIC, number of effective
#'   parameters, and total residual deviance (Spiegelhalter et al., 2002).}
#'   \item{table_exp_survival}{A data-frame with the posterior median, and 95\%
#'   credible interval of the expected survival time.}
#'   \item{'Between-trial variance'}{The posterior median and 95\% credible 
#'   interval of \emph{tau}. When a fixed-effect model has been performed, 
#'   \code{model_fppm} does not return this element.}
#'   \item{'Leverage plot'}{The leverage plot. See 'Details' and 'Value' in 
#'   \code{\link{leverage_plot}}.}
#'   \item{'Hazard ratio'}{A density plot on the hazard ratio of each possible 
#'   pairwise comparison for a range of time-points. The maximum time-point is 
#'   specified in the argument \code{max_time} of the \code{\link{model_fppm}} 
#'   function. See 'Details' and 'Value' in \code{\link{plot_hr_surv}}.}
#'   \item{'Overall survival'}{A plot on the survival probability of each 
#'   interventions in the network for a range of time-points. 
#'   The maximum time-point is specified in the argument \code{max_time} of the 
#'   \code{\link{model_fppm}} function. See 'Details' and 'Value' in 
#'   \code{\link{plot_hr_surv}}.}
#'
#' @details The function exports \code{table_effect_size} and
#'   \code{table_model_assessment} to separate 'xlsx' files (via the
#'   \code{\link[writexl:write_xlsx]{write_xlsx}} function) to the working
#'   directory of the user.
#'
#'   \code{fppm_plot} can be used only for a network of interventions. In the
#'   case of two interventions, the execution of the function will be stopped
#'   and an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{model_fppm}}, \code{\link{plot_hr_surv}},
#'   \code{\link[writexl:write_xlsx]{write_xlsx}}
#'
#' @references
#' Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A. Bayesian measures of
#' model complexity and fit. \emph{J R Stat Soc B} 2002;\bold{64}:583--616.
#' \doi{10.1111/1467-9868.00353}
#' 
#' Jansen JP. Network meta-analysis of survival data with fractional 
#' polynomials. \emph{BMC Med Res Methodol} 2011;\bold{11}:61.
#' \doi{10.1186/1471-2288-11-61}
#'
#' Kalyvas C, Papadimitropoulou K, Malbecq W, Spineli LM. Dealing with censoring
#' in a network meta-analysis of time-to-event data. 2022 (submitted)
#'
fppm_plot <- function (full, 
                       drug_names, 
                       control, 
                       time_title, 
                       save_xls) {
  
  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
  }
  
  # Matrix of possible pairwise comparisons (control versus experimental)
  comp12 <- t(combn(drug_names, 2))
  
  # Matrix of possible pairwise comparisons (experimental versus control)
  comp21 <- cbind(comp12[, 2], comp12[, 1])
  
  # Repeat comp21 three times (one per order)
  comp <- matrix(rep(t(comp21), 3), ncol = ncol(comp21), byrow = TRUE)
  colnames(comp) <- c("exper", "ctrl")
  
  # Effect size per order
  effect_size0 <- round(full$EM[, c("50%", "sd", "2.5%", "97.5%")], 3)
  colnames(effect_size0) <- c("median", "sd", "lower", "upper")
  effect_size <- data.frame(effect_size0, comp, rep(1:3, each = dim(comp12)[1]))
  colnames(effect_size)[7] <- c("model")

  # Model assessment parameters
  model_assessment <- full$model_assessment
  
  # Expected survival time
  exp_surv <- round(full$expect_Surv[, c("50%", "2.5%", "97.5%")], 3)
  rownames(exp_surv) <- drug_names
  colnames(exp_surv) <- c("median", "lower", "upper")
  
  # Between-trial variance (proportional hazards model)
  tau <- if (full$model == "RE") {
    round(full$tau[, c("50%", "2.5%", "97.5%")], 3)
  } else {
    "Not applicable"
  }
  
  # Leverage plot
  lev_fp <- leverage_plot(full, 
                          title = "Fractional polynomial")
  
  # Plot hazard ratio 
  plot_hr <- plot_hr_surv(full, 
                          type = "HR",
                          drug_names = drug_names,
                          control = control,
                          time_title = time_title)
  
  # Plot overall survival
  plot_surv <- plot_hr_surv(full, 
                            type = "Survival",
                            drug_names = drug_names,
                            control = control,
                            time_title = time_title)
  
  # Write the table with the results as .xlsx
  if (save_xls == TRUE) {
    write_xlsx(exp_surv, paste0("Table Expected Survival Time", ".xlsx"))
    write_xlsx(model_assessment,
               paste0("Table Model Assessment", ".xlsx"))
    write_xlsx(tau,
               paste0("Table Between-trial Variance", ".xlsx"))
  }
  
  # Obtain results
  results <- list(table_effect_size = 
                    knitr::kable(effect_size, 
                                 align = "llllccc", 
                                 caption = "Log hazard ratios"),
                  table_model_assessment = 
                    knitr::kable(model_assessment, 
                                 caption = "Model assessment parameters"),
                  table_exp_survival = 
                    knitr::kable(exp_surv, 
                                 align = "lll", 
                                 caption = "Expected survival time"),
                  'Between-trial variance' = tau,
                  'Leverage plot' = lev_fp,
                  'Hazard ratio' = plot_hr,
                  'Overall survival' = plot_surv)
  
  # Subset based on the model type
  results <- if (full$model == "RE") {
    results
  } else {
    results[c(-3)] 
  } 
  
  return(results)
}
