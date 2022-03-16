#' Determine the mean value of the distribution of the missingness parameter
#'
#' @description
#'   Generates the mean value of the distribution of the informative missingness 
#'   hazard ratio parameter (submitted) in the proper format depending on the 
#'   assumed structure of the missingness parameter. \code{\link{model_fppm}} 
#'   inherits \code{missingness_param_prior} through the argument 
#'   \code{mean_misspar}.
#'
#' @param assumption Character string indicating the structure of the
#'   informative missingness hazard ratio parameter.
#'   Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"},
#'   \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"},
#'   \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, or \code{"IND-UNCORR"}.
#'   The default argument is \code{"IDE-ARM"}. The abbreviations \code{"IDE"},
#'   \code{"HIE"}, and \code{"IND"} stand for identical, hierarchical and
#'   independent, respectively. \code{"UNCORR"} stands for uncorrelated.
#' @param mean_misspar A numeric value or a vector of two numeric values for the
#'   mean of the normal distribution of the informative missingness parameter
#'   in the logarithmic scale (see 'Details'). The default argument is 0 and 
#'   corresponds to the missing-at-random assumption.
#'
#' @return A value or vector to be passed to \code{\link{model_fppm}}.
#'
#' @details \code{\link{model_fppm}} considers the informative missingness 
#'   hazard ratio in the logarithmic scale (submitted).
#'
#'   When \code{assumption} is trial-specific (i.e., \code{"IDE-TRIAL"} or
#'   \code{"HIE-TRIAL"}), or independent (i.e., \code{"IND-UNCORR"}), only one 
#'   numeric value can be assigned to \code{mean_misspar} because the same 
#'   missingness scenario is applied to all trials and trial-arms of the 
#'   dataset, respectively. When \code{assumption} is \code{"IDE-ARM"} or 
#'   \code{"HIE-ARM"}, a maximum of two \emph{different} numeric values can be 
#'   assigned as a vector to \code{mean_misspars}: the first value refers to the
#'   experimental arm, and the second value refers to the control arm of a 
#'   trial. The first value is considered for all non-reference interventions 
#'   and the second value is considered for the reference intervention of the 
#'   network (i.e., the intervention with identifier equal to one). 
#'   This is necessary to ensure transitivity in the assumptions for the 
#'   missingness parameter across the comparisons in the network 
#'   (Spineli, 2019b). When one numeric value is considered, the same 
#'   missingness scenario is applied to all interventions in the dataset.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{model_fppm}}
#'
#' @references
#' Kalyvas C, Papadimitropoulou K, Malbecq W, Spineli LM. Dealing with censoring
#' in a network meta-analysis of time-to-event data. 2022 (submitted)
#'
#' Spineli LM. Modeling missing binary outcome data while preserving
#' transitivity assumption yielded more credible network meta-analysis
#' results. \emph{J Clin Epidemiol} 2019b;\bold{105}:19--26.
#'
missingness_param_prior <- function(assumption, mean_misspar) {

  # Specifying the 'mean_misspar' for the missingness parameter
  if (!is.element(assumption, c("IDE-ARM",
                                "IDE-TRIAL",
                                "IDE-COMMON",
                                "HIE-ARM",
                                "HIE-TRIAL",
                                "HIE-COMMON",
                                "IND-UNCORR"))) {
    stop("Insert 'IDE-ARM', 'IDE-TRIAL', 'IDE-COMMON', 'HIE-ARM', 
         'HIE-TRIAL', 'HIE-COMMON', or 'IND-UNCORR' ", call. = FALSE)
  } else if (is.element(assumption, c("HIE-ARM", "IDE-ARM")) &
             missing(mean_misspar)) {
    mean_misspar <- rep(0.0001, 2)
  } else if (is.element(assumption, c("IDE-TRIAL",
                                      "IDE-COMMON",
                                      "HIE-TRIAL",
                                      "HIE-COMMON",
                                      "IND-UNCORR")) &
             missing(mean_misspar)) {
    mean_misspar <- 0.0001
  } else if (is.element(assumption, c("HIE-ARM", "IDE-ARM")) &
             (length(mean_misspar) != 2)) {
    stop("'mean_misspar' must be a vector of two numeric values", call. = FALSE)
  } else if (is.element(assumption, c("HIE-ARM", "IDE-ARM")) &
             (length(mean_misspar) == 2)) {
    mean_misspar <- as.vector(mean_misspar)
    mean_misspar[1] <- ifelse(mean_misspar[1] == 0, 0.0001, mean_misspar[1])
    mean_misspar[2] <- ifelse(mean_misspar[2] == 0, 0.0001, mean_misspar[2])
  } else if (is.element(assumption, c("IDE-TRIAL",
                                      "IDE-COMMON",
                                      "HIE-TRIAL",
                                      "HIE-COMMON",
                                      "IND-UNCORR")) &
             (length(mean_misspar) > 1)) {
    stop("'mean_misspar' must be a scalar", call. = FALSE)
  } else if (is.element(assumption, c("IDE-TRIAL",
                                      "IDE-COMMON",
                                      "HIE-TRIAL",
                                      "HIE-COMMON",
                                      "IND-UNCORR")) &
             (length(mean_misspar) == 1)) {
    mean_misspar <- ifelse(mean_misspar == 0, 0.0001, mean_misspar)
  }

  return(mean_misspar)
}
