#' Determine the prior distribution for the heterogeneity parameter
#'
#' @description
#'   Generates the prior distribution (weakly informative or empirically-based)
#'   for the heterogeneity parameter.
#'   \code{\link{model_fppm}} inherits \code{heterogeneity_param_prior} via the
#'   argument \code{heter_prior}.
#'
#' @param model Character string indicating the analysis model with values
#'   \code{"RE"}, or \code{"FE"} for the random-effects and fixed-effect model,
#'   respectively. The default argument is \code{"RE"}.
#' @param heter_prior A list of three elements with the following order:
#'   1) a character string indicating the distribution with
#'   (currently available) values \code{"halfnormal"}, or \code{"uniform"}; 
#'   2) two numeric values that refer to the parameters of the selected 
#'   distribution. For \code{"halfnormal"}, these numbers refer to zero and the 
#'   scale parameter (equal to 4 or 1 being the corresponding precision of the 
#'   scale parameter 0.5 or 1). For \code{"uniform"}, these numbers refer to the
#'   minimum and maximum value of the distribution.
#'
#' @return A value to be passed to \code{\link{mode_fppm}}.
#'
#' @details
#'   The names of the (current) prior distributions follow the JAGS syntax.
#'   The users may refer to Dias et al., (2013) to determine the minimum and
#'   maximum value of the uniform distribution, and to Friede et al., (2017)
#'   to determine the mean and precision of the half-normal distribution.
#'   When \code{model} is \code{"FE"}, \code{heterogeneity_param_prior}
#'   is ignored in \code{\link{model_fppm}}.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{model_fppm}}
#'
#' @references
#' Friede T, Roever C, Wandel S, Neuenschwander B. Meta-analysis of two studies
#' in the presence of heterogeneity with applications in rare diseases.
#' \emph{Biom J} 2017;\bold{59}(4):658--671.
#' \doi{10.1002/bimj.201500236}
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
#' making 2: a generalized linear modeling framework for pairwise and
#' network meta-analysis of randomized controlled trials.
#' \emph{Med Decis Making} 2013;\bold{33}(5):607--617.
#' \doi{10.1177/0272989X12458724}
#'
heterogeneity_param_prior <- function(model, heter_prior) {

  # Specifying the prior distribution for the between-trial parameter
  if (model == "RE" &
      missing(heter_prior)) {
    stop("The argument 'heter_prior' needs to be defined", call. = FALSE)
  } else if (model == "FE" &
             missing(heter_prior)) {
    list(NA, NA, NA)
  } else if (model == "FE") {
    list(NA, NA, NA)
  } else if (model == "RE" &
             heter_prior[[1]] != "halfnormal" &
             heter_prior[[1]] != "uniform") {
    stop("Insert 'halfnormal', or 'uniform'", call. = FALSE)
  } else if (model == "RE" &
             heter_prior[[1]] == "halfnormal") {
    as.numeric(c(0, heter_prior[[3]], 1))
  } else if (model == "RE" &
             heter_prior[[1]] == "uniform") {
    as.numeric(c(0, heter_prior[[3]], 2))
  } 
}
