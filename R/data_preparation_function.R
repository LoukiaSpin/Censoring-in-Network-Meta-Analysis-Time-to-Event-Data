#' Prepare the dataset in the proper format for R2jags
#'
#' @description
#'   \code{data_preparation} prepares the dataset in the proper format for
#'   R2jags and returns a list of elements that \code{\link{model_fppm}} 
#'   inherits via the arguments \code{data_points} and \code{data_trial}.
#'
#' @param data_points A data-frame in long format with trial-arm data for each 
#'   time-point. See 'Format' below.
#' @param data_trial A data-frame of the one-trial-per-row format with arm-level
#'   data for the interventions studied in each arm of every trial. 
#'   In essence, an intervention identifier in each arm.
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
#' @return A list of vectors on the following elements to be passed
#'   to \code{\link{model_fppm}}:
#'   \item{I}{The pseudo-data-frame \code{I} (see 'Details').}
#'   \item{na}{The number of intervention arms in each trial.}
#'   \item{nt}{The number of interventions in the network.}
#'   \item{ns}{The number of trials in the network.}
#'   \item{dpoints}{The number of data-points that equal the sum of 
#'    time-points per arm across the trials.}
#'
#' @details \code{data_preparation} prepares the data for the Bayesian analysis.
#'   It checks whether the element \strong{m} exists in the \code{data_points} 
#'   argument (See 'Format' in \code{\link{model_fppm}}). If this element is 
#'   missing, \code{data_preparation} creates a pseudo-data-frame for \strong{m}
#'   that has the zero value for each time-point and trial-arm, and \code{NA} 
#'   for the unobserved trial-arms for the corresponding time-point. It all
#'   creates the pseudo-data-frame \code{I} that is identical with the 
#'   pseudo-data-frame for \code{m}. If the element \strong{m} exists in the 
#'   \code{data_points} and has values only for some time-points and trial-arms,
#'   the pseudo-data-frame for \strong{m} is identical to \strong{m} for the 
#'   corresponding time-points and trial-arms, and the pseudo-data-frame 
#'   \code{I} has the value one for these corresponding data-points. 
#'   Both pseudo-data-frames aim to retain the trials without information on 
#'   missing participant outcome data.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \href{https://CRAN.R-project.org/package=R2jags}{R2jags},
#'   \code{\link{model_fppm}}
#'
data_preparation <- function(data_points, data_trial) {

  # Intervention studied in each arm of every trial
  time <- if (is.null(data_points$time)) {
    stop("data_points: 'time' is missing", call. = FALSE)
  } else {
    data_points$time
  }
  r <- if (is.null(data_points$r)) {
    stop("data_points: the number of events, 'r', is missing", 
         call. = FALSE)
  } else {
    data_points$r
  }
  m <- if (is.null(data_points$m)) {
    message("Missing participant outcome data have *not* been collected")
    rep(NA, dim(data_points)[1])
  } else {
    data_points$m
  }
  c <- if (is.null(data_points$c)) {
    stop("data_points: the number of completers, 'c', is missing", 
         call. = FALSE)
  } else {
    data_points$c
  }
  n <- if (is.null(data_points$n)) {
    stop("data_points: the number randomised, 'n', is missing", 
         call. = FALSE)
  } else {
    data_points$n
  }
  a <- if (is.null(data_points$a)) {
    stop("data_points: the intervention arm, 'a', is missing", 
         call. = FALSE)
  } else {
    data_points$a
  }
  s <- if (is.null(data_points$s)) {
    stop("data_points: the study identifier, 's', is missing", 
         call. = FALSE)
  } else {
    data_points$s
  }
  dt <- if (is.null(data_points$dt)) {
    stop("data_points: the length of time-interval, 'dt', is missing", 
         call. = FALSE)
  } else {
    data_points$dt
  }
  treat <- if (dim(data_trial %>% select(starts_with("t")))[2] == 0) {
    stop("data_trial: the matrix on interventions, 't', is missing", 
         call. = FALSE)
  } else {
    data_trial %>% select(starts_with("t"))
  }
  # Total number of included trials
  ns <- length(data_trial[, 1])
  # Number of interventions investigated in every trial
  na <- apply(treat, 1, function(x) length(which(!is.na(x))))
  # Total number of interventions
  nt <- length(table(as.matrix(treat)))
  # The intervention with identifier '1' is the reference of the network
  dpoints <-  dim(data_points)[1]

  # Indicate the trials with fully reported missing outcome data (for each arm)
  # If a trial reports missing participant outcome data partially or not at all,
  # insert 'NA'.
  I <- m_new <- rep(NA, dpoints)
  for (i in 1:dpoints) { 
    I[i] <- if (is.na(m[i])) {
      0
    } else if (!is.na(m[i])) {
      1
    }
 
    m_new[i] <- if (is.na(m[i])) {
      0
    } else if (!is.na(m[i])) {
      m[i]
    }
  }

  # List of arguments
  results <- list(time = time,
                  r = r,
                  mod = m_new,
                  I = I,
                  c = c,
                  n = n,
                  a = a,
                  s = s,
                  dt = dt,
                  t = treat,
                  na = na,
                  nt = nt,
                  ns = ns,
                  dpoints = dpoints)

  return(results)
}
