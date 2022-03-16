#' Leverage plot 
#'
#' @description Plots the leverage against the square root of the
#'   posterior mean of residual deviance of the trial-arms for each time-point 
#'   under the network meta-analysis model with fractional polynomials.
#'
#' @param full An object of S3 class \code{\link{model_fppm}}.
#'   See 'Value' in \code{\link{model_fppm}}.
#' @param title A title to indicate the model.
#'
#' @return A scatterplot of the leverage against the square root of the
#'   posterior mean of residual deviance of the trial-arms for each time-point. 
#'   The green, yellow, and red curves correspond to the parabola
#'   \eqn{x^2 + y = k} with \eqn{k} = 1, 2, and 3, respectively. The data-points
#'   correspond to trial-arms and time-points. The number of points equals the number 
#'   of rows of \code{data_points}. Data points found outside the yellow 
#'   parabola are linked with a number. The number refers to the position of the
#'   trial-arm and time-point in the dataset (see 'Arguments' and 'Value' in 
#'   \code{\link{data_preparation}}).  These data-points contribute more than 1 
#'   to the deviance information criterion and, hence, the model's poor fit.
#'
#' @details \code{leverage_plot} is integrated in the \code{\link{results_plot}}
#'   function to return the leverage plot alongside other necessary plots. 
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{data_preparation}}, \code{\link{model_fppm}},
#'   \code{\link{results_plot}}
#'
leverage_plot <- function(full, title) {

  # Posterior mean of deviance contribution (observed outcomes)
  dev_o <- full$dev_o

  # Leverage (observed outcomes)
  lev_o <- full$lev_o

  # The sign of the difference between observed and fitted outcome
  sign_o <- full$sign_o

  # A function to extract numbers from a character.
  # Source: http://stla.github.io/stlapblog/posts/numextract.html
  numextract <- function(string) {
    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  }

  # Create a vector of trial-arms using the 'paste0' function
  trial_arm <- as.numeric(numextract(names(dev_o)))


  # Prepare the data-frame for the leverage plot (using ggplot2)
  prepare_lev <- round(data.frame(sign_o * sqrt(as.vector(dev_o)), lev_o),
                       2)
  colnames(prepare_lev) <- c("signed_dev_o", "lev_o")

  # Keep only trial-arms that exceed the parabola y + (x)^2 = c at c = 2.
  poor_o <- ifelse(prepare_lev$lev_o > 2 - (prepare_lev$signed_dev_o^2) |
                     prepare_lev$lev_o < - (2 - (prepare_lev$signed_dev_o^2)),
                   trial_arm, NA)
  poor_fit_o <- data.frame(prepare_lev[!is.na(poor_o), 1:2],
                           poor_o[!is.na(poor_o)])
  colnames(poor_fit_o) <- c("signed_dev", "leverage", "poor")

  # Leverage plot for observed outcomes
  observed <- ggplot(data = prepare_lev, aes(x = signed_dev_o, y = lev_o)) +
               geom_point(size = 2, colour = "black") +
               geom_smooth(aes(x = signed_dev_o,
                               y = 1 - (signed_dev_o^2)),
                           method = "loess",
                           formula = "y ~ x",
                           colour = "#009E73",
                           linetype = 2) +
               geom_smooth(aes(x = signed_dev_o,
                               y = 2 - (signed_dev_o^2)),
                           method = "loess",
                           formula = "y ~ x",
                           colour = "orange",
                           linetype = 2) +
               geom_smooth(aes(x = signed_dev_o,
                               y = 3 - (signed_dev_o^2)),
                           method = "loess",
                           formula = "y ~ x",
                           colour = "#D55E00",
                           linetype = 2) +
               geom_text_repel(data = poor_fit_o,
                               aes(x = signed_dev,
                                   y = leverage,
                                   label = poor),
                         color = "blue",
                         fontface = "bold",
                         hjust = "right",
                         size = 3.8,
                         max.overlaps = Inf,
                         nudge_x = -0.1,
                         direction = "y") +
               labs(x = expression(
                 "" %+-% sqrt("Posterior mean of the residual deviance")),
                    y = "Leverage (each data point)") +
               coord_cartesian(xlim = c(min(prepare_lev$signed_dev_o),
                                        max(prepare_lev$signed_dev_o)),
                               ylim = c(0,
                                        max(3 - (prepare_lev$signed_dev_o^2))),
                               expand = TRUE) +
               ggtitle(title) +
               theme_classic() +
               theme(axis.title.x = element_text(color = "black", size = 12),
                     axis.title.y = element_text(color = "black", size = 12),
                     axis.text.x = element_text(color = "black", size = 12),
                     axis.text.y = element_text(color = "black", size = 12),
                     plot.title = element_text(color = "black", size = 11,
                                               face = "bold"))

  return(observed)
}
