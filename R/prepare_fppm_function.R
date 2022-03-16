#' WinBUGS code for Bayesian network meta-analysis with fractional polynomials
#'
#' @description
#'   The WinBUGS code, as written by Jansen (2011) to run a one-stage
#'   Bayesian network meta-analysis with fractional polynomials, extended to 
#'   incorporate the pattern-mixture model for aggregate (binary) missing 
#'   participant outcome data. 
#'
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
#'   independent, respectively. \code{"UNCORR"} stands for uncorrelated.
#'
#' @return An R character vector object to be passed to \code{\link{model_fppm}}
#'   through the \code{\link[base:textConnection]{textConnection}} function as 
#'   the argument \code{object}.
#'
#' @details \code{prepare_fppm} creates the model in the JAGS dialect of the 
#'   BUGS language. The output of this function constitutes the argument
#'   \code{model.file} of the \code{\link[R2jags:jags]{jags}} function (in the
#'   R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}) via the
#'   \code{\link[base:textConnection]{textConnection}} function.
#'
#' @author {Chrysostomos Kalyvas}, {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_fppm}}, \code{\link[R2jags:jags]{jags}}, 
#'   \code{\link[base:textConnection]{textConnection}}
#'
#' @references
#' Jansen JP. Network meta-analysis of survival data with fractional 
#' polynomials. \emph{BMC Med Res Methodol} 2011;\bold{11}:61.
#' \doi{10.1186/1471-2288-11-61}
#'
#' Kalyvas C, Papadimitropoulou K, Malbecq W, Spineli LM. Dealing with censoring
#' in a network meta-analysis of time-to-event data. 2022 (submitted)
#'
prepare_fppm <- function(model, assumption) {
  
  stringcode <- "model {
                    for (i in 1:dpoints) {
                      timen1[i] <- equals(P1, 0)*log(time[i]) + (1 - equals(P1, 0))*pow(time[i], P1)
                      timen2[i] <- (1 - equals(P2, P1))*(equals(P2, 0)*log(time[i]) + 
                      (1 - equals(P2, 0))*pow(time[i], P2)) + equals(P2, P1)*equals(P2, 0)*log(time[i])*log(time[i]) + 
                      (1 - equals(P2, 0))*pow(time[i], P2)*log(time[i])
                      r[i] ~ dbin(p_o[i], c[i])
                      q[i] <- q0[i]*I[i]
                      m[i] ~ dbin(q0[i], n[i]) 
                      q0[i] ~ dunif(0, 1)
                      p_o[i] <- 1 - exp(-h_o[i]*dt[i])
                      h_o[i] <- h[i]/(1 - q[i]*(1 - theta[i]))
                      log(h[i]) <- Beta[s[i], a[i], 1] + Beta[s[i], a[i], 2]*timen1[i] + Beta[s[i], a[i], 3]*timen2[i]
                      rhat_o[i] <- p_o[i]*c[i]
                      dev_o[i] <- 2*(r[i]*(log(r[i]) - log(rhat_o[i])) + (c[i] - r[i])*(log(c[i] - r[i]) - log(c[i] - rhat_o[i])))
                    }
                    totresdev_o <- sum(dev_o[])
                    for (i in 1:ns) {
                      for (k in 1:3) {
                         mu[i, k] ~ dnorm(0, 0.0001)
                      }
                    }\n"
                      
  stringcode <- if (model == "FE") {
    paste(stringcode, "for (i in 1:ns) {
                         for (j in 1:na[i]) {
                           Beta[i, j, 1] <- mu[i, 1] + d[t[i, j], 1] - d[t[i, 1], 1]
                           Beta[i, j, 2] <- mu[i, 2] + d[t[i, j], 2] - d[t[i, 1], 2]
                           Beta[i, j, 3] <- mu[i, 3] + d[t[i, j], 3] - d[t[i, 1], 3]
                         }
                       }
                       for (k in 1:3) {
                         for (t in 1:(ref - 1)) {
                           EM_ref[t, k] <- d[t, k]
                         }
                         for (t in (ref + 1):nt) {
                           EM_ref[t, k] <- d[t, k]
                         }
                       }
                       for (k in 1:3) {
                         for (t in 1:(nt - 1)) {
                           for (c in (t + 1):nt) {
                             EM[c, t, k] <- d[c, k] - d[t, k]
                           }
                         }
                       }\n")
  } else {
    paste(stringcode, "for (i in 1:ns) {
                         w[i, 1] <- 0
                         delta[i, 1] <- 0
                         for (j in 1:na[i]) {
                           Beta[i, j, 1] <- mu[i, 1] + delta[i, j]
                           Beta[i, j, 2] <- mu[i, 2] + d[t[i, j], 2] - d[t[i, 1], 2]
                           Beta[i, j, 3] <- mu[i, 3] + d[t[i, j], 3] - d[t[i, 1], 3]
                         }
                         for (j in 2:na[i]) {
                           delta[i, j] ~ dnorm(md[i, j], precd[i, j])
                           md[i, j] <- d[t[i, j], 1] - d[t[i, 1], 1] + sw[i, j]
                           precd[i, j] <- 2*prec*(j - 1)/j
                           w[i, j] <- delta[i, j] - (d[t[i, j], 1] - d[t[i, 1], 1])
                           sw[i, j] <- sum(w[i, 1:(j - 1)])/(j - 1)
                         }
                       } 
                       prec <- pow(tau, -2)
                       tau_a ~ dnorm(0, heter_prior[2])I(0, )
                       tau_b ~ dunif(0, heter_prior[2])
                       tau <- tau_a*equals(heter_prior[3], 1) + tau_b*equals(heter_prior[3], 2)
                       for (k in 1:3) {
                         for (t in 1:(ref - 1)) {
                           EM_ref[t, k] <- d[t, k]
                           pred_ref[t, k]~ dnorm(EM_ref[t, k], prec)
                         }
                         for (t in (ref + 1):nt) {
                           EM_ref[t, k] <- d[t, k]
                           pred_ref[t, k]~ dnorm(EM_ref[t, k], prec)
                         }
                       }
                       for (k in 1:3) {
                         for (t in 1:(nt - 1)) {
                           for (c in (t + 1):nt) {
                             EM[c, t, k] <- d[c, k] - d[t, k]
                             EM_pred[c, t, k] ~ dnorm(EM[c, t, k], prec)
                           }
                         }
                       }\n")
  }

  stringcode <- paste(stringcode, "for (k in 1:3) {
                                     d[ref, k] <- 0
                                     for (t in 1:(ref - 1)) {
                                       d[t, k] ~ dnorm(0, 0.0001)
                                     }
                                     for (t in (ref + 1):nt) {
                                       d[t, k] ~ dnorm(0, 0.0001)
                                     }
                                   }
                                   sorted1 <- rank(d[, 1])
                                   sorted2 <- rank(d[, 2])
                                   sorted3 <- rank(d[, 3])
                                   for (t in 1:nt) {
                                     order1[t] <- (nt + 1 - sorted1[t])*equals(D, 1) + sorted1[t]*(1 - equals(D, 1))
                                     order2[t] <- (nt + 1 - sorted2[t])*equals(D, 1) + sorted2[t]*(1 - equals(D, 1))
                                     order3[t] <- (nt + 1 - sorted3[t])*equals(D, 1) + sorted3[t]*(1 - equals(D, 1))
                                     for (l in 1:nt) {
                                       effectiveness1[t, l] <- equals(order1[t], l)
                                       effectiveness2[t, l] <- equals(order2[t], l)
                                       effectiveness3[t, l] <- equals(order3[t], l)
                                       cumeffectiveness1[t, l] <- sum(effectiveness1[t, 1:l])
                                       cumeffectiveness2[t, l] <- sum(effectiveness2[t, 1:l])
                                       cumeffectiveness3[t, l] <- sum(effectiveness3[t, 1:l])
                                     }
                                     SUCRA1[t] <- sum(cumeffectiveness1[t, 1:(nt - 1)])/(nt - 1)
                                     SUCRA2[t] <- sum(cumeffectiveness2[t, 1:(nt - 1)])/(nt - 1)
                                     SUCRA3[t] <- sum(cumeffectiveness3[t, 1:(nt - 1)])/(nt - 1)
                                   }
                                   for (m in 1:maxt) {
                                     time1[m] <- equals(P1, 0)*log(m) + (1 - equals(P1, 0))*pow(m, P1)
                                     time2[m] <- (1 - equals(P2, P1))*(equals(P2, 0)*log(m) + (1 - equals(P2, 0))*pow(m, P2)) + 
                                     equals(P2, P1)*equals(P2, 0)*log(m)*log(m) + (1 - equals(P2, 0))*pow(m, P2)*log(m)
                                   }
                                   for (k in 1:(nt - 1)) {
                                     for (c in (k + 1):nt) {
                                       for (m in 1:maxt) {
                                         log_HR[c, k, m] <- (d[c, 1] - d[k, 1]) + (d[c, 2] - d[k, 2])*time1[m] + (d[c, 3] - d[k, 3])*time2[m]
                                       }
                                     }
                                   }
                                   for (k in 1:3) {
                                     mu_mean[k] <- mean(mu[1:ns, k])
                                     for (n in 1:nt) {
                                       beta[n, k] <- mu_mean[k] + d[n, k]
                                     }
                                   }
                                   for (n in 1:nt) {
                                     for (m in 1:maxt) {
                                       log(hazard[n, m]) <- beta[n, 1] + (beta[n, 2]*time1[m]) + (beta[n, 2]*time2[m])
                                       cum_h[n, m] <- sum(hazard[n, 1:m])
                                       T[n, m] <- 1 - exp(-cum_h[n, m])
                                       Surv[n, m] <- 1 - T[n, m]
                                     }
                                     expect_Surv[n] <- sum(Surv[n, 1:maxt])  # Median survival time
                                   }\n")
  
  stringcode <- if (assumption == "IDE-COMMON") {
    paste(stringcode, "for (i in 1:dpoints) {
                         log(theta[i]) <- phi
                       }
                       phi ~ dnorm(meand_phi, precd_phi)\n")
  } else if (assumption == "IDE-TRIAL") {
    paste(stringcode, "for (i in 1:dpoints) {
                         log(theta[i]) <- phi[s[i]]
                       }
                       for (i in 1:ns) {
                         phi[i] ~ dnorm(meand_phi, precd_phi)
                       }\n")
  } else if (assumption == "IDE-ARM") {
    paste(stringcode, "for (i in 1:dpoints) {
                         log(theta[i]) <- phi[t[s[i], a[i]]]
                       }
                       phi[ref] ~ dnorm(meand_phi[2], precd_phi)
                       for (t in 1:(ref - 1)) {
                         phi[t] ~ dnorm(meand_phi[1], precd_phi)
                       }
                       for (t in (ref + 1):nt) {
                         phi[t] ~ dnorm(meand_phi[1], precd_phi)
                       }\n")
  } else if (assumption == "HIE-COMMON") {
    paste(stringcode, "for (i in 1:dpoints) {
                         log(theta[i]) <- phi[i]
                         phi[i] ~ dnorm(mean_phi, prec_phi)
                       }
                       mean_phi ~ dnorm(meand_phi, precd_phi)
                       prec_phi <- pow(sd_phi, -2)
                       sd_phi ~ dunif(0, psi_phi)
                       psi_phi <- pow(precd_phi, -2)\n")
  } else if (assumption == "HIE-TRIAL") {
    paste(stringcode, "for (i in 1:dpoints) {
                         log(theta[i]) <- phi[i]       
                         phi[i] ~ dnorm(mean_phi[s[i]], prec_phi[s[i]])
                       }
                       for (i in 1:ns) {
                         mean_phi[i] ~ dnorm(meand_phi, precd_phi)
                         prec_phi[i] <- pow(sd_phi[i], -2)
                         sd_phi[i] ~ dunif(0, psi_phi)
                       }
                       psi_phi <- pow(precd_phi, -2)\n")
  } else if (assumption == "HIE-ARM") {
    paste(stringcode, "for (i in 1:dpoints) {
                         log(theta[i]) <- phi[i]
                         phi[i] ~ dnorm(mean_phi[t[s[i], a[i]]], prec_phi[t[s[i], a[i]]])
                       }
                       mean_phi[ref] ~ dnorm(meand_phi[2], precd_phi)
                       prec_phi[ref] <- pow(sd_phi[ref], -2)
                       sd_phi[ref] ~ dunif(0, psi_phi)
                       for (t in 1:(ref - 1)) {
                         mean_phi[t] ~ dnorm(meand_phi[1], precd_phi)
                         prec_phi[t] <- pow(sd_phi[t], -2)
                         sd_phi[t] ~ dunif(0, psi_phi)
                       }
                       for (t in (ref + 1):nt) {
                         mean_phi[t] ~ dnorm(meand_phi[1], precd_phi)
                         prec_phi[t] <- pow(sd_phi[t], -2)
                         sd_phi[t] ~ dunif(0, psi_phi)
                       }
                       psi_phi <- pow(precd_phi, -2)\n")
  } else if (assumption == "IND-UNCORR") {
    paste(stringcode, "for (i in 1:dpoints) {log(theta[i]) <- phi[i]
                         phi[i] ~ dnorm(meand_phi, precd_phi)
                       }\n")
  }
  
  stringcode <- paste(stringcode, "\n}")
  
  return(stringcode)
}
