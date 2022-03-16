## Hazard ratio and Overall Survival from 1 to max time.
plot_hr_surv <- function (full, 
                          type, 
                          drug_names, 
                          control, 
                          time_title) {
  
  drug_names <- if (missing(drug_names)) {
    stop("The argument 'drug_names' must be defined", .call = FALSE)
  } else {
    drug_names
  }
  
  type <- if (missing(type)) {
    stop("The argument 'type' must be defined", .call = FALSE)
  } else {
    type
  }
  
  control <- if (is.null(control)) {
    stop("The argument 'control' must be defined", .call = FALSE)
  } else if (!is.element(control, drug_names)) {
    stop("The argument 'control' must be a substet of 'drug_names'", 
         .call = FALSE)
  } else if (is.element(control, drug_names)) {
    control
  } 
 
  time_title <- if (missing(time_title)) {
    stop("The argument 'time_title' must be defined", .call = FALSE)
  } else {
    time_title
  }
  
  time_max <- full$max_time
  
  # NUmber of interventions
  num_treat <- length(drug_names)
  
  # Matrix of possible pairwise comparisons (control versus experimental)
  comp12 <- t(combn(drug_names, 2))
  
  # Matrix of possible pairwise comparisons (experimental versus control)
  comp21 <- cbind(comp12[, 2], comp12[, 1])
  
  # Experimental versus control arm
  temp21 <- round(full$log_HR[, c("50%", "2.5%", "97.5%")], 4)
  
  # Control versus experimental arm
  temp12 <- cbind(temp21[, 1] * (-1), temp21[, 3] * (-1), temp21[, 2] * (-1))
  
  # Repeat comp21 as many times as time_max
  comp_time21 <- matrix(rep(t(comp21), time_max), ncol = ncol(comp21), byrow = TRUE)
  
  # Repeat comp12 as many times as time_max
  comp_time12 <- matrix(rep(t(comp12), time_max), ncol = ncol(comp12), byrow = TRUE)
  
  # Merge temp21 with tem12, and similarly for the comp_times
  temp_both <- cbind(rbind(temp21, temp12), rbind(comp_time21, comp_time12))
  
  # Vector of the range 1 to time_max
  count1 <- rep(rep(1:time_max, each = dim(comp21)[1]), 2)
  
  # Prepare dataset for ggplot2 (Hazard ratio)
  df1 <- as.data.frame(cbind(temp_both, count1))
  colnames(df1) <- c("median", "lower", "upper", "exper", "ctrl", "time")
  
  # Subset defined by 'control'
  df1_sub <- subset(df1, ctrl == control)
  df1_sub$exper_new <- paste(df1_sub$exper, "vs", df1_sub$ctrl)
  
  # Prepare dataset for ggplot2 (Survival
  surv_data <- data.frame(full$Surv[, c("50%", "2.5%", "97.5%")], 
                          rep(drug_names, time_max), 
                          rep(1:time_max, each = num_treat))
  colnames(surv_data)  <- c("median", "lower", "upper", "treat", "time")
  
  # Create the plot
  if (type == "HR") {
    ggplot(df1_sub, aes(x = as.numeric(time))) + 
      stat_smooth(aes(y = exp(as.numeric(median)), colour = exper_new), formula = y ~ s(x, k = time_max), method = "gam", se = FALSE) + 
      stat_smooth(aes(y = exp(as.numeric(lower)), colour = exper_new), formula = y ~ s(x, k = time_max), method = "gam", se = FALSE, linetype = "twodash") +
      stat_smooth(aes(y = exp(as.numeric(upper)), colour = exper_new), formula = y ~ s(x, k = time_max), method = "gam", se = FALSE, linetype = "twodash") +
      labs(x = time_title, y = "Hazard Ratio") +
      scale_color_brewer(palette = "Set1") +
      theme_bw() +  
      theme(legend.title = element_blank(), legend.position = "bottom")
  } else {
    ggplot(surv_data, aes(x = as.numeric(time))) + 
      stat_smooth(aes(y = as.numeric(median), colour = treat), formula = y ~ s(x, k = time_max), method = "gam", se = FALSE) + 
      stat_smooth(aes(y = as.numeric(lower), colour = treat), formula = y ~ s(x, k = time_max), method = "gam", se = FALSE, linetype = "twodash") +
      stat_smooth(aes(y = as.numeric(upper), colour = treat), formula = y ~ s(x, k = time_max), method = "gam", se = FALSE, linetype = "twodash") +
      labs(x = time_title, y = "Overall Survival") +
      scale_color_brewer(palette = "Set1") +
      theme_bw() +  
      theme(legend.title = element_blank(), legend.position = "bottom")
    
  }
}