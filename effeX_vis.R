# plot simultaneous confience interval plots
# TODO: check why the confidence interval for some studies (e.g. ukb fc r rest memory) doesn't include d at the ends
# TODO: shade the intervals red if they include zero (i.e. not significant)

plot_sim_ci <- function(v_d_clean, i) {
  data <- v_d_clean[[i]]
  # remove na
  na_idx <- is.na(data$d) | is.na(data$sim_ci_lb) | is.na(data$sim_ci_ub)
  data$d <- data$d[!na_idx]
  data$sim_ci_lb <- data$sim_ci_lb[!na_idx]
  data$sim_ci_ub <- data$sim_ci_ub[!na_idx]
  sorted_indices <- order(data$d)
  sorted_d <- data$d[sorted_indices]
  sorted_upper_bounds <- data$sim_ci_ub[sorted_indices]
  sorted_lower_bounds <- data$sim_ci_lb[sorted_indices]
  
  plot(sorted_d, type = "l", ylim = c(min(sorted_lower_bounds, na.rm = TRUE), max(sorted_upper_bounds, na.rm = TRUE)),
       main = names(v_d_clean[i]), xlab = "Edges/Voxels", ylab = "Cohen's d")
  
  polygon(c(1:length(sorted_d), rev(1:length(sorted_d))), 
          c(sorted_upper_bounds, rev(sorted_lower_bounds)), 
          col = rgb(0.5, 0.5, 0.5, alpha = 0.3), border = NA)
}

