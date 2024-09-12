plot_sim_ci <- function(data, study_name, study) {
  
  study_details <- study[study$name == study_name,]
  # remove na
  na_idx <- is.na(data$d) | is.na(data$sim_ci_lb) | is.na(data$sim_ci_ub)
  data$d <- data$d[!na_idx]
  data$sim_ci_lb <- data$sim_ci_lb[!na_idx]
  data$sim_ci_ub <- data$sim_ci_ub[!na_idx]
  # sort data from smallest to largest d
  sorted_indices <- order(data$d)
  sorted_d <- data$d[sorted_indices]
  # sort confidence intervals by the same order
  sorted_upper_bounds <- data$sim_ci_ub[sorted_indices]
  sorted_lower_bounds <- data$sim_ci_lb[sorted_indices]
  
  # downsample data for plotting
  downsample <- length(sorted_indices) %/% 100
  if (downsample == 0) {
    downsample = 1
  }
  sorted_d <- sorted_d[seq(1, length(sorted_d), by = downsample)]
  sorted_upper_bounds <- sorted_upper_bounds[seq(1, length(sorted_upper_bounds), by = downsample)]
  sorted_lower_bounds <- sorted_lower_bounds[seq(1, length(sorted_lower_bounds), by = downsample)]
  
  
  # for coloring of confidence intervals:
  below_zero <- sorted_upper_bounds < 0
  below_cross_idx <- which(diff(below_zero) == -1) # the last TRUE before switch
  
  above_zero <- sorted_lower_bounds > 0
  above_cross_idx <- (which(diff(above_zero) == 1)) + 1 # the last FALSE before switch to true
  
  n_title <- paste0("n = ", data$n)
  
  # calculate the percent of edges/voxels with confidence intervals that don't overlap with zero:
  percent_below_zero <- sum(sorted_upper_bounds < 0) / length(sorted_upper_bounds)
  percent_above_zero <- sum(sorted_lower_bounds > 0) / length(sorted_lower_bounds)
  
  
  # if there are no values below zero, set the index to 1
  if (length(below_cross_idx) == 0) {
    below_cross_idx = 1
  } 
  
  # if there are no values above zero, set the index to the end
  if (length(above_cross_idx) == 0) {
    above_cross_idx = length(above_zero)
  } 
  
  # plot a line for d
  par(mar=c(2, 4, 5, 2))
  plot(sorted_d, type = "l", ylim = c(min(sorted_lower_bounds, na.rm = TRUE), max(sorted_upper_bounds, na.rm = TRUE)),
       xlab = "Edges/Voxels", ylab = "Cohen's d", axes = FALSE)
  # add a horizontal line at y = 0
  abline(h = 0, col = "#ba2d25", lty = 3)
  axis(2, las = 1)  # Add left axis with labels parallel to the axis (las = 1)
  legend("topleft", inset = c(-0.1, -0.5),
         legend = c(
           bquote(bold("Dataset:")), 
           paste(study_details[,"dataset"], "  "),
           bquote(bold("Map Type:")), 
           paste(study_details$map_type, "  "),
           bquote(bold("Test type:")), 
           paste(study_details$orig_stat_type, "  "),
           bquote(bold("Component1:")), 
           paste(study_details$test_component_1, "  "),
           bquote(bold("Component2:")), 
           paste(study_details$test_component_2, "  "),
           bquote(bold("Sample Size:")),
           paste(n_title)
         ), 
         bty = "n", ncol = 6, cex = 1, text.width = c(15, 15, 15, 15, 25, 20), x.intersp = 0.0, xpd = TRUE)
  legend("bottomright", inset = c(0, -0.2), legend = c(bquote(bold("Maximum conservative effect size: ")), 
                                                       ifelse((abs(max(data$sim_ci_lb, na.rm = TRUE)) > abs(min(data$sim_ci_ub, na.rm = TRUE))), 
                                                              ifelse((max(data$sim_ci_lb, na.rm = TRUE) > 0),
                                                                     round(abs(max(data$sim_ci_lb, na.rm = TRUE)), 2), 0),
                                                              ifelse((min(data$sim_ci_ub, na.rm = TRUE) < 0), round(abs(min(data$sim_ci_ub, na.rm = TRUE)), 2), 0))), xjust = 1, yjust = 1, col = 2, bty = "n", cex = 1, x.intersp = 0, xpd = TRUE)
  
  
  # plot and shade the cofidence intervals:
  # green for intervals that are entirely below zero
  polygon(c(1:below_cross_idx, rev(1:below_cross_idx)), 
          c(sorted_upper_bounds[1:below_cross_idx], rev(sorted_lower_bounds[1:below_cross_idx])), 
          col = rgb(177/255, 207/255, 192/255, alpha = 0.5), border = NA)
  
  # red for intervals that include zero
  polygon(c(below_cross_idx:above_cross_idx, rev(below_cross_idx:above_cross_idx)), 
          c(sorted_upper_bounds[below_cross_idx:above_cross_idx], rev(sorted_lower_bounds[below_cross_idx:above_cross_idx])), 
          col = rgb(237/255, 185/255, 185/255, alpha = 0.5), border = NA)
  
  # green for intervals that are entirely above zero
  polygon(c(above_cross_idx:length(above_zero), rev(above_cross_idx:length(above_zero))), 
          c(sorted_upper_bounds[above_cross_idx:length(above_zero)], rev(sorted_lower_bounds[above_cross_idx:length(above_zero)])), 
          col = rgb(177/255, 207/255, 192/255, alpha = 0.5), border = NA)
}
