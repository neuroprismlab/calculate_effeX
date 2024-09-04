# calculate simultaneous confidence intervals and add to d_maps

calc_sim_ci <- function(d_maps, alpha = 0.05, num_sdx_r2d = 2, output_file = 'sim_ci_data') {
  
  output_path = file.path(intermediate_dir, paste0(output_file, '_', Sys.Date(), '.RData'))
  
  for (i in 1:length(d_maps)) {
    
    for (t in names(d_maps[[i]])) {

      this_r <- d_maps[[i]][[t]]$b.standardized
      this_alpha_corrected <- alpha / length(this_r)
      this_n <- d_maps[[i]][[t]]$n[1]
      
      z_95 <- qnorm(1 - this_alpha_corrected) # e.g., 0.05 = 1.96
      
      r_ci_lb <- tanh(atanh(this_r) - z_95 / sqrt(this_n-3))
      r_ci_ub <- tanh(atanh(this_r) + z_95 / sqrt(this_n-3))
      
      d_ci_lb <- num_sdx_r2d * r_ci_lb / (1 - r_ci_lb ^ 2) ^ (1/2)
      d_ci_ub <- num_sdx_r2d * r_ci_ub / (1 - r_ci_ub ^ 2) ^ (1/2)
      
      # add sim CI to d
      d_maps[[i]][[t]]$sim_ci_lb <- d_ci_lb
      d_maps[[i]][[t]]$sim_ci_ub <- d_ci_ub
    }
    
  }
  
  # save results
  save(sim_ci = d_maps, file = output_path)
  
  # return d
  return(list(d_maps = d_maps))
}