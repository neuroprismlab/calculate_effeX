# checker to check dimensions of all fields of d_maps
# input: sim_ci (list)
# output: d_maps (with proper dimensions)

checker <- function(d_maps, output_file = 'checked_d_maps') {
  
  output_path = file.path(intermediate_dir, paste0(output_file, '_', Sys.Date(), '.RData'))
  
  for (i in 1:length(d_maps)) {
    
    for (t in names(d_maps[[i]])) {
      
      d <- d_maps[[i]][[t]]$d
      sim_ci_lb <- d_maps[[i]][[t]]$sim_ci_lb
      sim_ci_ub <- d_maps[[i]][[t]]$sim_ci_ub
      
      if (dim(d)[1] > 1) {
        d_maps[[i]][[t]]$d <- t(d)
      }
      
      if (length(sim_ci_lb)[1] > 1) {
        d_maps[[i]][[t]]$sim_ci_lb <- t(sim_ci_lb)
      }
      
      if (length(sim_ci_ub)[1] > 1) {
        d_maps[[i]][[t]]$sim_ci_ub <- t(sim_ci_ub)
      }
    }}
  
  save(d_maps, file = output_path)
  
  return(d_maps)
}
