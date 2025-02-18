# checker to check dimensions of all fields of d_maps
# input: sim_ci (list)
# output: d_maps (with proper dimensions)

checker <- function(d_maps, output_file = 'checked_d_maps') {
  
  output_path = file.path(intermediate_dir, paste0(output_file, '_', Sys.Date(), '.RData'))
  
  for (i in 1:length(d_maps)) {
    
    for (t in names(d_maps[[i]])) {

      d <- d_maps[[i]][[t]]$d
      sim_ci_lb <- unlist(d_maps[[i]][[t]]$sim_ci_lb)
      sim_ci_ub <- unlist(d_maps[[i]][[t]]$sim_ci_ub)
      if (grepl("motion.regression", t)) {
        d.fullres <- d_maps[[i]][[t]]$d.fullres
        sim_ci_lb.fullres <- d_maps[[i]][[t]]$sim_ci_lb.fullres
        sim_ci_ub.fullres <- d_maps[[i]][[t]]$sim_ci_ub.fullres
      }
      
      # transpose if needed
      
      if (dim(d)[1] > 1) {
        d_maps[[i]][[t]]$d <- t(d)
      }
      
      if (length(sim_ci_lb)[1] > 1) {
        d_maps[[i]][[t]]$sim_ci_lb <- t(sim_ci_lb)
      }
      
      if (length(sim_ci_ub)[1] > 1) {
        d_maps[[i]][[t]]$sim_ci_ub <- t(sim_ci_ub)
      }
      
      # repeat for regression case -  # TODO: could also simplify + combine w above
      if (grepl("motion.regression", t)) {
        
        if (dim(d.fullres)[1] > 1) {
          d_maps[[i]][[t]]$d.fullres <- t(d.fullres)
        }
        
        if (length(sim_ci_lb.fullres)[1] > 1) {
          d_maps[[i]][[t]]$sim_ci_lb.fullres <- t(sim_ci_lb.fullres)
        }
        
        if (length(sim_ci_ub.fullres)[1] > 1) {
          d_maps[[i]][[t]]$sim_ci_ub.fullres <- t(sim_ci_ub.fullres)
        }
      }
      
    }}
  
  save(d_maps, file = output_path)
  
  return(d_maps)
}
