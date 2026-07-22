# checker to check dimensions of all fields of d_maps
# input: sim_ci (list)
# output: d_maps (with proper dimensions)
checker <- function(d_maps, output_file = 'checked_d_maps', int_dir = intermediate_dir) {
  
  output_path = file.path(int_dir, paste0(output_file, '_', Sys.Date(), '.RData'))
  
  for (i in 1:length(d_maps)) {
    
    for (t in names(d_maps[[i]])) {

      # d
      d <- d_maps[[i]][[t]]$d
      sim_ci_lb <- unlist(d_maps[[i]][[t]]$sim_ci_lb)
      sim_ci_ub <- unlist(d_maps[[i]][[t]]$sim_ci_ub)
      if (grepl("motion.regression", t)) {
        d.fullres <- d_maps[[i]][[t]]$d.fullres
        sim_ci_lb.fullres <- d_maps[[i]][[t]]$sim_ci_lb.fullres
        sim_ci_ub.fullres <- d_maps[[i]][[t]]$sim_ci_ub.fullres
      }
      
      # transpose if needed
      
      if (!is.null(dim(d))) {
        if (dim(d)[1] > 1) {
          d_maps[[i]][[t]]$d <- t(d)
        }
      }
      
      # FIX: added is.null guard so scalar sim_ci_lb/ub (e.g. from multivariate tests) don't error
      if (!is.null(dim(sim_ci_lb)) && dim(sim_ci_lb)[1] > 1) {
        d_maps[[i]][[t]]$sim_ci_lb <- t(sim_ci_lb)
      }
      
      if (!is.null(dim(sim_ci_ub)) && dim(sim_ci_ub)[1] > 1) {
        d_maps[[i]][[t]]$sim_ci_ub <- t(sim_ci_ub)
      }

      # r_sq
      r_sq <- d_maps[[i]][[t]]$r_sq
      r_sq_sim_ci_lb <- unlist(d_maps[[i]][[t]]$r_sq_sim_ci_lb)
      r_sq_sim_ci_ub <- unlist(d_maps[[i]][[t]]$r_sq_sim_ci_ub)
      if (grepl("motion.regression", t)) {
        r_sq.fullres <- d_maps[[i]][[t]]$r_sq.fullres
        r_sq_sim_ci_lb.fullres <- d_maps[[i]][[t]]$r_sq_sim_ci_lb.fullres
        r_sq_sim_ci_ub.fullres <- d_maps[[i]][[t]]$r_sq_sim_ci_ub.fullres
      }

      if (!is.null(dim(r_sq))) {
        if (dim(r_sq)[1] > 1) {
          d_maps[[i]][[t]]$r_sq <- t(r_sq)
        }
      }

      # FIX: added is.null guard
      if (!is.null(dim(r_sq_sim_ci_lb)) && dim(r_sq_sim_ci_lb)[1] > 1) {
        d_maps[[i]][[t]]$r_sq_sim_ci_lb <- t(r_sq_sim_ci_lb)
      }
      
      if (!is.null(dim(r_sq_sim_ci_ub)) && dim(r_sq_sim_ci_ub)[1] > 1) {
        d_maps[[i]][[t]]$r_sq_sim_ci_ub <- t(r_sq_sim_ci_ub)
      }
      
      # repeat for regression case -  # TODO: could also simplify + combine w above
      if (grepl("motion.regression", t)) {
        
        if (!is.null(dim(d.fullres))) {
          if (dim(d.fullres)[1] > 1) {
            d_maps[[i]][[t]]$d.fullres <- t(d.fullres)
          }
        }
        
        # FIX: added is.null guard
        if (!is.null(dim(sim_ci_lb.fullres)) && dim(sim_ci_lb.fullres)[1] > 1) {
          d_maps[[i]][[t]]$sim_ci_lb.fullres <- t(sim_ci_lb.fullres)
        }
        
        if (!is.null(dim(sim_ci_ub.fullres)) && dim(sim_ci_ub.fullres)[1] > 1) {
          d_maps[[i]][[t]]$sim_ci_ub.fullres <- t(sim_ci_ub.fullres)
        }

        # r_sq
        if (!is.null(dim(r_sq.fullres))) {
          if (dim(r_sq.fullres)[1] > 1) {
            d_maps[[i]][[t]]$r_sq.fullres <- t(r_sq.fullres)
          }
        }
        
        # FIX: added is.null guard
        if (!is.null(dim(r_sq_sim_ci_lb.fullres)) && dim(r_sq_sim_ci_lb.fullres)[1] > 1) {
          d_maps[[i]][[t]]$r_sq_sim_ci_lb.fullres <- t(r_sq_sim_ci_lb.fullres)
        }
        
        if (!is.null(dim(r_sq_sim_ci_ub.fullres)) && dim(r_sq_sim_ci_ub.fullres)[1] > 1) {
          d_maps[[i]][[t]]$r_sq_sim_ci_ub.fullres <- t(r_sq_sim_ci_ub.fullres)
        }
      }
      
    }}
  
  save(d_maps, file = output_path)
  
  return(d_maps)
}