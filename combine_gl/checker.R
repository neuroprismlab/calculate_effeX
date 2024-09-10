# checker to check dimensions of all fields of data
# input: sim_ci (list)
# output: data (with proper dimensions)

checker <- function(sim_ci, output_file = 'checked_data') {
  
  output_path = file.path(intermediate_dir, paste0(output_file, '_', Sys.Date(), '.RData'))
  
  data <- sim_ci$d_maps
  
  for (i in 1:length(data)) {
    
    for (t in names(data[[i]])) {
      
      d <- data[[i]][[t]]$d
      sim_ci_lb <- data[[i]][[t]]$sim_ci_lb
      sim_ci_ub <- data[[i]][[t]]$sim_ci_ub
      
      if (dim(d)[1] > 1) {
        data[[i]][[t]]$d <- t(d)
      }
      
      if (dim(sim_ci_lb)[1] > 1) {
        data[[i]][[t]]$sim_ci_lb <- t(sim_ci_lb)
      }
      
      if (dim(sim_ci_ub)[1] > 1) {
        data[[i]][[t]]$sim_ci_ub <- t(sim_ci_ub)
      }
    }}
  
  save(data, file = output_path)
  
  return(data)
}