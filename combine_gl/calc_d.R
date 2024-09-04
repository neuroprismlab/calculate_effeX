### NEW CALC D SCRIPT

##############################################
#
# Convert effect maps to d
#
# In: 
# - cleaned_data$study: table with study attributes
# - cleaned_data$data: list of lists with effect maps
#
# Output: 
# - d_maps: list with effect maps converted to d, stored as attribute "d"
#
##############################################

calc_d <- function(study, data, num_sdx_r2d = 2, output_file = 'd_maps') {
  # num_sdx_r2d: number of standard deviations in X to use for Maya's r-to-d conversion
  
  output_path = file.path(intermediate_dir, paste0(output_file, '_', Sys.Date(), '.RData'))
  
  d_maps <- data
  
  # loop through each effect map
  for (i in 1:length(data)) {
    effect_name <- names(data)[i]
    study_idx <- which(study$name == effect_name)

    # for each analysis type within this study, convert d
    for (t in names(data[[i]])) {
      # convert to d
      this_r <- data[[i]][[t]]$b.standardized
      this_d <- num_sdx_r2d * this_r / ((1 - this_r^2) ^ (1/2))
      # add d to results list
      d_maps[[i]][[t]]$d <- this_d
    }
  }
  
  save(d_maps, file = output_path)
  return(d_maps)
}
