##############################################
#
# Convert effect maps to d
#
# In: 
# - cleaned_data$study: table with study attributes
# - cleaned_data$effect_map: list with effect maps
#
# Output: 
# - effect_map: list with effect maps converted to d, stored as attribute "d"
#
##############################################

 calc_d <- function(study, effect_map, num_sdx_r2d = 2) {
  # num_sdx_r2d: number of standard deviations in X to use for Maya's r-to-d conversion
  
  # TODO: MAKE SURE THAT THE INDEXING IS CORRECT FOR STUDY VS. EFFECT_MAP!!!

  # loop through each effect map
  for (i in 1:length(effect_map)) {
    effect_name <- toupper(names(effect_map))[i]
    study_idx <- which(toupper(study$name) == effect_name)
    # Check if effect map is a t value by checking if study$orig_stat_type is equal to "t"
    if (study$orig_stat_type[study_idx] == "t") {
      # convert to d
      this_t <- effect_map[[i]]$orig_stat
      this_n <- effect_map[[i]]$n[1]
      this_n_groups <- 1
      this_d <- this_t / sqrt(this_n)
      # add d to effect_maps_d
      effect_map[[i]]$d <- this_d
    }
    # check if effect map is a correlation by checking if study$orig_stat_type is equal to "r"
    else if (study$orig_stat_type[study_idx] == "r") {
      # convert to d
      this_r <- effect_map[[i]]$orig_stat
      this_n <- effect_map[[i]]$n[1]
      this_n_groups <- 1
      this_d <- num_sdx_r2d * this_r / ((1 - this_r^2) ^ (1/2))
      # add d to results list
      effect_map[[i]]$d <- this_d
    }
    # check if effect map is a d value by checking if study$orig_stat_type is equal to "d"
    else if (study$orig_stat_type[study_idx] == "d") {
      # Do nothing
      this_d <- effect_map[[i]]$orig_stat
      # add d to results list
      effect_map[[i]]$d <- this_d
    }
    # check if effect map is a t2 value by checking if study$orig_stat_type is equal to "t2"
    else if (study$orig_stat_type[study_idx] == "t2") {
      # convert to d
      this_t2 <- effect_map[[i]]$orig_stat
      this_n1 <- effect_map[[i]]$n1[1]
      this_n2 <- effect_map[[i]]$n2[1]
      this_n <- this_n1 + this_n2
      this_n_groups <- 2
      this_d <- this_t2 * sqrt(1/this_n1 + 1/this_n2)

      # add d to results list
      effect_map[[i]]$d <- this_d
    }
    else {
      print("Error: could not convert to d. Check that orig_stat_type is one of: r, t, d, t2.")
    }
  }


    # # run all FC maps through triangle_to_square() function to ensure they are in square format
    # TODO: this only converts d and orig_stat to square now, not std or p
    # for (i in 1:length(effect_map)) {
    #     effect_map[[i]]$d <- triangle_to_square(effect_map[i])[[1]]
    #     effect_map[[i]]$orig_stat <- triangle_to_square(effect_map[i])[[2]]
    # }


  
    return(list(study = study, effect_map = effect_map))

 }