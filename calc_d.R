##############################################
#
# Convert effect maps to d
#
# In: 
# - cleaned_data$study: table with study attributes
# - cleaned_data$effect_map: list with effect maps
#
# Output: 
# - effect_map_d: list with effect maps converted to d
#
##############################################


## IN PROGRESS!!! ## Have not tested at all yet since Discovery is down



# # Convert effect maps to d

 calc_d <- function(study, effect_map) {
  

  # loop through each effect map
  for (i in 1:length(effect_map)) {
    # Check if effect map is a t value by checking if study$orig_stat_type is equal to "t"
    if (study$orig_stat_type[i] == "t") {
      # convert to d
      
      # add d to results list

    }
    # check if effect map is a correlation by checking if study$orig_stat_type is equal to "r"
    else if (study$orig_stat_type[i] == "r") {
      # convert to d

      # add d to results list
    }
    # check if effect map is a d value by checking if study$orig_stat_type is equal to "d"
    else if (study$orig_stat_type[i] == "d") {
      # Do nothing

      # add d to results list
    }
    # check if effect map is a t2 value by checking if study$orig_stat_type is equal to "t2"
    else if (study$orig_stat_type[i] == "t2") {
      # convert to d

      # add d to results list
    }
  }
  
    return(results)

 }