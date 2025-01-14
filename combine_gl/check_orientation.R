# check_orientation

# takes data and brain_masks
# loops through each FC study, check if that study's brain mask is upper triangle
# if it is, leave it, if it's lower triangle, convert to upper,
# if it's neither, warn the user

# check if every study has a brain mask, give warnings for which ones don't
# TODO:

check_orientation <- function(data, brain_masks) {

  # warn of any studies that don't have brain masks
  for (d in names(data)) {
    if (!(tolower(d) %in% tolower(names(brain_masks)))) {
      print(paste0("the study ", d, " does not have a brain mask and therefore the orientation of the data cannot be checked. Please add a brain mask."))
    }
  }

  for (i in seq(1, length(brain_masks))) {
    
    if (grepl('_fc_', names(brain_masks[i]))) { # for FC only
      if (all(brain_masks[[i]]$mask[upper.tri(brain_masks[[i]]$mask)] == 0)) {
        print(paste0('the study ', names(brain_masks[i]), ' uses lower triangle'))
        
        # for the entries of data matching this study, flip all FC-shaped vars
        for (v in names(data[[names(brain_masks[i])]])) {
          for (a in names(data[[names(brain_masks[i])]][[v]])) {
            if (length(data[[names(brain_masks[i])]][[v]][[a]]) == sum(brain_masks[[i]]$mask == 1)) {
              print(paste0('flipping triangle for ', names(data[[names(brain_masks[i])]][[v]][a])))
              data[[names(brain_masks[i])]][[v]][[a]] <- lower_to_upper_triangle(data[[names(brain_masks[i])]][[v]][[a]])
            }
          }
        }
      }
    }
  }
  return(data)
}