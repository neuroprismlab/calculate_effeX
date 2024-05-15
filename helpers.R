# Helper functions for calculate_effex


####### Helpers for clean_data.R ########

## Triangle Helper Function:

# take in an effect map, 
# check if it's FC or act, if it's FC:
    # check if it's a triangle or square, 
    # convert to square if it's a triangle
    # return as a numeric vector (1D)
# if it's act: return the effect map as a vector

# input: an effect map as a list with d (either square or triangle).. e.g. effect_maps[1]
# output: an effect map as a vector (square) in a format that plots nicely as
# a symmetric matrix, or as a vector of the activation map if act

triangle_to_square <-  function(effect_map) {
  # check if the map is an FC study:
  if (grepl("fc", names(effect_map))) {
    # if the map is a full matrix (square):
    if (sqrt(length(effect_map[[1]]$d)) %% 1 == 0) {
      # do nothing
      square_map <- c(effect_map[[1]]$d)
    } else if ((((-1 + sqrt(1 + 8 * length(effect_map[[1]]$d))) / 2) + 1) %% 1 == 0) {
      # else if the map is half a matrix (triangle): 
      n_nodes <- ((-1 + sqrt(1 + 8 * length(effect_map[[1]]$d))) / 2) + 1
      # make a mask for the upper triangle (triumask)
      triumask <- upper.tri(matrix(1, nrow = n_nodes, ncol = n_nodes))
      # make a new matrix called upper that holds that upper mask
      upper <- triumask
      # fill the upper triangle with the map data
      upper[triumask] <- effect_map[[1]]$d
      # transpose that matrix so that upper now contains data in the lower triangle
      upper <- t(upper)
      # fill the upper triangle (now empty) with the map data
      upper[upper.tri(upper)] <- effect_map[[1]]$d
      square <- upper
      # transform matrix back to vector
      square_map <- c(square)
    }
    return(square_map)
  }
  else {
    return(c(effect_map[[1]]$d))
  }
}