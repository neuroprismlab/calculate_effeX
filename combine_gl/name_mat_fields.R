# Function to take a matlab struct file loaded into R with readMat and properly
# name the struct fields

library(R.matlab)

# Load the .mat file
path <- '/work/neuroprism/effect_size/data/group_level/29-Aug-2024abcd_fc_t2_rest_sex.mat'
mat_data <- readMat(path)

# Define a function to automatically assign names to the sublists
assign_names <- function(mat_list) {
  # Extract the top-level names
  main_names <- names(mat_list)
  
  # Loop over each element in the main list
  for (i in seq_along(mat_list)) {
    sublist <- mat_list[[i]]
    
    # Check if the sublist is a list and contains further nested lists (MATLAB structs)
    if (is.list(sublist)) {
      # Extract the field names from the attributes (different method)
      field_names <- names(sublist)
      
      # If the field names are missing, try to extract them from the attributes
      if (is.null(field_names)) {
        field_names <- attributes(sublist)$dimnames[[1]]
      }
      
      # If field names are found, assign them to the sublist
      if (!is.null(field_names)) {
        names(mat_list[[i]]) <- field_names
      }
      
      # Recursively apply the function to any nested lists
      mat_list[[i]] <- assign_names(mat_list[[i]])
    }
  }
  
  # Ensure the main list is named
  if (!is.null(main_names)) {
    names(mat_list) <- main_names
  }
  
  return(mat_list)
}

# Apply the function to your loaded data
results_named <- assign_names(mat_data)

# Inspect the results
str(results_named)

