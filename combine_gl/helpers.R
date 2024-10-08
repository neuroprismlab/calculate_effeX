# Helper functions for calculate_effex


####### Helpers for clean_data.R ########

## Triangle Helper Function:

# take in an effect map, 
# check if it's FC or act, if it's FC:
    # check if it's a triangle or square, 
    # convert to square if it's a triangle
    # return as a numeric vector (1D)
# if it's act: return the effect map as a vector

# input: an effect map as a list with d and orig_stat (either square or triangle).. e.g. effect_maps[1]
# output: an effect map as a vector (square) in a format that plots nicely as
# a symmetric matrix, or as a vector of the activation map if act


triangle_to_square <-  function(effect_map) {
    
    # check if the map is an FC study:
    if (grepl("fc", names(effect_map))) {

        # for d
        # if d is a full matrix (square):
        if (sqrt(length(effect_map[[1]]$d)) %% 1 == 0) {
            # do nothing
            square_d <- c(effect_map[[1]]$d)
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
            square_d <- c(square)
        }

        # for orig_stat
        # if the map is a full matrix (square):
        if (sqrt(length(effect_map[[1]]$orig_stat)) %% 1 == 0) {
            # do nothing
            square_orig_stat <- c(effect_map[[1]]$orig_stat)
        } else if ((((-1 + sqrt(1 + 8 * length(effect_map[[1]]$orig_stat))) / 2) + 1) %% 1 == 0) {
            # else if the map is half a matrix (triangle): 
            n_nodes <- ((-1 + sqrt(1 + 8 * length(effect_map[[1]]$orig_stat))) / 2) + 1
            # make a mask for the upper triangle (triumask)
            triumask <- upper.tri(matrix(1, nrow = n_nodes, ncol = n_nodes))
            # make a new matrix called upper that holds that upper mask
            upper <- triumask
            # fill the upper triangle with the map data
            upper[triumask] <- effect_map[[1]]$orig_stat
            # transpose that matrix so that upper now contains data in the lower triangle
            upper <- t(upper)
            # fill the upper triangle (now empty) with the map data
            upper[upper.tri(upper)] <- effect_map[[1]]$orig_stat
            square <- upper
            # transform matrix back to vector
            square_orig_stat <- c(square)
        }

        return(list(square_d, square_orig_stat))
    }

    # if the map is an activation study, just return the maps as vectors (1D)
    else {
        return(list(c(effect_map[[1]]$d), c(effect_map[[1]]$orig_stat)))
    }
}



########### HELPERS FOR CHANGING A SQUARE MATRIX TO A TRIANGLE MATRIX ###########

# FIRST NEED TO REORDER ACCORDING TO MAPPING IF APPLICABLE, THEN TAKE LOWER TRIANGLE, THEN PLOT

# 1. Convert to square matrix if not already and reorder by mapping if applicable
mapping <- function(data, map_path = NA) {

    # takes data (e.g. effect_map[[1]]$orig_stat) (COULD BE A SQUARE OR A TRIANGLE) and the path to a mapping file (e.g. 268 note network mapping)
    # returns a square matrix that is reorganized according to the map provided (if provided)
    # if it's a triangle...
    if (sqrt(length(data)) %% 1 != 0) {
        nrow = ((-1 + sqrt(1 + 8 * length(data))) / 2) + 1
        mat <- matrix(0, nrow = nrow, ncol = nrow)
        mat[upper.tri(mat)] <- data
        # reflect the triangle to get a filled in square matrix (so that when we reorganize it it doesn't get messed up)
        mat <- mat + t(mat)
        # now mat is a full square matrix 

    } else { # if it's a square...
        nrow = sqrt(length(data))
        # turn effect map into a square matrix
        mat <- matrix(data = data, nrow = nrow)
        # set the diagonal to zero
        diag(mat) <- 0

    }

    # if map is provided:
    if (!is.na(map_path)) {
        # load map
        mapping <- read.csv(map_path, header = TRUE)

        # reorder based on mapping
        ordered_mat <- mat[mapping$oldroi, mapping$oldroi]
    } else {
        # if no map provided, don't change the order
        ordered_mat <- mat
    }
    
    # return the ordered matrix (or original if no map provided) as a square matrix with zeros in diagonal
    return(ordered_mat)
} 

# 2. Take lower triangle
lower_triangle <- function(ordered_mat) {

    # make a mask for the lower triangle (trilmask)
    trilmask <- lower.tri(matrix(1, nrow = nrow(ordered_mat), ncol = ncol(ordered_mat)), diag = TRUE)
    # set the lower triangle to NA
    ordered_mat[trilmask] <- NA
    
    return(ordered_mat)
}

# 3. Plot the lower triangle 
plot_lower_triangle <- function(triangle_ordered, map_path = NA) {
    
    # triangle_ordered to matrix
    # triangle_ordered <- matrix(triangle_ordered, nrow = sqrt(length(triangle_ordered)))

    # melt the matrix for ggplot
    triangle_melted <- melt(triangle_ordered)
    colnames(triangle_melted) <- c("Var1", "Var2", "value")

    heatmap_plot <- ggplot(triangle_melted, aes(Var1, Var2, fill = value)) +
      
      geom_tile() +
      scale_fill_gradient2(limits = c(min(triangle_melted$value), max(triangle_melted$value)),
        low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme_minimal() +
      theme(axis.title.x = element_text(margin = margin(t = 10)),
            axis.title.y = element_text(margin = margin(r = 10)),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.margin = margin(.5, .5, .5, .5, "lines"),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

    if (!is.na(map_path)) {
        # load mapping
        mapping <- read.csv(map_path, header = TRUE)

        for (i in 1:(nrow(mapping) - 1)) {
            if (mapping$category[i] != mapping$category[i + 1]) {
            heatmap_plot <- heatmap_plot + geom_vline(xintercept = i, color = "black", size = 0.3) +
                geom_hline(yintercept = i, color = "black")
            }
        }
        
        # Calculate the positions of the labels
        label_positions <- c(1, which(mapping$category[-1] != mapping$category[-length(mapping$category)]) + 1, length(mapping$category) + 1)
        label_positions <- (label_positions[-1] + label_positions[-length(label_positions)]) / 2
        label_strings <- mapping$label[label_positions]
        
        # Add labels to each mapping category
        heatmap_plot <- heatmap_plot + annotate("text", x = label_positions, y = -6, label = label_strings, angle = 90, hjust = 1, vjust=0.5, size=3.5) + coord_cartesian(clip="off")
        heatmap_plot <- heatmap_plot + annotate("text", x = -10, y = label_positions, label = label_strings, angle = 0, hjust = 0.5, vjust=1, size=3.5)

        # Add axis labels to the heatmap
        heatmap_plot <- heatmap_plot + labs(x = "Network", y = "Network")
        
    }
    print(heatmap_plot)
}

get_triangle <- function(triangle_ordered) {
    return(triangle_ordered[!is.na(triangle_ordered)])
}

# 4. Combine all functions
# returns a triangle as a vector without NA values, so just the one triangle
square_to_triangle <- function(effect_map, map_path = NA, show_plot = TRUE) {
    ordered_mat <- mapping(effect_map, map_path)
    triangle_ordered <- lower_triangle(ordered_mat)
    if (show_plot == TRUE) {
        plot_lower_triangle(triangle_ordered, map_path)
    }
    result <- get_triangle(triangle_ordered)
    return(result) #TODO: save the result as a file
}






#######e From a triangle, plot a full square matrix
# Useful for QC!

plot_full_mat <- function(triangle_ordered, mapping_path, export = FALSE, export_path = NA, show_plot = FALSE) {
    # takes an ordered triangle vector (without NAs) and plots the full matrix
    if (export == TRUE) {png(export_path)}
    
    nrow = (((-1 + sqrt(1 + 8 * length(triangle_ordered))) / 2) + 1)

    # mirror the triangle across the x = y line to get full matrix
    # first fill in half the matrix with the triangle data
    mat <- matrix(0, nrow = nrow, ncol = nrow)
    mat[upper.tri(mat)] <- triangle_ordered
    full_mat <- mat + t(mat) #- diag(diag(triangle_ordered))

    # melt the matrix for ggplot
    melted <- melt(full_mat)
    colnames(melted) <- c("Var1", "Var2", "value")

    heatmap_plot <- ggplot(melted, aes(Var1, Var2, fill = value)) +
      
      geom_tile() +
      scale_fill_gradient2(limits = c(min(melted$value), max(melted$value)),
        low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme_minimal() +
      theme(axis.title.x = element_text(margin = margin(t = 10)),
            axis.title.y = element_text(margin = margin(r = 10)),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.margin = margin(.5, .5, .5, .5, "lines"),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

    if (!is.na(mapping_path)) {
        # load mapping
        mapping <- read.csv(mapping_path, header = TRUE)

        
        for (i in 1:(nrow(mapping) - 1)) {
            if (mapping$category[i] != mapping$category[i + 1]) {
            heatmap_plot <- heatmap_plot + geom_vline(xintercept = i, color = "black", size = 0.3) +
                geom_hline(yintercept = i, color = "black")
            }
        }
        
        # Calculate the positions of the labels
        label_positions <- c(1, which(mapping$category[-1] != mapping$category[-length(mapping$category)]) + 1, length(mapping$category) + 1)
        label_positions <- (label_positions[-1] + label_positions[-length(label_positions)]) / 2
        label_strings <- mapping$label[label_positions]
        
        # Add labels to each mapping category
        heatmap_plot <- heatmap_plot + annotate("text", x = label_positions, y = -6, label = label_strings, angle = 90, hjust = 1, vjust=0.5, size=3.5) + coord_cartesian(clip="off")
        heatmap_plot <- heatmap_plot + annotate("text", x = -10, y = label_positions, label = label_strings, angle = 0, hjust = 0.5, vjust=1, size=3.5)

        # Add axis labels to the heatmap
        heatmap_plot <- heatmap_plot + labs(x = "Network", y = "Network")
    }
    
    if (show_plot == TRUE) {
        print(heatmap_plot)
    }
    
    # export plot as png if export = TRUE
    if (export == TRUE) {
        ggsave(export_path, plot = heatmap_plot, device = "png", width = 10, height = 10, dpi = 200)
        dev.off()
    }
}



##### Helpers for adding phenotypic categories to study dataframe

# add phenotypic categories to study dataframe

add_phen <- function(study, effect_maps, phen_file = "/work/neuroprism/effect_size/data/helper_data/phen.csv") {
    # load phenotypic data file (phen_file) from data directory
    phen <- read.csv(phen_file, header = TRUE)
    # merge phenotypic data with study data
    # TODO: check to see if the data merges properly, and if the study names are in the same format in both (e.g. capitalization, _ vs. .)
    phen_study <- left_join(study, phen, by = "name")

    for (i in 1:dim(phen_study)[1]) {
        name <- phen_study[i, "name"]
        len <- length(effect_maps[[which(toupper(names(effect_maps)) == name)]]$orig_stat)
        if (phen_study[i,"map_type"] == "FC") {
            if ((((-1 + sqrt(1 + 8 * len)) / 2) + 1) == 55) {
            ref <- "UKB_55"
            }
            else if ((((-1 + sqrt(1 + 8 * len)) / 2) + 1) == 268) {
            ref <- "Shen_268"
            }
            else {
            stop(paste0("Unknown parcellation found, please add this parcellation. Length of effect map: ", len, ". study name: ", name))
            }
        }
        
        else if (phen_study[i, "map_type"] == "ACT") {
            ref <- "Voxel"
        }
        # add ref column
        phen_study$ref[i] <- ref
    }

    return(phen_study)
}


#########################
# NEW HELPER FUNCTIONS FOR NEW DATA FORMAT 03/09/2024

# Function to take a matlab struct file loaded into R with readMat and properly
# name the struct fields

library(R.matlab)

# very special correction for readMat bug with nested lists
# WARNING: sadly, this is very idiosyncratic to the way the data was saved
# Unfortunately we aren't authors of the original readMat function, but it'd
# be nice to have a more elegant workaround
readMat_act_correction <- function(mat_struct) {
  mat_struct[[1]][[1]][[7]][[24]] <- mat_struct[[1]][[1]][[8]]
  mat_struct[[1]][[1]][[8]] <- mat_struct[[1]][[2]]
  mat_struct[[1]][[2]] <- mat_struct[[2]]
  mat_struct[[3]] <- NULL
  mat_struct[[2]] <- NULL
  return(mat_struct)
}

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

#######
# QC plot that takes triangle and brain mask


plot_mat_w_mask <- function(triangle, mask, mapping_path = "/work/neuroprism/effect_size/data/helper_data/map268_subnetwork.csv", export = FALSE, export_path = NA, show_plot = TRUE) {
  # takes an ordered triangle vector (without NAs) and plots the full matrix
  if (export == TRUE) {png(export_path)}
  
  nrow = (((-1 + sqrt(1 + 8 * length(triangle))) / 2) + 1)
  
  # mirror the triangle across the x = y line to get full matrix
  # first fill in half the matrix with the triangle data
  full_mat <- mask
  full_mat[mask==1] <- triangle
  
  # melt the matrix for ggplot
  melted <- melt(full_mat)
  colnames(melted) <- c("Var1", "Var2", "value")
  
  heatmap_plot <- ggplot(melted, aes(Var1, Var2, fill = value)) +
    
    geom_tile() +
    scale_fill_gradient2(limits = c(min(melted$value), max(melted$value)),
                         low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    theme(axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(.5, .5, .5, .5, "lines"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  
  if (!is.na(mapping_path)) {
    # load mapping
    mapping <- read.csv(mapping_path, header = TRUE)
    
    
    for (i in 1:(nrow(mapping) - 1)) {
      if (mapping$category[i] != mapping$category[i + 1]) {
        heatmap_plot <- heatmap_plot + geom_vline(xintercept = i, color = "black", size = 0.3) +
          geom_hline(yintercept = i, color = "black")
      }
    }
    
    # Calculate the positions of the labels
    label_positions <- c(1, which(mapping$category[-1] != mapping$category[-length(mapping$category)]) + 1, length(mapping$category) + 1)
    label_positions <- (label_positions[-1] + label_positions[-length(label_positions)]) / 2
    label_strings <- mapping$label[label_positions]
    
    # Add labels to each mapping category
    heatmap_plot <- heatmap_plot + annotate("text", x = label_positions, y = -6, label = label_strings, angle = 90, hjust = 1, vjust=0.5, size=3.5) + coord_cartesian(clip="off")
    heatmap_plot <- heatmap_plot + annotate("text", x = -10, y = label_positions, label = label_strings, angle = 0, hjust = 0.5, vjust=1, size=3.5)
    
    # Add axis labels to the heatmap
    heatmap_plot <- heatmap_plot + labs(x = "Network", y = "Network")
  }
  
  if (show_plot == TRUE) {
    print(heatmap_plot)
  }
  
  # export plot as png if export = TRUE
  if (export == TRUE) {
    ggsave(export_path, plot = heatmap_plot, device = "png", width = 10, height = 10, dpi = 200)
    dev.off()
  }
}