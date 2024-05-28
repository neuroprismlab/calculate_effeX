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


####### Helpers for square_to_triangle ########

# take in an effect map,
# check if it's FC or act, if it's FC:
    # check if it's a triangle or square, 
    # convert to triangle if it's a square
    # return as a numeric vector (1D)
# if it's act: return the effect map as a vector

# input: an effect map as a list with d and orig_stat (either square or triangle).. e.g. effect_maps[1]
# output: an effect map as a vector (triangle)

# square_to_triangle <-  function(effect_map) {
    
#     # check if the map is an FC study:
#     if (grepl("fc", names(effect_map))) {

#         # for d
#         # if d is a full matrix (square):
#         if (sqrt(length(effect_map[[1]]$orig_stat)) %% 1 == 0) {
#             square_d <- matrix(effect_map[[1]]$orig_stat, nrow = sqrt(length(effect_map[[1]]$orig_stat)))
#             # make a mask for the upper triangle (triumask)
#             triumask <- upper.tri(matrix(1, nrow = nrow(square_d), ncol = ncol(square_d)), diag = TRUE)
#             # set the upper triangle to NA
#             square_d[triumask] <- NA
#             # transform matrix back to vector
#             triangle_d <- c(square_d)
#         } else if ((((-1 + sqrt(1 + 8 * length(effect_map[[1]]$d))) / 2) + 1) %% 1 == 0) {
#             # else if the map is half a matrix (triangle): 
#             triangle_d <- c(effect_map[[1]]$d)
#         }

#         # for orig_stat
#         # if the map is a full matrix (square):
#         if (sqrt(length(effect_map[[1]]$orig_stat)) %% 1 == 0) {
#             # do nothing
#             square_orig_stat <- matrix(effect_map[[1]]$orig_stat, nrow = sqrt(length(effect_map[[1]]$orig_stat)))
#             # make a mask for the upper triangle (triumask)
#             triumask <- upper.tri(matrix(1, nrow = nrow(square_orig_stat), ncol = ncol(square_orig_stat)))
#             # fill the upper triangle with the map data
#             upper <- square_orig_stat
#             upper[triumask] <- effect_map[[1]]$orig_stat
#             # transform matrix back to vector
#             triangle_orig_stat <- c(upper[upper.tri(upper)])
#         } else if ((((-1 + sqrt(1 + 8 * length(effect_map[[1]]$orig_stat))) / 2) + 1) %% 1 == 0) {
#             # else if the map is half a matrix (triangle):




########### HELPERS FOR CHANGING A SQUARE MATRIX TO A TRIANGLE MATRIX ###########

# FIRST NEED TO REORDER ACCORDING TO MAPPING IF APPLICABLE, THEN TAKE LOWER TRIANGLE, THEN PLOT

# 1. Reorder by mapping if applicable
mapping <- function(effect_map, map_path) {
    # turn effect map into a square matrix
    mat <- matrix(data = effect_map[[1]]$orig_stat, nrow = sqrt(length(effect_map[[1]]$orig_stat)))

    # load map
    mapping <- read.csv(map_path, header = TRUE)

    # reorder based on mapping
    ordered_mat <- mat[mapping$oldroi, mapping$oldroi]

    # return the ordered matrix as a matrix
    return(ordered_mat)
}

# 2. Take lower triangle
lower_triangle <- function(ordered_mat) {
    # if result_as_vector is TRUE, then return the flattened matrix as a vector
    # if include_NA is TRUE, then return the vector with NAs included if results_as_vector is TRUE
    # if include_NA is FALSE, then return the vector without NAs 
    # cannot return matrix if includ_NA is FALSE

    # make a mask for the upper triangle (triumask)
    trilmask <- lower.tri(matrix(1, nrow = nrow(ordered_mat), ncol = ncol(ordered_mat)), diag = TRUE)
    # set the upper triangle to NA
    ordered_mat[trilmask] <- NA
    
    # if (result_as_vector == TRUE) {
    #     if (include_NA == TRUE) {
    #         result = c(ordered_mat)
    #     } else {
    #         result = c(ordered_mat[!is.na(ordered_mat)])
    #     }
    # } else {
    #     if (include_NA == FALSE) {
    #         stop("Cannot return a matrix with include_NA == FALSE")
    #     }
    #     result = ordered_mat
    # }
    return(ordered_mat)
}

# 3. Plot the lower triangle 
plot_lower_triangle <- function(triangle_ordered, mapping_path, show_plot = TRUE) {
    # load mapping
    mapping <- read.csv(mapping_path, header = TRUE)

    # triangle_ordered to matrix
    triangle_ordered <- matrix(triangle_ordered, nrow = sqrt(length(triangle_ordered)))

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
        
    if (show_plot == TRUE) {
        print(heatmap_plot)
    }
}

# 4. Combine all functions
square_to_triangle <- function(effect_map, map_path, show_plot = TRUE) {
    ordered_mat <- mapping(effect_map, map_path)
    triangle_ordered <- lower_triangle(ordered_mat)
    plot_lower_triangle(triangle_ordered, map_path, show_plot)
    return(triangle_ordered)
}

#######e From a triangle, plot a full square matrix

plot_full_mat <- function(triangle_ordered, mapping_path) {
    # load mapping
    mapping <- read.csv(mapping_path, header = TRUE)

    # mirror the triangle across the x = y line to get full matrix
    # first replace the NAs in the lower triangle with 0
    triangle_ordered[is.na(triangle_ordered)] <- 0
    full_mat <- triangle_ordered + t(triangle_ordered) #- diag(diag(triangle_ordered))

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
        
    print(heatmap_plot)
    }
