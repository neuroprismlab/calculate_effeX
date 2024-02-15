#######################
# Plot matrices
#######################

library(lattice)    
library(gridExtra)
library(scales) # for oob, to squish values above and below the color limits to the limit so they're not plotted as grey
source("/Users/stephanienoble/Library/CloudStorage/GoogleDrive-stephanie.noble@yale.edu/My Drive/Lab/More/Software/scripts/R/myscripts/effect_size/atlas_tools/setpaths_atlas_tools")

# Set plot params
n_plot_cols__by_study <- 4
n_plot_cols__by_cat <- 4
pdf_height <- 20
pdf_width <- 13
png_width <- 1600
png_height <- 600
panel_title_fontsize <- 0.8
colaxis_d_thresh <- 0.2 # 1.2
# colaxis <- c(-colaxis_d_thresh,colaxis_d_thresh)
col_palette_bands <- 300
color.palette <- colorRampPalette(c("blue", "white", "red"))(col_palette_bands)


# Dcoeff plots

mat_to_plot__d<-list()
for (this_study in study_info$study) {
    if (study_info[this_study,"n_edges"]==35778) { # only for Shen atlas

        this_d<-d_master__reordered[[this_study]]
        this_n_nodes<-study_info[this_study,"n_nodes"]

        triumask <- upper.tri(matrix(1, this_n_nodes, this_n_nodes))

        mat_to_plot__d[[this_study]]<-triumask
        mat_to_plot__d[[this_study]][triumask] <- this_d

        # flip so we can see them all
        mat_to_plot__d[[this_study]] <- mat_to_plot__d[[this_study]] + t( mat_to_plot__d[[this_study]])

        # TODO: under construction. DON'T UNCOMMENT! Now doing this before this script
        # if (!grepl("hcp", this_study)) {
        #     # reorder_matrix_by_atlas
        #     mat_to_plot__d[[this_study]] <- reorder_matrix_by_atlas(mat_to_plot__d[[this_study]])
        # }

    } else {
        # otherwise skip non-Shen plots and do activation plots in Bisweb
        mat_to_plot__d[[this_study]]<-triumask
        mat_to_plot__d[[this_study]][triumask] <- 0
    }
}

# mat_to_plot__d2<-mat_to_plot__d[1:8]
# plots <- lapply(names(mat_to_plot__d2), 
#                 function(x)
#                     levelplot(mat_to_plot__d2[[x]], 
#                     col.regions = color.palette,  
#                     at = seq(-colaxis_d_thresh, colaxis_d_thresh, 
#                     length.out = col_palette_bands), 
#                     scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), 
#                     main=list(label=x,cex=panel_title_fontsize)) )

# for each study, by stat
for (this_stat in desired_stat_order) {
    this_matrices_by_study_filename <- paste(results_dir, "/matrices_by_study_", this_stat, ".pdf", sep="")
    pdf(file=this_matrices_by_study_filename, height = pdf_height, width = pdf_width)  # Adjust the dimensions as per your requirement
    # which names contain this_stat
    this_idx <- grep(paste("_",this_stat,"_",sep=""), names(mat_to_plot__d))

    plots <- lapply(names(mat_to_plot__d[this_idx]),
        function(x) draw_atlas_boundaries(mat_to_plot__d[[x]]) +
            ggtitle(x) +
            theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,-1,0))) +
            theme(legend.position = "none") ) 
    # add legend to last plot
    # plots[[length(plots)]] <- plots[[length(plots)]] + 
    #     theme(legend.position = "right", 
    #         plot.margin = margin(0, 0, 0, 0, "cm")) +
    #     guides(fill = guide_colorbar(title = expression(italic("d")), title.hjust = 0.5))
    
    grid.arrange(grobs = plots, ncol = ceiling(length(plots)/5))
    # grid.arrange(grobs = plots, ncol = n_plot_cols__by_study)
    dev.off()
}


# Meta-analysis dcoeff plots

mat_to_plot<-vector(mode = "list", length=dim(cat_beta_matrix)[2])
names(mat_to_plot)<-colnames(cat_beta_matrix)
for (this_idx in 1:dim(cat_beta_matrix)[2]) {
    this_beta<-cat_beta_matrix[,this_idx]
    this_n_nodes <- ((-1+sqrt(1+8*length(this_beta)))/2)+1
    triumask <- upper.tri(matrix(1, this_n_nodes, this_n_nodes))
    mat_to_plot[[this_idx]]<-triumask
    mat_to_plot[[this_idx]][triumask] <- this_beta
}

png(file=matrices_by_cat_filename, width=png_width, height=png_height)  # Adjust the dimensions as per your requirement
plots <- lapply(names(mat_to_plot),
    function(x) draw_atlas_boundaries(mat_to_plot[[x]]) +
        ggtitle(gsub("nonbrain_code", "", x)) +
        theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,-1,0))) +
        theme(legend.position = "none") )
# add legend to last plot
plots[[length(plots)]] <- plots[[length(plots)]] + 
    theme(legend.position = "right", 
          plot.margin = margin(0, 0, 0, 2, "cm")) +
    guides(fill = guide_colorbar(title = expression(italic("d")), title.hjust = 0.5))
grid.arrange(grobs = plots, ncol = n_plot_cols__by_cat)
dev.off()
# plots <- lapply(names(mat_to_plot),
#     function(x)
#         ggcorrplot(mat_to_plot[[x]],outline.color = "#00000000") +
#         scale_fill_gradient2(limits = c(-colaxis_d_thresh, colaxis_d_thresh),low = "blue", high =  "red", oob = squish) + 
#         ggtitle(x) +
#         theme(panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             plot.title = element_text(hjust = 0.5, margin=margin(0,0,-1,0)))
#         )
# https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots

# # another simple option
# png(file=matrices_by_cat_filename, width=png_width, height=png_height)  # Adjust the dimensions as per your requirement
# plots <- lapply(mat_to_plot, levelplot)
# grid.arrange(grobs = plots, ncol = n_plot_cols__by_cat)
# dev.off()

# # another option, using simpler ggplot (see https://blog.aicry.com/r-heat-maps-with-ggplot2/index.html)
# # this ends up being flipped x-y
# longData<-melt(mat_to_plot[[1]])
# longData<-longData[longData$value!=0,]
# ggplot(longData, aes(x = Var2, y = Var1)) + 
#     geom_raster(aes(fill=value)) + 
#     coord_equal() +
#     scale_fill_gradient2(limits = c(-colaxis_d_thresh, colaxis_d_thresh),low = "blue", high =  "red", oob = squish) +
#     theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# # Plot matrices

# # load d, study
# # make functions: structure_data, draw_atlas_boundaries, summarize_matrix_by_atlas, bipolar, etc.

# # Set params
# n_nodes <- 268

# # For matrices
# study_names <- names(d)
# fullmask <- matrix(1, n_nodes, n_nodes)

# do_plot <- 0
# do_summary <- 1

# if (do_plot) {
#   ## Matrices
#   # Group some study types together
#   # Note: ukb and hcp_act need unique plots
#   study_name_ids_age <- c(1, 13, 16) # age
#   study_name_ids_gender <- c(4, 7, 9, 15, 27) # t2 gender
#   study_name_ids_psych <- c(3, 8, 14, 21) # r cognition / psychological
#   study_name_ids_all <- 1:27

#   study_name_ids <- study_name_ids_all
#   it <- 1
#   par(mfrow = c(2, ceiling(length(study_name_ids) / 2)))

#   for (i in study_name_ids) {
#     study_name <- study_names[i]
#     triside <- "upper"
#     if (study$dataset[i] == "hcp" && study$map_type[i] == "fc" && study$stat_type[i] != "d") {
#       triside <- "lower"
#     }
#     if (length(d[[study_name]]) == n_nodes^2) {
#       m <- structure_data(d[[study_name]], ismatrix = TRUE, mask = fullmask)
#     } else {
#       m <- structure_data(d[[study_name]], ismatrix = TRUE, triangleside = triside)
#     }

#     plot(m, main = study_name)
#     # summarize_matrix_by_atlas(m)
#     it <- it + 1
#   }

#   colormap(bipolar(), 0.1)
#   dev.new(width = 50, height = 10)
#   dev.off()
# }




# Add boundaries
# x_bound <- 5
# y_bound <- 5

# rect(0.5, 0.5, nrow(adj_matrix) + 0.5, ncol(adj_matrix) + 0.5, border = "red", lwd = 2)


## OLD 

# Old setup

    # grid_rows <- 2
    # n_plot_cols__by_study <- length(matrix_plot_study_ids)/grid_rows
    # par(mfrow = c(grid_rows, n_plot_cols__by_study))
    # plot.new()
    # row_index <- ceiling(i / n_plot_cols__by_study)
    # col_index <- i %% n_plot_cols__by_study
    # if (col_index == 0)
    #     col_index <- n_plot_cols__by_study

    # xleft <- (col_index - 1) / n_plot_cols__by_study
    # xright <- col_index / n_plot_cols__by_study
    # ytop <- 1 - (row_index - 1) / grid_rows
    # ybottom <- 1 - row_index / grid_rows
        

# Old plotter in for loop
        # Plot matrix
        # matrix_filename<-paste0(results_dir, this_study, '_matrix.png')
        # png(matrix_filename, width = 500, height = 500, units = "px", res = 300)
        # plt[[i]]<-levelplot(t2, col.regions = color.palette, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)),  at = seq(colaxis[1], colaxis[2], length.out = 100))
        # dev.off()

        # par(mar = c(1,1,1,1))
        # image(t(apply(2*t2,2,rev)),col = colorRampPalette(c("blue", "white", "red"))(256), zlim = c(colmin, colmax), asp=1) # because "image" plot is rotated...
        # legend("bottom", legend = c(colmin, colmax), fill = colorRampPalette(c("blue", "white", "red"))(256), bty = "n", title = "Values")
        # grid.arrange(plt1,plt2, nrow=grid_rows)


        # Add color legend
        # legend(x = 1, y = 5, legend = list(at = seq(colmin, colmax, length.out = 5), col = color.palette(5)), 
        #              bty = "n", title = "Values")
