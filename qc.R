# Run QC on effect size FC matrices and output the plots to QC folder

fc_qc <- function(cleaned_data) {
    for (s in 1:dim(cleaned_data$study)[1]) {
        name <- cleaned_data$study$name[s]
        map_idx <- which(toupper(names(cleaned_data$effect_map)) == name)
        if ((cleaned_data$study$map_type[s] == "FC") & (cleaned_data$study$ref[s] == "Shen_268")) {
            plot_full_mat(cleaned_data$effect_map[[map_idx]]$orig_stat, "/work/neuroprism/effect_size/data/helper_data/map268_subnetwork.csv", export = TRUE, export_path = paste0("/work/neuroprism/effect_size/data/QC/", cleaned_data$study$name[s], "_orig_stat_mat.png"), show_plot = FALSE)
        }
        else if ((cleaned_data$study$map_type[s] == "FC") & (cleaned_data$study$ref[s] == "UKB_55")) {
            plot_full_mat(cleaned_data$effect_map[[map_idx]]$orig_stat, mapping_path = NA, export = TRUE, export_path = paste0("/work/neuroprism/effect_size/data/QC/", cleaned_data$study$name[s], "_orig_stat_mat.png"), show_plot = FALSE)
        }
    }
}