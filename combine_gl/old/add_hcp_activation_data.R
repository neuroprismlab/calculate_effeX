## add HCP task-activation data to d_clean
## needs to be run locally for now... should be run on Discovery in the future

# load d_clean
load("~/Google Drive/My Drive/NeuroPRISM/effect_size_shiny_git/effect_size_shiny/data/effect_maps_clean.Rdata")
# loads effect_map and study

# TODO put this into a for loop:
for (task in c("emotion", "gambling", "relational", "social", "wm")) {
    task_l = task
    task_u = toupper(task_l) # upper-case

    n = case_when(
        task_l == "emotion" ~ 1022, # TODO: check these numbers since they differ from the n's in the hcp helper files...
        task_l == "gambling" ~ 1057,
        task_l == "relational" ~ 1016,
        task_l == "social" ~ 1027,
        task_l == "wm" ~ 1058
    )

    for (i in 1:length(task_l)) {

        # load nifti
        data_files <- list.files(path = "/Users/neuroprism/Library/CloudStorage/GoogleDrive-halleeninet@gmail.com/.shortcut-targets-by-id/17uYR-Ubbo9n0459awrNyRct0CL4MwAXl/Hallee-Steph share/visualize_effects_app/data", full.names = TRUE)
        pattern <- paste0(task_u[i], ".*\\.nii\\.gz")
        matching_file <- grep(pattern, data_files, value = TRUE)
        nifti <- readnii(matching_file)

        nifti_data <- img_data(nifti)

        list <- as.list(nifti_data)

        # remove zero values (likely those outside of the brain, but could theoretically include some in the brain)
        nonzero_list <- list[list!=0]

        # convert to numeric
        nonzero_numeric <- unlist(nonzero_list)

        # create new list entry to add to cleaned_data$effect_map
        new_data <- list(d = nonzero_numeric, 
            orig_stat = nonzero_numeric,
            n = n)

        new_study_data <- data.frame(
            basefile = basename(matching_file),
            folder = dirname(matching_file),
            name = sub('\\.nii\\.gz$', '', basename(matching_file)),
            ext = ".nii.gz",
            dataset = "hcp",
            map_type = "act",
            orig_stat_type = "d",
            var1 = task_u,
            var2 = "rest")

    # combine new data with cleaned_data
    effect_map[[length(effect_map) + 1]] <- new_data
    effect_map <- `names<-`(effect_map, `[<-`(names(effect_map), (length(effect_map)), new_study_data$name))
    study <- rbind(study, new_study_data)

    }
}


# export d_clean result
save(effect_map, study, file = '/Users/neuroprism/Library/CloudStorage/GoogleDrive-halleeninet@gmail.com/My Drive/NeuroPRISM/effect_size_shiny_git/effect_size_shiny/data/effect_maps_clean_hcp.Rdata')
