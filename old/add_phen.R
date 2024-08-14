# add phenotypic categories to study dataframe

add_phen <- function(study, effect_maps, phen_file = "phen.csv", data_dir = "data/", output_file = "phen_study.RData") {
    # load phenotypic data file (phen_file) from data directory
    phen <- read.csv(paste0(data_dir, phen_file), header = TRUE)
    # merge phenotypic data with study data
    # TODO: check to see if the data merges properly, and if the study names are in the same format in both (e.g. capitalization, _ vs. .)
    phen_study <- merge(study, phen, by = "name")

    for (i in 1:dim(phen_study)[1]) {
        name <- phen_study[i, "name"]
        len <- length(effect_maps[[which(toupper(names(effect_maps)) == name)]]$d)
        if (phen_study[i,"map_type"] == "FC") {
            if (len == 35778) {
                stop(paste0("Triangle matrix found, please ensure all matrices are converted to full squares with the triangle_to_square helper function. Study name: ", name))
            }
            else if (sqrt(len) == 55) {
            ref <- "UKB_55"
            }
            else if (sqrt(len) == 268) {
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

    save(phen_study, file = output_file)

    return(phen_study)
}
