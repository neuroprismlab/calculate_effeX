# add phenotypic categories to study dataframe

add_phen <- function(study, phen_file = "phen.csv", data_dir = "data/", output_file = "phen_study.Rdata") {
    # load phenotypic data file (phen_file) from data directory
    phen <- read.csv(paste0(data_dir, phen_file), header = TRUE)
    # merge phenotypic data with study data
    # TODO: check to see if the data merges properly, and if the study names are in the same format in both (e.g. capitalization, _ vs. .)
    phen_study <- merge(study, phen, by = "name")

    for (i in 1:dim(phen_study)[1]) {
        name <- phen_study[i, "name"]
        len <- length(effect_map[[which(toupper(names(effect_map)) == name)]]$orig_stat)
        if (phen_study[i,"map_type"] == "FC") {
            if (len == 35778) {
            parc <- "Shen_268_node"
            }
            else if (len == 1485) {
            parc <- "UKB_55_node"
            }
            else if (len == 71824) {
            parc <- "Shen_268_node_full"
            }
            else {
            stop(paste0("Unknown parcellation found, please add this parcellation. Length of effect map: ", len, ". study name: ", name))
            }
        }
        
        else if (phen_study[i, "map_type"] == "ACT") {
            parc <- NA
        }
        # add parc column
        phen_study$parc[i] <- parc
    }

    save(phen_study, file = output_file)

    return(phen_study)
}