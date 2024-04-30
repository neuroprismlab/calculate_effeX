# add phenotypic categories to study dataframe

add_phen <- function(study, phen_file = "phen.csv", data_dir = "data/") {
    # load phenotypic data file (phen_file) from data directory
    phen <- read.csv(paste0(data_dir, phen_file), header = TRUE)
    # merge phenotypic data with study data
    # TODO: check to see if the data merges properly, and if the study names are in the same format in both (e.g. capitalization, _ vs. .)
    phen_study <- merge(study, phen, by = "name")

    return(phen_study)
}