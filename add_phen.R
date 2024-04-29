# add phenotypic categories to study dataframe

add_phen <- function(study, phen_file = "phen.csv") {
    # load phenotypic data file (phen_file)
    phen <- read.csv(phen_file, header = TRUE)
    # merge phenotypic data with study data
    # TODO: check to see if the data merges properly, and if the study names are in the same format in both (e.g. capitalization, _ vs. .)
    phen_study <- merge(study, phen, by = "study")