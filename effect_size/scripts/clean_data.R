# Combine and clean all group level studies
# cleaned_data <- clean_data(data_dir, script_dir, 'clean_data', intermediate_dir, testing=TRUE)

clean_data <- function(data_dir = data_dir,
                       script_dir = script_dir,
                       output_file = 'clean_data',
                       intermediate_dir = intermediate_dir,
                       testing = FALSE) {
  
  output_path = file.path(intermediate_dir, paste0(output_file, '_', Sys.Date(), '.RData'))
  
  print(data_dir)
  # Get list of all files
  all_items <- list.files(data_dir, pattern = "\\.mat$", full.names = TRUE)
  
  # Filter out directories
  mat_files <- all_items[file.info(all_items)$isdir == FALSE]
  
  # Get only the file names
  mat_file_names <- basename(mat_files)
  
  # if (testing) {
  #   #mat_file_names <- mat_file_names[grep("^04-Sep-2024hcp_fc_t2_REST_Gender", mat_file_names)]
  #   mat_file_names <- mat_file_names[grepl("Sep", mat_file_names)]
  #   #mat_file_names <- mat_file_names[-grep("^04-Sep-2024hcp_fc_t2_REST_Gender", mat_file_names)]
  # }
  
  # create empty dataframe for study info
  column_names <- c("basefile", "folder", "name", "ext", "dataset", "map_type",
                    "orig_stat_type", "test_component_1", "test_component_2", 
                    "category", "ref") 
  # TODO: add level_map
  study <- setNames(data.frame(matrix(ncol = length(column_names), nrow = 0)), column_names)
  
  # create empty list for stat_maps (list of lists)
  stat_maps <- list()
  
  # create empty list for brain masks
  brain_masks <- list()
  
  # loop through all file names
  for (file in mat_file_names) {
    print(paste("loading", file))
    # load the file
    mat_struct = readMat(file.path(data_dir, file))
    if (grepl('_act_', file)) { # special correction for readMat bug (nested lists)
      mat_struct <- readMat_act_correction(mat_struct)
    }
    
    data <- assign_names(mat_struct)
    
    # extract study info
    study_info <- data$results$study.info
    
    # fill in study info for this file
    # only add test_component_2 if it exists
    
    test_component_2 <- NA
    
    # Check if there are more than one test components and assign accordingly
    if (length(study_info$test.components) > 1) {
      test_component_2 <- study_info$test.components[[2]][[1]]
    }
    
    # TODO: make sure test component 2 is the task for activation studies, currently component 2 is rest
    
    name = sub("\\.mat$", "", file)
    
    ref = ifelse(study_info$map == "act", "voxel", ifelse(study_info$dataset == "ukb", "ukb_55", ifelse(grepl("_rbc", study_info$dataset), "schaefer_200", "shen_268")))
    
    new_row <- data.frame(basefile = file, folder = data_dir, name = name, ext = ".mat",
                          dataset = study_info$dataset, map_type = study_info$map,
                          orig_stat_type = study_info$test,
                          test_component_1 = study_info$test.components[[1]][[1]],
                          test_component_2 = test_component_2,
                          category = study_info$category,
                          ref = ref
    )
    
    # add to study dataframe

    study <- rbind(study, new_row)
    
    
    ## Now create stat_maps as list of lists
    # each study has one list, and within each study there is a list for each 
    # combination of motion and pooling
    stat_maps[[name]] <- data$results$data
    
    # add brain mask to brain_masks
    brain_masks[[name]]$mask <- data$results$study.info$mask
    
  }
  
  names(brain_masks) <- tolower(names(brain_masks))
  names(stat_maps)   <- tolower(names(stat_maps))
  
  if (testing) {
    save(study, stat_maps, brain_masks, file = output_path)
  }
  return(list(study = study, data = stat_maps, brain_masks = brain_masks))
}
