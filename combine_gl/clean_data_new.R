# step 1: clean and combine data
# step 2: QC
# step 3: calculate effect size
# step 4: calculate sim CIs

# STEP 1
# Input:
#   file names: (dataset)_(map_type)_(contributor)_(date)_(extra_description).mat
#   - study_info
#       - dataset
#       - <test_components> (e.g., {'condition_label', 'score_label'}
#       - map
#       - test
#       - level_map (optional)
#       - brain_mask
#       - category
#   - data
#       - <pooling strategy>
#            - <motion strategy>
#                 - b_standardized
#                 - p
#                 - std_brain
#                 - std_score
#                 - n         ------------ NaN if two-sample
#                 - n1        ------------ NaN if one-sample
#                 - n2        ------------ NaN if one-sample
#                 - pooling_method (e.g., 'net')
#                 - motion_method (e.g., 'regression')
#
# Output: 
#   file_name: combined_maps_(date).RData
#
#   - data
#       - <name>
#           - <pooling>
#               - <motion>
#		        			 - b_standardized
#                  - p
#                  - std_brain 
#                  - std_score 
#                  - n         ------------ NaN if two-sample
#                  - n1        ------------ NaN if one-sample
#                  - n2        ------------ NaN if one-sample
#                  - pooling_method (e.g., 'net')
#                  - motion_method (e.g., 'regression')
#		    - <name>
#			    - ...
#
#   - study_info (table)
#       - basefile
#       - folder
#       - name
#       - ext
#       - dataset
#       - map_type
#       - orig_stat_type
#       - test_component_1
#       - test_component_2

clean_data <- function(data_dir = '/work/neuroprism/effect_size/data/group_level/',
                       script_dir = '/home/h.shearer/hallee/calculate_effeX/combine_gl',
                       output_file = paste0('/work/neuroprism/effect_size/data/combined_gl/output/combined_maps_', Sys.Date(), '.RData'),
                       testing = FALSE) {
  
  source(file.path(script_dir, 'helpers.R'))
  
  # Get list of all files
  # Get list of all .mat files and directories
  all_items <- list.files(data_dir, pattern = "\\.mat$", full.names = TRUE)
  
  # Filter out directories
  mat_files <- all_items[file.info(all_items)$isdir == FALSE]
  
  # Get only the file names
  mat_file_names <- basename(mat_files)
  
  if (testing) {
    mat_file_names <- mat_file_names[grep("^03-Sep", mat_file_names)]
  }
  
  # create empty dataframe for study info
  column_names <- c("basefile", "folder", "name", "ext", "dataset", "map_type",
                    "orig_stat_type", "test_component_1", "test_component_2") 
  # TODO: add level_map
  study <- setNames(data.frame(matrix(ncol = length(column_names), nrow = 0)), column_names)
  
  # create empty list for stat_maps (list of lists)
  stat_maps <- list()
  
  # loop through all file names
  for (file in mat_file_names) {
    print(c("loading file ", file))
    # load the file
    mat_struct = readMat(file.path(data_dir, file))
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
    
    name = sub("\\.mat$", "", file)
    
    new_row <- data.frame(basefile = file, folder = data_dir, name = name, ext = ".mat",
                          dataset = study_info$dataset, map_type = study_info$map,
                          orig_stat_type = study_info$test,
                          test_component_1 = study_info$test.components[[1]][[1]],
                          test_component_2 = test_component_2
    )
    
    # add to study dataframe
    study <- rbind(study, new_row)
    
    
    ## Now create stat_maps as list of lists
    # each study has one list, and within each study there is a list for each 
    # combination of motion and pooling
    
    stat_maps[[name]] <- data$results$data
    
    # add brain mask to stat_maps as well
    stat_maps[[name]]$brain_mask <- data$results$study.info$mask
    
  }
  save(study, stat_maps, file = output_file)
  return(list(study = study, data = stat_maps))
}