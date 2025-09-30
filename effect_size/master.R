# Master script

# step 1: clean and combine data
# step 2: QC
# step 3: calculate effect size
# step 4: calculate sim CIs

# Input:
#   file names: (date)_(dataset)_(map_type)_(test_type)_(condition)_(outcome).mat
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
#                 - stat
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
#   file_name: combined_data_(date).RData
#
#   - data
#       - <name>
#           - <pooling>_<motion>
#		        			 - stat
#                  - p
#                  - std_brain 
#                  - std_score 
#                  - n         ------------ NaN if two-sample
#                  - n1        ------------ NaN if one-sample
#                  - n2        ------------ NaN if one-sample
#                  - pooling_method (e.g., 'net')
#                  - motion_method (e.g., 'regression')
#                  - d
#                  - ci_lb
#                  - ci_ub
#           - <pooling>_<motion>
#                 ...
#		    - <name>
#			      - ...
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
#       - category
# 
#   - brain_masks
#       - <name>
#           - mask
#       - <name>
#           ...


master <- function(data_dir, script_dir, intermediate_dir, output_dir, num_sdx_r2d = 2, alpha = 0.05, final_output_file = 'braineffex_data') {

  library(oro.nifti)

  # Source everything else in the present directory
  script_filenames <- list.files(script_dir, pattern = "\\.R$", full.names = TRUE)
  script_filenames <- script_filenames[basename(script_filenames) != 'master.R'] # skip master.R
  for (file in script_filenames) {
    source(file)
  }

  # set final output path
  final_output_path = file.path(output_dir, paste0(final_output_file, '_', Sys.Date(), '.RData'))
  
  # load individual group-level datasets, combine, and clean
  cleaned_data <- clean_data(data_dir, script_dir, 'clean_data', intermediate_dir, testing=TRUE)

  # extract data from cleaned_data
  study <- cleaned_data$study
  brain_masks <- cleaned_data$brain_masks
  data <- cleaned_data$data

  # check orientation of FC data with masks
  oriented_list <- check_orientation(data, brain_masks)
  brain_masks <- oriented_list$brain_masks
  data <- oriented_list$data

  # calculate Cohen's d and simultaneous confidence intervals for each study
  d_maps <- calc_d(study, data, output_dir = intermediate_dir)

  # checker to check dimensions (and probably more things eventually)
  data <- checker(d_maps, int_dir = intermediate_dir)

  # load template data and phen_keys
  template <- readNIfTI("data/template_nifti.nii.gz") # assumes MNI - TODO: get actual ref
  template <- template@.Data

  # combine all data together into variable v for saving purposes
  v <- list(study = study, brain_masks = brain_masks, data = data, template = template)

  # run meta-analysis and add results to v
  v <- meta_analysis(v, v$brain_masks, grouping_var = "category")

  # save the final results
  save(v, file = final_output_path)
  cat(paste('Results saved in',final_output_path))
}

