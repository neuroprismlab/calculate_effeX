##############################################
#
# Set params and simple strings
#
##############################################

###### General ######

testing <- 1

# What to skip (for summary only?)
skip_nii <- 1  # may want to skip if connecting to remote data source
skip_hbn <- 0  # problematic - too few subjects - may want to skip
skip_nii__vis <- 0
filter__vis <- 1 # problematic - too few subjects - may want to skip
make_figs<-0
save_figs<-1
do_log<-0
do_meta_all_edges<-0


###### Loading and cleaning ######

# required fields to load for each statistic type
req_fields <- list(
    d = c("d", "n"), # only pre-calculated for t activation map niftis
    t = c("stats", "n"),
    t2 = c("stats", "n1", "n2"),
    r = c("r", "n")
)

###### Statistics ######

# effect size conversion params
num_sdx_r2d <- 1  ## number of standard deviations in X to use for Maya's r-to-d conversion

# stats: CI params
# note that we're using 0.025 for two-tailed; Matlab default for t-test is p=0.05 
this_alpha<-0.025 # statistical threshold, not plotting transparency

# effect size classification and characteristics
conventional_d_thresh <- c(0.2,0.5,0.8) # Cohen's d thresholds for conventional effect size classification
ttest_types <- c("one.sample","two.sample")
desired_stat_order <- c("r","t2","t","t (act)")


###### Plotting ######

# plot params
ylim_ci_plot <- c(-1.2,1.2) #c(-0.8, 0.8)
ylim_ci_plot_smaller <- c(-0.4,0.4)
breaks_ci_plot_smaller <- c(-0.3, 0, 0.3)
multiplot_width <- 2.8
multiplot_height <- 57

forest_img_height <- 15
forest_img_width <- 12
forest_img_pointsize <- 18
box_border_thickness<- 1.5 # TODO: double check
legend_closeness <- 0.8
# unique_cols_by_dataset<-c("#003f5c","#7a5195","#ff764a","#bc5090","#ef5675","#374c80","#ffa600")


###### Directory and Filenames ######

# directories and file strings
data_dir <- '/Users/stephanienoble/Documents/data/mnt'
atlas_tools_dir <- '/Users/stephanienoble/Library/CloudStorage/GoogleDrive-smn33@yale.edu/My Drive/Lab/More/Software/scripts/R/myscripts/effect_size/atlas_tools'
# atlas_tools_dir <- '/Users/stephanienoble/Library/CloudStorage/GoogleDrive-stephanie.noble@yale.edu/My Drive/Lab/More/Software/scripts/R/myscripts/effect_size/atlas_tools'

exts <- c('mat$', 'nii$', 'nii.gz$')

if (Sys.getenv("USER") == "steph" && Sys.info()["nodename"]=="StephsMBP") { # old laptop
    local_project_dir <- '/Users/steph/Library/CloudStorage/GoogleDrive-smn33@yale.edu/My Drive/Lab/More/Tasks-Ongoing/K99/Effect_Size/'
    # local_project_dir <- '/Users/steph/Library/CloudStorage/GoogleDrive-stephanie.noble@yale.edu/My Drive/Lab/More/Tasks-Ongoing/K99/Effect_Size/'
} else {
    local_project_dir <- '/Users/stephanienoble/Library/CloudStorage/GoogleDrive-smn33@yale.edu/My Drive/Lab/More/Tasks-Ongoing/K99/Effect_Size/'
    # local_project_dir <- '/Users/stephanienoble/Library/CloudStorage/GoogleDrive-stephanie.noble@yale.edu/My Drive/Lab/More/Tasks-Ongoing/K99/Effect_Size/'
}

local_effect_map_dir <- paste0(local_project_dir, 'data/')
results_dir <- paste0(local_project_dir, 'results/')


# set skip strings for loading
if (skip_nii) {
    skip_act_str <- '_no_act'
} else {
    skip_act_str <- ''
}

if (skip_hbn) {
    datasets_to_skip <- c('hbn')
    skip_datasets_str <- paste0('_no_', datasets_to_skip)
} else {
    datasets_to_skip <- NULL
    skip_datasets_str <- ''
}
if (testing) {
    testing_str <- '_testing'
} else {
    testing_str <- ''
}
skip_str <- paste0(skip_act_str, skip_datasets_str,testing_str)

# set skip strings for visualization
if (skip_nii__vis) {
    skip_act_str__vis <- '_no_act'
} else {
    skip_act_str__vis <- ''
}

if (filter__vis) {
    skip_studies_str__vis <- '_filtered' # TODO: move studies to skip here
} else {
    skip_studies_str__vis <- ''
}
skip_str__vis <- paste0(skip_act_str__vis, skip_studies_str__vis)


## Output filenames ##

# original loaded data

full_dataset_filename <- paste0(local_effect_map_dir, 'full_data', skip_str) # TODO: where is this used? is this supposed to be for the precursor file?

effect_maps_precursor_filename <- paste0(local_effect_map_dir, 'effect_maps_precursor', skip_str, '.RData')

# effect maps

effect_maps_filename <- paste0(local_effect_map_dir, 'effect_maps', skip_str, '.RData')

## Plot filenames ##

# hist
combined_filename <- paste0(results_dir, 'esz_hist', skip_str)
individual_tasks_filename <- paste0(results_dir, 'eszindvid_hist', skip_str)
hist_by_study_filename <- paste0(results_dir, 'hist_by_study', skip_str__vis, '.png')
hist_by_stat_filename <- paste0(results_dir, 'hist_by_stat', skip_str__vis, '.png')
# CI ribbon 
d_ci_ribbon_filename <- paste0(results_dir, 'd_ci_ribbon_by_study', skip_str__vis, '.png')
d_ci_ribbon_small_filename <- paste0(results_dir, 'd_ci_ribbon_by_study__small', skip_str__vis, '.png')
d_conservative_filename <- paste0(results_dir, 'd_conservative_by_study', skip_str__vis, '.png')
d_conservative_small_filename <- paste0(results_dir, 'd_conservative_by_study__small', skip_str__vis, '.png')
# matrix
matrices_by_study_filename <- paste0(results_dir, 'matrices_by_study', skip_str__vis, '.pdf')
matrices_by_cat_filename <- paste0(results_dir, 'matrices_by_cat', skip_str__vis, '.png')
# meta-analysis
forest__by_stat_filename <- paste0(results_dir, 'forest__by_stat', skip_str__vis, ".pdf") # TODO: put in setparams
funnel__by_stat_filename <- paste0(results_dir, 'funnel__by_stat', skip_str__vis, ".pdf") # TODO: put in setparams

## Log filenames ##

# general
log_dir <- paste0(results_dir, 'logs/')
study_dataset_info_filename=paste0(log_dir,"study_dataset_info.png")
n_for_conventional_dthresh_filename=paste0(log_dir,"n_for_conventional_dthresh.png")
# perc ci_contains_zero
conservative_effect_summary_by_statistic_filename<-paste0(log_dir,"conservative_effect_summary_by_statistic",skip_str__vis,".txt")
# perc_zero_by_stat_filename=paste0(log_dir,"perc_zero_by_stat",skip_str__vis,".txt")
# max_conservative_effect_by_stat_filename=paste0(log_dir,"max_conservative_effect_by_stat",skip_str__vis,".txt")
# perc d above threshold
d_perc_gt_thresh1_filename=paste0(log_dir,"perc_gt_dthresh_",conventional_d_thresh[1],".png")
d_perc_gt_thresh2_filename=paste0(log_dir,"perc_gt_dthresh_",conventional_d_thresh[2],".png")
# typical esz from meta 
mdl__by_stat_filename <- paste0(log_dir, 'mdl__by_stat', skip_str__vis, ".txt")
mdl__by_stat__within_pairwise_filename <- paste0(log_dir, 'mdl_by_stat__within_pairwise', skip_str__vis, ".txt")
mdl__by_stat__between_pairwise_filename <- paste0(log_dir, 'mdl_by_stat__between_pairwise', skip_str__vis, ".txt")

