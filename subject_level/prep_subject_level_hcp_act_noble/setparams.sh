#######################################################  
#
# Set paths and params
#
# hcp data example path: s3://hcp-openaccess/HCP_1200/<subID>/MNINonLinear/Results/tfMRI_<task>/tfMRI_<task>_hp200_s<?>_level2vol.feat/cope<cope>.feat/stats/tstat1.nii.gz
#
#######################################################  

# Params

test=0
get_motion=0
get_stats=1

# task-contrast task_contrast_pairs

task_contrast_pairs=("SOCIAL 6") 
#task_contrast_pairs=("SOCIAL 6" "WM 20" "RELATIONAL 4" "GAMBLING 6" "EMOTION 3")
motion_method_params=("none" "regression" "threshold")

# Directories

# remote (AWS)
source_dir='s3://hcp-openaccess/HCP_1200' 

# local data dir
target_dir='/work/neuroprism/data_shared/hcp/'
stat_dir="$target_dir/activation"
motion_dir="$target_dir/motion"
mean_motion_dir="${motion_dir}_average"

# result effect map dir
#results_dir='/work/neuroprism/effect_size/data/group_level/intermediates/hcp/'
#results_dir='/work/neuroprism/effect_size/data/group_level_intermediates/hcp_act/'


# temp dir
tmp_dir="$target_dir/temp"



