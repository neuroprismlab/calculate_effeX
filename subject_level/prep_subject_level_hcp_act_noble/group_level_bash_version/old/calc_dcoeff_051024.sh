#!/bin/bash
####################################
#
# calculate t-stat / d for HCP task-based activation
#
# requirements: FSL
#
# prereq:
#   1. copy task files, masks, and motion data locally: get_hcp_data_for_dcoeff.sh
#   2. combine RL and LR motion: average_RL_LR_motion.sh
#
####################################


# Set paths and params

export PATH=$(dirname "$0"):$PATH
source setparams_activation.sh


# Setup: make temp files

temp_log="${tmp_dir}/log"
group_mask="${tmp_dir}/group_mask.nii.gz"
stats_4D_pos="${tmp_dir}/all_sub_stats.nii.gz"
stats_4D_neg="${tmp_dir}/all_sub_stats__neg.nii.gz"
masks_4D="${tmp_dir}/all_sub_masks.nii.gz"

mkdir -p $tmp_dir
mkdir -p $results_dir

# Calculate

for task_contrast_pair in "${task_contrast_pairs[@]}"; do
    
    read -r task cope <<< "$task_contrast_pair"

    echo "Getting filenames"

    # Define group results prefix
    group_stat_prefix_pos="${results_dir}/${task}_cope${cope}_pos"
    
    # Get all statistic files, then find/replace to create mask and motion filenames
    stat_filenames=($stat_dir/${task}_$cope/*/MNINonLinear/Results/tfMRI_${task}/*/cope${cope}.feat/stats/tstat1.nii.gz)
    mask_filenames=("${stat_filenames[@]/tstat1.nii.gz/mask.nii.gz}")
    
    # motion, which has a slightly different path
    base_dirs_motion=("${stat_filenames[@]%%/MNINonLinear/*}")
    base_dirs_motion=("${base_dirs_motion[@]/activation\/${task}_$cope/motion_average}")
    motion_filenames=()
    for base_dir in "${base_dirs_motion[@]}"; do
        motion_filenames+=("$base_dir/MNINonLinear/Results/tfMRI_${task}/Movement_RelativeRMS_mean.txt")
    done

    #echo ${stat_filenames[2]}
    #echo ${motion_filenames[2]}

    echo "Calculating stats"

    # Create group mask

    [[ -f ${group_mask} ]] && rm ${group_mask} # remove any previous analysis
    fslmerge -t ${masks_4D} ${mask_filenames[@]}
    fslmaths ${masks_4D} -Tmin ${group_mask}
    # OLD: #3dMean -prefix ${group_mask} -mask_inter ${mask_filenames[@]} &> ${temp_log}

    # Do t-test
    
    [[ -f ${stats_4D_pos} ]] && rm ${stats_4D_pos} # remove any previous analysis
    [[ -f ${group_stat_prefix_pos}_tstat1.nii.gz ]] && rm ${group_stat_prefix_pos}_tstat1.nii.gz # remove any previous analysis
    fslmerge -t ${stats_4D_pos} ${stat_filenames[@]}
    randomise -i ${stats_4D_pos} -m ${group_mask} -1 -o ${group_stat_prefix_pos} >> ${temp_log}
    
#    # For reference: to test in negative direction - TODO: remove "_pos" if not using this
#    #group_stat_prefix_neg="${results_dir}/${task}_cope${cope}_neg"
#    #fslmaths ${stats_4D_pos} -mul -1 ${stats_4D_neg}
#    #randomise -i ${stats_4D_neg} -m ${group_mask} -1 -o ${group_stat_prefix_neg} >> ${temp_log}

    # Calculate ground truth Cohen's d
    [[ -f ${group_stat_prefix_pos}_dcoeff.nii.gz ]] && rm ${group_stat_prefix_pos}_dcoeff.nii.gz # remove any previous analysis
    sqrt_n=$(echo "sqrt(${#stat_filenames[@]}-1)" | bc)
    fslmaths ${group_stat_prefix_pos}_tstat1.nii.gz -div $sqrt_n ${group_stat_prefix_pos}_dcoeff.nii.gz
    # OLD: 3dcalc -a ${group_stat_prefix_pos}_tstat1.nii.gz -expr 'a/'"$sqrt_n" -prefix ${group_stat_prefix_pos}_dcoeff.nii.gz

   
done




