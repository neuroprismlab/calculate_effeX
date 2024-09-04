#!/bin/bash
#set -x
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

# Setup: define temp filenames and create temp dirs

temp_log="${tmp_dir}/log"
group_mask="${tmp_dir}/group_mask.nii.gz"
stats_4D="${tmp_dir}/all_sub_stats.nii.gz"
masks_4D="${tmp_dir}/all_sub_masks.nii.gz"
design_matrix="${tmp_dir}/design.mat"
design_contrast="${tmp_dir}/design.con"

mkdir -p $tmp_dir
mkdir -p $intermediates_dir

# Remove any previous temp intermediates from tempdir
rm "${tmp_dir}/"*.*
# TODO: confirm that this doesn't remove log
# TODO: check whether intermediates already exist in $group_stat_dir and prompt whether to replace 
# TODO: ls files in this dir that will be removed and ask confirmation first
#[[ -f ${stats_4D} ]] && rm ${stats_4D} # remove any previous analysis
#[[ -f ${group_stat_dir}_tstat1.nii.gz ]] && rm ${group_stat_dir}_tstat1.nii.gz # remove any previous analysis
#[[ -f ${group_mask} ]] && rm ${group_mask} # remove any previous analysis
#[[ -f ${group_stat_dir}_dcoeff.nii.gz ]] && rm ${group_stat_dir}_dcoeff.nii.gz # remove any previous analysis


for task_contrast_pair in "${task_contrast_pairs[@]}"; do
    read -r task cope <<< "$task_contrast_pair"

    echo "Getting subject statistic filenames for task ${task} cope ${cope}."
    
    # Get all statistic filenames - subsequently used to find/replace to create mask & motion filenames

    stat_filenames=($stat_dir/${task}_$cope/*/MNINonLinear/Results/tfMRI_${task}/*/cope${cope}.feat/stats/tstat1.nii.gz)

    
    # Calculate stats, first filtering by motion or including as confound as specified

    for motion_method in "${motion_method_params[@]}"; do

        # Define group intermediates and results dir and files
    
        group_stat_dir="${intermediates_dir}/${task}_cope${cope}_motion-${motion_method}"
        motion_log="${group_stat_dir}/all_motion.txt"
        missing_motion_log="${group_stat_dir}/missing_motion.txt"
        missing_masks_log="${group_stat_dir}/missing_masks.txt"
        high_motion_log="${group_stat_dir}/high_motion.txt"
        results_effect_size="${results_dir}/hcp_act_d_${task}_rest__v2__motion-${motion_method}.nii.gz"
        
        mkdir -p $group_stat_dir/

        # Motion-specific setup

        if [ "$motion_method" != "none" ]; then
            echo "Doing motion-specific setup: ${motion_method}."
            # motion, which has a slightly different path
            base_dirs_motion=("${stat_filenames[@]%%/MNINonLinear/*}")
            base_dirs_motion=("${base_dirs_motion[@]/activation\/${task}_$cope/motion_average}")

            # we are going to remove subjects who are missing motion
            
            new_stat_filenames=()
            motion=()
            #motion_subids=()

            [[ -f ${motion_log} ]] && rm ${motion_log} # remove any previous logs since we append data
            [[ -f ${missing_motion_log} ]] && rm ${missing_motion_log} # remove any previous logs since we append data
                 
            for base_dir in "${base_dirs_motion[@]}"; do
                file=("$base_dir/MNINonLinear/Results/tfMRI_${task}/Movement_RelativeRMS_mean.txt")
                this_subid="$(basename "$base_dir")"
                #motion_subids+="$this_subid" # TODO: not sure this is useful for anything

                if [[ -f "$file" ]]; then
                    motion+=("$(< "$file")")
                    new_stat_filenames+=("${stat_filenames[$i]}")
                    echo "$this_subid ${motion[-1]}" >> $motion_log # recording motion and subIDs 
                else
                    echo "$this_subid" >> $missing_motion_log
                    #motion+=("NaN")
                fi
            done 

            # new filenames with motion=NaN files removed
            
            stat_filenames=("${new_stat_filenames[@]}")

            # Do motion strategy-specific setup
            if [ "$motion_method" == "regression" ]; then
                # Regression: create design matrix with confounds
                printf "1 %s\n" "${motion[@]}" > $design_matrix
                printf "1 0" > $design_contrast
                Text2Vest $design_matrix $design_matrix
                Text2Vest $design_contrast $design_contrast


            elif [ "$motion_method" == "threshold" ]; then
                
                # we are going to remove subjects who are high motion
                new_stat_filenames=()
                
                [[ -f ${high_motion_log} ]] && rm ${high_motion_log} # remove any previous logs since we append data

                for ((i=0; i<${#motion[@]}; i++)); do
                    if (( $(echo "${motion[$i]} <= $low_motion_threshold" | bc -l) )); then
                        new_stat_filenames+=("${stat_filenames[$i]}")
                    else
                        echo "${stat_filenames[$i]%%/MNINonLinear/*}" >> $high_motion_log # note: this is another way of getting subids
                    fi
                done

                # new filenames with high motion files removed

                stat_filenames=("${new_stat_filenames[@]}")
            
            fi

        fi        
        
        echo "Calculating stats."

        # Set masks filenames based on stats filenames and check each mask exists
        mask_filenames=("${stat_filenames[@]/tstat1.nii.gz/mask.nii.gz}") # TODO: should be able to define after the motion-based thresholding of stats files, right before the stats calculation, so we don't need to remove all these files too
        for i in "${!mask_filenames[@]}"; do
            if [[ ! -f "${mask_filenames[$i]}" ]]; then
                echo "${mask_filenames[$i]}" >> $missing_masks_log
                unset 'mask_filenames[$i]'
            fi
        done
        mask_filenames=("${mask_filenames[@]}")

        # Create group mask (this is in the loop in case subjects are removed for thresholding)
        
        fslmerge -t ${masks_4D} ${mask_filenames[@]}
        fslmaths ${masks_4D} -Tmin ${group_mask}

        # Calculate t-stats
        
        group_stat_prefix="${group_stat_dir}/ttest1"
        fslmerge -t ${stats_4D} ${stat_filenames[@]}
        if [ "$motion_method" == "regression" ]; then
            randomise -i ${stats_4D} -m ${group_mask} -d $design_matrix -t $design_contrast -o ${group_stat_prefix} >> ${temp_log}
        else
            randomise -i ${stats_4D} -m ${group_mask} -1 -o ${group_stat_prefix} >> ${temp_log}
        fi

        # Calculate Cohen's d
        
        sqrt_n=$(echo "sqrt(${#stat_filenames[@]}-1)" | bc)
        fslmaths ${group_stat_prefix}_tstat1.nii.gz -div $sqrt_n ${group_stat_prefix}_dcoeff.nii.gz
        
        mv ${group_stat_prefix}_dcoeff.nii.gz $results_effect_size
        
    done

done




