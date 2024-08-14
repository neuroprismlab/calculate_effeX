#!/bin/bash
#set -x
####################################  
#
# average RL and LR motion
#
####################################  


# Set paths and params

export PATH=$(dirname "$0"):$PATH
source setparams_activation.sh


# Setup

mkdir -p "${mean_motion_dir}"


# Average

for task_contrast_pair in "${task_contrast_pairs[@]}"; do
    
    # note that we are only using the specified task, not the contrast
    read -r task cope <<< "$task_contrast_pair"

    motion_filenames_RL=($motion_dir/*/MNINonLinear/Results/tfMRI_${task}_RL/Movement_RelativeRMS_mean.txt)
   
    echo "Averaging data for each of ${#motion_filenames_RL[@]} subjects ($task task)."

    for ((i=0; i<${#motion_filenames_RL[@]}; i++)); do
   
        this_file_RL=${motion_filenames_RL[$i]}
        this_file_LR=${motion_filenames_RL[$i]//RL/LR}

        # check if the LR file exists
        if [[ -f $this_file_LR ]]; then
            
            base_dir=$(dirname "$this_file_RL")
            avg_dir=${base_dir/$motion_dir/${mean_motion_dir}}
            avg_dir=${avg_dir/_RL/}
            this_file_avg="$avg_dir/Movement_RelativeRMS_mean.txt"

            mkdir -p "$avg_dir"
            paste <(awk '{print $1}' "$this_file_RL") <(awk '{print $1}' "$this_file_LR") | awk '{print ($1+$2)/2}' > "$this_file_avg"
        
        fi
    done

    echo "Done with $task task."

done


