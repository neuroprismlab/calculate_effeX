#!/bin/bash
#######################################################
#
# Estimates HCP group-level activation maps
# (OLS only, no mixed effects)
# Author: Stephanie Noble
#
# Prerequisites:
#   1. FSL
#   2. Appx. 4-8 GB mem
## 
#######################################################

############## 0. SETUP ##############

# Set paths and params

# add current dir to path and current config file
export PATH=$(dirname "$0"):$PATH
source setparams.sh


# Setup

mkdir -p $tmp_dir

if [ $test -eq 1 ]; then
    dry_run_str='--dryrun'
    test_report_str='would be '
else
    dry_run_str=''
    test_report_str=''
fi

#subjects=($(aws s3 ls $source_dir | awk '{print $2}'))



############## 1. Estimate group-level effects ##############

# TODO: test

for task_contrast_pair in "${task_contrast_pairs[@]}"; do                                                 
    
    read -r task cope <<< "$task_contrast_pair"

    this_task_dir="$stat_dir/${task}_${cope}"                        
    results_prefix="$results_dir/${task}_${cope}"
    merged_copes="${tmp_dir}/merged_copes"
    group_mask="${results_prefix}_group_mask.nii.gz"

    # Merge subject-level stat files
    fslmerge -t "${merged_copes}" "${this_task_dir}/*/*/*/*/*/*/stats/tstat1.nii.gz"

    # Make group mask
    if [[ -f "${group_mask}" ]]; then
        rm "${group_mask}"
    fi
    3dMean -prefix "${group_mask}" -mask_inter "${this_task_dir}/*/*/*/*/*/*/stats/mask.nii.gz"

    # Run randomise
    randomise -i "${merged_copes}" -m "${group_mask}" -o "${results_prefix}" -1 -n 5000 -D                                                 
done


# Save in format expected for EffeX group-level
# TODO
# want: tstat (results), pval (results), mask (results), study info incl n - get from here
# TODO: apply mask (results) to tstat (results) and pval (results), create study_info incl n (get from here)
# there, as for the other maps, we will do the conversion: t -> Cohen's D: D=t/sqrt(N) for a one-sample tt
