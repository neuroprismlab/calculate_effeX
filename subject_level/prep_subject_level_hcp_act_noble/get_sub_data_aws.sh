#!/bin/bash
#######################################################
#
# Download HCP task activation stats + motion from AWS
# Author: Stephanie Noble
#
# Prerequisites:
#   1. if needed, create dataset account and get free tier AWS account (need for HCP but not for openneuro)
#   2. install awscli - specify local install and bin directories
#   3. if needed, run aws configure to specify public and secret keys
# 
# Usage:
#   0. if you're copying a non-negligible amount of data: connect to a compute node 
#   1. set paths
#   2. run with test=1
#   3. when ready, set test=0
#
# Next, follow up with average_RL_LR_motion.sh
#
#######################################################


############## 0. SETUP ##############

# Set paths and params

# add current dir to path and current config file
export PATH=$(dirname "$0"):$PATH
source setparams.sh


# Setup

mkdir -p $stat_dir
mkdir -p $motion_dir
mkdir -p $tmp_dir

if [ $test -eq 1 ]; then
    dry_run_str='--dryrun'
    test_report_str='would be '
else
    dry_run_str=''
    test_report_str=''
fi

#subjects=($(aws s3 ls $source_dir | awk '{print $2}'))


############## 1. Copy data ##############

# Copy stats, masks, and motion
# IMPORTANT NOTE: Yes, it would be more organized to sync everything to the same folder, but it takes an insane amount of time to sync (and sometimes does not complete) when there are existing local subdirs :(
# So we make separate folders for each data type

#for sub_ID in "${subjects[@]}"; do
for task_contrast_pair in "${task_contrast_pairs[@]}"; do
    
    read -r task cope <<< "$task_contrast_pair"
   
    this_task_dir="$stat_dir/${task}_${cope}/"
    mkdir -p $this_task_dir
   
    if [ $get_stats -eq 1 ]; then

        aws s3 sync $source_dir/ $this_task_dir/ --exclude '*' --include "*/tfMRI_${task}/*/cope${cope}.feat/stats/mask.nii.gz" $dry_run_str
        #aws s3 sync $source_dir/ $this_task_dir/ --exclude '*' --include "*/tfMRI_${task}/*/cope${cope}.feat/stats/tstat1.nii.gz" --include "*/tfMRI_${task}/*/cope${cope}.feat/stats/mask.nii.gz" $dry_run_str
        #aws s3 sync $source_dir/ $this_task_dir/ --exclude '*' --include "*/tfMRI_${task}/*/cope${cope}.feat/stats/varcope1.nii.gz" --include "*/tfMRI_${task}/*/cope${cope}.feat/stats/mask.nii.gz" $dry_run_str
        # TODO: check and run the following command to copy varcopes as well to run FLAME (mixed effects)
        # aws s3 sync $source_dir/ $this_task_dir/ --exclude '*' --include "*/tfMRI_${task}/*/cope${cope}.feat/stats/tstat1.nii.gz" --include "*/tfMRI_${task}/*/cope${cope}.feat/stats/mask.nii.gz" --include "*/tfMRI_${task}/*/cope${cope}.feat/stats/varcope.nii.gz" $dry_run_str
        #aws s3 sync $source_dir/ $this_task_dir/ --exclude '*' --include "*/tfMRI_${task}/*/cope${cope}.feat/stats/mask.nii.gz" $dry_run_str
    fi

    if [ $get_motion -eq 1 ]; then
        aws s3 sync $source_dir/ $motion_dir/ --exclude '*' --include '*/tfMRI_${task}/Movement_RelativeRMS_mean.txt' $dry_run_str
    fi
done

printf "\nData ${test_report_str}copied to ${stat_dir}\n\n"


# Save in format expected for EffeX sub-level
# TODO
# want: sub-level stats files?, study info incl n, mask?? (results), study info incl n - get from here
# TODO: apply mask (results) to tstat (results) and pval (results), create study_info incl n (get from here)
# there, as for the other maps, we will do the conversion: t -> Cohen's D: D=t/sqrt(N) for a one-sample tt
