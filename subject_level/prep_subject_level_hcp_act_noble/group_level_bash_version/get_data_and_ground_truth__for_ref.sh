#!/bin/bash
##########################################################
#
# Estimates HCP Cohen's d activation map from sub-level data
# Derived from cluster_power_failure/empirical_study/get_data_and_ground_truth.sh
#
# Prerequisites: 
#   1. FSL
#   2. Appx. 4-8 GB mem
#
##########################################################



############## 0. SETUP ##############

# Set paths and params

# add original sub-level dir config file to path
sublevel_dir = '/work/neuroprism/effect_size/scripts/subject_level/prep_subject_level_hcp_act_noble/'
export PATH=$sublevel_dir:$PATH
source setparams.sh


# Set up - special for full dataset
nSubs_subset=$nSubs_total
subNames_subset="$subNamesWithInput"
permDir="$dataDir_localRepository"
permOutputsDir="${permDir}"
dataDir_lowerLevel="${dataDir_localRepository_lowerLevel}"
logfile="${permDir}/log"
templogfile="${permDir}/tmp"
emptyImg="${permOutputsDir}/emptyimg.nii.gz"
doRandomise=false # this overwrites the config setting - set false for ground truth calculation bc faster to use FLAME to calculate t-statistic

# ...additional things for FLAME setup
if [ "$doRandomise" == "false" ]; then
    designFile_Pos="$permDir/design.fsf"
    designFile_Neg="$permDir/design_Neg.fsf"
    outputFile_Pos="$permOutputsDir/${processedSuffix}_Pos"
    outputFile_Neg="$permOutputsDir/${processedSuffix}_Neg"
    resultImg_Pos="${outputFile_Pos}${resultImgSuffix}"
    resultImg_Neg="${outputFile_Neg}${resultImgSuffix}"
fi

mkdir -p $dataDir_lowerLevel 

# Get names of all HCP subjects (unless already done)
printf "\nGetting subject data.\n"
if [[ -f $subNames ]]; then
    read -p "Subject names file $subNames exists - overwrite?"
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm $subNames
        s3cmd -c $hcpConfigFile ls $cloudDataDir/ > $subNames
        sed -e 's#                       DIR   '$cloudDataDir'##g' -i $subNames
        sed -e 's#/##g' -i $subNames
    else
        printf "Using existing subject names file.\n"
    fi  
fi

# Get names of all subjects that have this specific task data (unless already done)
skip='no'
if [[ -f $subNamesWithInput ]]; then
    read -p "File  $subNamesWithInput exists - overwrite?" 
    if [[ $REPLY =~ ^[Yy]$ ]]; then    
        rm $subNamesWithInput
    else skip='yes'
    fi
else
    touch $subNamesWithInput
fi

if [ $skip == 'no' ]; then
    while read subject; do
        s3cmd -c $hcpConfigFile ls $cloudDataDir/$subject/$cloudDataDir_contd/$inputFileSuffix >> $subNamesWithInput
    done < $subNames 
    sed -e 's#                       DIR   '$cloudDataDir/'##g' -i $subNamesWithInput
    sed -e 's#'/$cloudDataDir_contd/$inputFileSuffix/'##g' -i $subNamesWithInput
fi


############ DOWNLOAD DATA ############

# Check whether data exists / whether to remove existing data
if [ ! "$(ls -A ${dataDir_lowerLevel})" ] ; then # check whether empty
    get_data=true
else
    printf "Data exists in ${dataDir_lowerLevel}."
    read -rep $'\nEnter \"delete all previous data\" to delete all previously downloaded data and download again and any other character to return. Response?   ' response
    if [[ "$response" =~ "delete all previous data" ]]; then
        printf "Okay, deleting previously created data for these jobs.\n"
        rm -r ${dataDir_lowerLevel}
        mkdir ${dataDir_lowerLevel}
    else
        printf "Okay, keeping previous data.\n"
        get_data=false
    fi
fi

# Copy data locally (unless already done)
if [ "$get_data" = true ] ; then
    printf "Getting data.\n"
    while read subject; do
        printf "Downloading subject $subject \n"
        mkdir -p $dataDir_lowerLevel/${subject}_$inputFileSuffix/
        s3cmd -c $hcpConfigFile get --recursive --skip-existing $cloudDataDir/$subject/$cloudDataDir_contd/$inputFileSuffix/ $dataDir_lowerLevel/${subject}_$inputFileSuffix/ || exit
    done < $subNamesWithInput
fi


########## ESTIMATE "GROUND TRUTH" #########

# Create second level results (unless already done)
if [ ! -f $groundTruthTstat ]; then
    echo $groundTruthTstat
    printf "Processing (second level)...\n"
    if [ $doRandomise = true ]; then
        . $scriptsDir/do_second_level__randomise.sh 
    else
        . $scriptsDir/do_second_level__FLAME.sh 
    fi
else
    printf "Result file $groundTruthTstat already exists, moving on.\n"
fi

# Create ground truth map
# t -> Cohen's D: D=t/sqrt(N) for a one-sample t-test
if [ ! -f $groundTruthDcoeff ]; then
    sqrt_n=$(echo "sqrt($nSubs_subset)" | bc)
    3dcalc -a $groundTruthTstat -expr 'a/'"$sqrt_n" -prefix $groundTruthDcoeff
else
    printf "Using existing $groundTruthDcoeff .\n"
fi


############## CLEANUP ##############

# Rename results folder
nSubs_total=$(wc -l < $subNamesWithInput)
dataDir_localRepository_new="$dataMasterDir/GroupSize$nSubs_total"
if [ ! -d $dataDir_localRepository_new ]; then
    mv $dataDir_localRepository $dataDir_localRepository_new
fi

printf "Finished.\n\n"



