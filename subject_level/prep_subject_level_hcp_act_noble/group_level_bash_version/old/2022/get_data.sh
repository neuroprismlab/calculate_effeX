#!/bin/bash
#set -x
##########################################################
#
# Downloads all cope files
# Usage: get_data.sh config/cfg.sh
#
##########################################################


############## SETUP ##############

clear
[[ ! -z $1 && -f $1 ]] && source $1 || { echo "Error: Config file needed." ; exit 1 ; }

# add s3cmd to path
if [[ ${PATH} != *"${s3cmdDir}"* ]]; then
export PATH=$PATH:${s3cmdDir}
fi

mkdir -p $localDir_lowerLevel 


######### CREATE SUBJECT NAMES FILE ###########

# Get names of all HCP subjects (unless already done)
printf "\nGetting subject data.\n"
skip='no'
if [[ -f $subNames ]]; then
    read -p "Subject names file $subNames exists - overwrite?"
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm $subNames
    else
        skip='yes'
        printf "Using existing subject names file.\n"
    fi  
fi

if [ $skip == 'no' ]; then
    s3cmd -c $s3cmdConfigFile ls $remoteDir_prefix/ > $subNames
    sed -e 's#                       DIR   '$remoteDir_prefix'##g' -i $subNames
    sed -e 's#/##g' -i $subNames
fi

# Get names of all subjects that have this specific task data (unless already done)
skip='no'
if [[ -f $subNamesWithData ]]; then
    read -p "Subjects-with-data file $subNamesWithData exists - overwrite?" 
    if [[ $REPLY =~ ^[Yy]$ ]]; then    
        rm $subNamesWithData
    else 
        skip='yes'
        printf "Using existing subjects-with-data file.\n"
    fi
else
    touch $subNamesWithData
fi

if [ $skip == 'no' ]; then
    while read subject; do
        s3cmd -c $s3cmdConfigFile ls $remoteDir_prefix/$subject/$remoteDir_suffix/$copeFileSuffix >> $subNamesWithData
    done < $subNames 
    sed -e 's#                       DIR   '$remoteDir_prefix/'##g' -i $subNamesWithData
    sed -e 's#'/$remoteDir_suffix/$copeFileSuffix/'##g' -i $subNamesWithData
fi



############ COPY DATA ############

# Check whether any data exists locally / whether to remove existing data
if [ ! "$(ls -A ${localDir_lowerLevel})" ] ; then # check whether empty
    get_data=true
else
    printf "Data exists in ${localDir_lowerLevel}."
    read -rep $'\nEnter \"delete all previous data\" to delete all previously downloaded data and download again and any other character to return. Response?   ' response
    if [[ "$response" =~ "delete all previous data" ]]; then
        printf "Okay, deleting previously created data for these jobs.\n"
        rm -r ${localDir_lowerLevel}
        mkdir ${localDir_lowerLevel}
    else
        printf "Okay, keeping previous data.\n"
        get_data=false
    fi
fi


# Copy data locally for each subject (unless data exists)
if [ "$get_data" = true ] ; then
    printf "Getting data.\n"
    while read subject; do
        printf "Downloading subject $subject \n"
        mkdir -p $localDir_lowerLevel/${subject}_$copeFileSuffix/
        echo s3cmd -c $s3cmdConfigFile get --recursive --skip-existing $remoteDir_prefix/$subject/$remoteDir_suffix/$copeFileSuffix/ $localDir_lowerLevel/${subject}_$copeFileSuffix/ || exit
    done < $subNamesWithData
fi


########## NEW: SUMMARIZE BY NETWORKS ########

## Summarize subjects by network (unless already done)
#if [ $doNetworkSummary = true ]; then
#   
#    mkdir $outputNetDir
#    mkdir $outputNetDir_Subs
#    printf "Creating network-level summaries for subject "
#    
#    while read subject; do
#        printf "$subject "
#        thisSubLowerLevel="${localDir_lowerLevel}/${subID}_cope${copeNum}.feat/stats/cope1.nii.gz"
#        thisSubNetworkSummary="$outputNetDir_Subs/${subject}_$netSuffix" # TODO
#        if [ ! -f $thisSubNetworkSummary ]; then
#            printf "1\n$thisSubLowerLevel" > $tmpImgFilename
#            . $scriptsDir/CompROI.sh -meanonly $atlas "$roiList" $tmpImgFilename $thisSubNetworkSummary
#        else
#            printf "(skipping--already exists) "
#        fi
#    done < $subNamesWithData
#fi
#
printf "Finished.\n\n"



