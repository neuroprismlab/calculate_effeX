# Paths and parameters for Cluster Power Failure scripts

################### USER-DEFINED PARAMETERS ###################

testing=false

# Script directories
scriptsDir="/mridata2/home2/smn33/scripts/effect_size_repo/"
s3cmdDir="/mridata2/home2/smn33/scripts/s3cmd-2.0.2/"
s3cmdConfigFile="$scriptsDir/config/hcp_access_S1200" # S3 access config file 

# Data directory
localDir_project="/data15/mri_group/smn33_data/effect_size_repo" # local copy of data - TODO

# Data params
dataset="hcp"
hcpReleaseNo="1200"
task="SOCIAL"
# task-cope mapping: SOCIAL_cope6; WM_cope20; GAMBLING_cope6; RELATIONAL_cope4; EMOTION_cope3
case $task in
    'SOCIAL')
        copeNum="6" ;;
    'WM')
        copeNum="20" ;;
    'GAMBLING')
        copeNum="6" ;;
    'RELATIONAL')
        copeNum="4" ;;
    'EMOTION')
        copeNum="3" ;;
    *)
        echo "Error: must specify task."
        exit
esac

# NEW: do network summary (for network-level inference, added for "Leveling Up")
#doNetworkSummary=false
#atlas="$scripts/new_net_scripts/shen_1mm_268_parcellation__in_subnetworks.nii.gz" # TODO - confirm there's 10 nets, should be
#roiList="$(seq 1 10)"

# Effect size files
groupLevelPrefix="group_level" # used for naming


################# FILENAMING FOR REMOTE AND LOCAL DATA #################

# On AWS: fixed filenaming for lower-level data
# full filename is $remoteDir_prefix/$subject/$remoteDir_suffix/$copeFileSuffix
remoteDir_prefix="s3://hcp-openaccess/HCP_${hcpReleaseNo}"
remoteDir_suffix="MNINonLinear/Results/tfMRI_$task/tfMRI_${task}_hp200_s4_level2vol.feat"
copeFileSuffix="cope${copeNum}.feat" #  also used for naming local data

# Local: lower-level data
localDir_dataset="$localDir_project/$dataset"
localDir_task="${localDir_dataset}/${task}_cope${copeNum}"
localDir_lowerLevel="$localDir_task/lower_level"

# Subject name filenames
subNames="$localDir_dataset/all_sub_file_names.txt" # all subjects
subNamesWithData="$localDir_task/file_names_with_task_data.txt" # subjects with this task / cope



################# EFFECT SIZE FILES AND DIRECTORIES #################

# Intermediate folder and processed string
testingStr=$( [ $testing = "true" ] && echo "_TESTING" || echo "" )
localDir_intermediates="$localDir_task/intermediates$testingStr"

# Design and script files for getting stats
designTemplate="$scriptsDir/config/design_templates/design_template.fsf" #FLAME
designFile="$localDir_intermediates/design.fsf"
stat_script="$scriptsDir/do_second_level__FLAME.sh"

# Effect size output files
statFilePrefix="${localDir_task}/${groupLevelPrefix}" # provided to FSL
tstatSuffix=".gfeat/cope1.feat/stats/tstat1.nii.gz" # file suffix conventions by FSL (see next line for how this comes together with provided prefix)
groupLevelTstat="${statFilePrefix}${tstatSuffix}"
groupLevelDcoeff="${localDir_task}/dcoeff.nii.gz"

# NEW: for network summary
#outputNetDir="$outputDir/nets"
#outputNetDir_Subs="$outputNetDir/subs"
#tmpImgFilename="$outputNetDir/tmp_filename" # TODO
#netSuffix="cope${copeNum}_by_net.nii.gz"





