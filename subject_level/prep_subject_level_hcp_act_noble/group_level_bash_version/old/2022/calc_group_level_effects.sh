#!/bin/bash
##########################################################
#
# Estimates Cohen's d map
# Usage: calc_group_level_effects.sh config/cfg.sh
# Make sure instance meets memory requirements
#
##########################################################

############## SETUP ##############

clear
[[ ! -z $1 && -f $1 ]] && source $1 || { echo "Error: Config file needed." ; exit 1 ; }


########## ESTIMATE EFFECT SIZES #########

# Create group level results (unless already done)
if [ ! -f $groupLevelTstat ]; then
    echo $groupLevelTstat
    printf "Processing (second level)...\n"
    . $stat_script
else
    printf "Result file $groupLevelTstat already exists, moving on. Delete this file if you want to re-calculate t-statistic file.\n"
fi

# Create Cohen's D map
# t -> Cohen's D: D=2*t/sqrt(DOF) , where DOF=(n-1) for a one-sample t-test
# TODO: generalize to one- or two-sample
if [ ! -f $groupLevelDcoeff ]; then
    nSubs=$(wc -l < $subNamesWithData) # count subs in dir
    sqrt_DOF=$(echo "sqrt($nSubs-1)" | bc)
    3dcalc -a $groupLevelTstat -expr '2*a/'"$sqrt_DOF" -prefix $groupLevelDcoeff
else
    printf "Effect size file already exists: $groupLevelDcoeff . Delete this file if you want to re-calculate d-coefficient.\n"
fi


########### NEW: ESTIMATE NETWORK-LEVEL EFFECT SIZES - in progress #########
#
#if [ $doNetworkSummary = true ]; then
#    # Create second level results (unless already done)
#    if [ ! -f $groupLevelTstatNet ]; then
#        echo $groupLevelTstatNet
#        printf "Processing (second level)...\n"
#        mkdir $localDirNet
#            . $scriptsDir/do_second_level__FLAME_net.sh # TODO: will have to make separate script
#    else
#        printf "Result file $groupLevelTstatNet already exists, moving on.\n"
#    fi
#
#    # Create effect size map
#    # t -> Cohen's D: D=2*t/sqrt(DOF) , where DOF=(n-1) for a one-sample t-test
#    if [ ! -f $groupLevelDcoeffNet ]; then
#        sqrt_DOF=$(echo "sqrt($nSubs-1)" | bc)
#        3dcalc -a $groupLevelTstatNet -expr '2*a/'"$sqrt_DOF" -prefix $groupLevelDcoeffNet
#    else
#        printf "Using existing $groupLevelDcoeffNet .\n"
#    fi
#fi
#

printf "Finished.\n\n"



