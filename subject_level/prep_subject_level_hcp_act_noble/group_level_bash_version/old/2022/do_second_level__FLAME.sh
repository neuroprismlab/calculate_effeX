#!/bin/bash
######################################################################
#
# This script is part of the Cluster Power Failure project
#
# Details: Runs second level parametric inference
# Usage: Called from get_data_and_ground_truth.sh
# Note: Check memory requirements for large datasets 
#
######################################################################

############# SETUP #############

# Set up design template
cp $designTemplate $designFile

sed -i "s#Xoutput_file#${statFilePrefix}#g" $designFile  
sed -i "s#Xnum_inputs#${nSubs}#g" $designFile

# ...continuing with multi-line edits
for ((subject=1; subject<=$nSubs; subject++)); do
    subID=$(sed "${subject}q;d" $subNamesWithData)
    sed -i "s#Xinput_files#set feat_files($subject) \"${localDir_lowerLevel}/${subID}_${copeFileSuffix}\"\nXinput_files#g" $designFile
    sed -i 's#XEV_vals#set fmri(evg'"${subject}"'.1) 1\nXEV_vals#g' $designFile
    sed -i "s#Xgroup_membership#set fmri(groupmem.$subject) 1\nXgroup_membership#g" $designFile
done

sed -i "s#Xinput_files##g" $designFile
sed -i "s#XEV_vals##g" $designFile
sed -i "s#Xgroup_membership##g" $designFile

############# ANALYSIS #############

# Run analysis and clean up
printf "\n++ Processing data - positive and negative contrasts... "
feat $designFile
printf "Done second level.\n"

## Copy and change design template for negative contrast
#cp $designFile $designFile_Neg
#sed -i 's#set fmri(evg\(.*\).1) 1#set fmri(evg\1.1) -1#g' $designFile_Neg
#sed -i "s#${statFilePrefix}#${statFilePrefix_Neg}#g" $designFile_Neg
# feat $designFile_Neg
#


