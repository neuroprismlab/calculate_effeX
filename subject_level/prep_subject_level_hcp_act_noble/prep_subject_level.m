%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepares study subject-level data and meta-data for Brain EffeX
%
% Outputs 3 structures:
%   - study_info
%       - dataset
%       - test
%       - map
%       - brain_mask
%   - brain_data
%       - <name of condition 1>
%           - sub_ids
%           - data -------------- last dim is n_sub long
%           - sub_ids_motion
%           - motion
%       - <name of condition 2>
%           ...
%   - outcome
%       - <name of outcome 1> ----- e.g., 'male_vs_female', 'EMOTION_vs_REST', 'gF_corr'
%           - sub_ids
%           - score ------------- continuous, integers (ordinal), binary (for t-test), or strings/factors (categorical)
%           - contrast ---------- cell array, one contrast per column, for t-test between brain data conditions
%           - category ---------- cognitive, biometric, etc.
%       - <name of outcome 2>
%           ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Set paths

target_dir='/work/neuroprism/data_shared/hcp';
stat_dir=[target_dir,'/activation'];
motion_dir=[target_dir,'/motion_average'];
%task_contrast_pairs = {'GAMBLING','6'}; % TODO: add the rest and loop
task_contrast_pairs={'SOCIAL','6'; 'WM','20'; 'RELATIONAL','4';'GAMBLING','6';'EMOTION','3'};
% TODO: this could be defined by the setparams file used for the bash script

out_dir = '/work/neuroprism/effect_size/data/subject_level/';
prefix = 's'; % to mark the file type


% 1. Study Info

study_info.dataset = 'hcp';
study_info.map = 'act'; % 'fc': functional connectivity, 'act': activation map
study_info.test = 't'; % 'r', 't', or 't2'

% 2. Brain

for task_idx=1:size(task_contrast_pairs,1)
    
    clearvars -except stat_dir motion_dir task_contrast_pairs out_dir prefix study_info brain_data outcome task_idx
    % TODO: encapsulate instead

    % 2.0. Setup paths
    this_condition=task_contrast_pairs{task_idx,1}; % the task condition under which brain data was collected, e.g., 'REST', 'GAMBLING'
    this_cope=task_contrast_pairs{task_idx,2};
    this_data_dir=[stat_dir,'/',this_condition,'_',this_cope];

    fprintf('Running task %s\n',this_condition);


    % 2.1. Masks

    % find all mask files

    mask_files = dir([this_data_dir,'/*/*/*/*/*/*/stats/mask.nii.gz']); % TODO: add actual subdir structure
    %mask_files = {mask_files.name}; % TODO: make sure this gets unique filenames for each mask

    % load masks

    % preallocate
    D = niftiread(fullfile(mask_files(1).folder, mask_files(1).name));
    mask_hdr = niftiinfo(fullfile(mask_files(1).folder, mask_files(1).name));
    mask_4D = zeros([size(D), numel(mask_files)]);
    mask_4D(:,:,:,1) = D;

    for i = 2:numel(mask_files)
        mask_4D(:,:,:,i) = niftiread(fullfile(mask_files(i).folder, mask_files(i).name));
    end

    % make mask based on intersection (all==1)

    mask = all(mask_4D, 4);

    fprintf('Finished group mask.\n')

    % 2.2 Data & sub IDs

    % find all stat files

    stat_files = dir([this_data_dir,'/*/*/*/*/*/*/stats/tstat1.nii.gz']);
    %stat_files = {stat_files.name}; % TODO: make sure this gets unique filenames for each stat

    % extract sub IDs

    subids_brain = arrayfun(@(x) str2num(extractBefore(extractAfter(x.folder, [this_data_dir,'/']), '/')), stat_files);
    %subids_brain = cell2mat(subids_brain);



    % load and mask stats

    % preallocate
    tmp = niftiread(fullfile(stat_files(1).folder, stat_files(1).name));
    D = tmp(mask);
    stats_masked = zeros([size(D), numel(stat_files)]);
    stats_masked(:,1) = D;

    for i = 2:numel(stat_files)
        tmp = niftiread(fullfile(stat_files(i).folder, stat_files(i).name));
        D(:,i) = tmp(mask);
    end

    fprintf('Finished concatenating subjects.\n')

    % add to master variables

    brain_data.(this_condition).sub_ids = subids_brain; 
    brain_data.(this_condition).data = D; % last dim is n_sub long
    brain_data.(this_condition).mask = mask;
    brain_data.(this_condition).mask_hdr = mask_hdr;
    %study_info.mask = mask; % TODO: decide whether to provide a mask per condition in data, provide one overall in study_info, or either
    %niftiwrite(mask, [out_dir,task_contrast_pairs{1},'_',task_contrast_pairs{2},'group_mask.nii']);



    % 3. Motion

    % find all motion files

    %if strcmp(this_condition,'REST')
    %    motion_files = dir([motion_dir,'/*/*/*/rfMRI_',this_condition,'/Movement_RelativeRMS_mean.txt']);
    %    %this_motion_path=[data_dir,'/motion/rfMRI_',this_condition,'1_allsubs_Movement_RelativeRMS_mean.txt'];
    %else
        motion_files = dir([motion_dir,'/*/*/*/tfMRI_',this_condition,'/Movement_RelativeRMS_mean.txt']);
    %    %this_motion_path=[data_dir,'/motion/tfMRI_',this_condition,'_allsubs_Movement_RelativeRMS_mean.txt'];
    %end    

    %motion_files = {motion_files.name}; % TODO: make sure this gets unique filenames for each stat
    %this_subids_motion_path=[data_dir,'/motion/subids'];

    % load motion

    for i = 1:numel(motion_files)
        motion(i) = load(fullfile(motion_files(i).folder, motion_files(i).name));
        subids_motion(i) = str2num(extractBefore(extractAfter(motion_files(i).folder, [motion_dir,'/']), '/')); 
    end

    fprintf('Finished loading motion.\n')

    % get motion for each subject with brain data, NaN is missing (prob some error if we have one without the other, esp motion without brain)
    %motion_filtered=nan(length(subids_brain), 1);
    %[common_subids, brain_idx, motion_idx] = intersect(subids_brain, subids_motion);
    %motion_filtered(brain_idx) = motion(motion_idx);


    % add to master variable
    brain_data.(this_condition).sub_ids_motion = subids_motion;
    brain_data.(this_condition).motion = motion;



    % 4. Outcome

    % load outcome data
    % [subs_outcome, sc, con,cat] = load('tmp_outcome_data_path');

    % add to master variable
    this_test=['test',num2str(task_idx)]; % save a separate field for each unique test
    score_label = ''; % summarize the test to be conducted, e.g., 'male_vs_female', 'EMOTION_vs_REST', 'gF_corr'
    outcome.(this_test).sub_ids = NaN;
    outcome.(this_test).score = NaN; % continuous, integers (ordinal), binary (for t-test), or strings/factors (categorical)
    outcome.(this_test).score_label = '';
    outcome.(this_test).reference_condition = '';
    outcome.(this_test).contrast = {this_condition}; % for t-test between brain_data conditions: cell array, one contrast per column
    outcome.(this_test).category = 'cognitive'; % cognitive, biometric, etc.
    
    fprintf('Done.\n')

end

% 5. Save Data

out_path = [out_dir,strjoin({prefix,study_info.dataset,study_info.map,'noble','1'},'_'),'.mat'];
save(out_path, 'study_info', 'brain_data', 'outcome')


