%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Prepares study subject-level data and meta-data for Brain EffeX
% HCP task matrices
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
%           - contrast ---------- 2D cell array for t-test between brain data conditions
%           - category ---------- demographic, cognitive, biometric, psychiatric, etc.
%       - <name of outcome 2>
%           ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Setup - USER-DEFINED

% add scripts to path: atlas tools
%addpath(genpath('~/scripts/NBS_benchmarking/NBS_benchmarking/support_scripts/atlas_scripts/'))

% set paths and strings
data_dir='~/project/HCP_S1200/';
out_dir = '/work/neuroprism/effect_size/data/subject_level/';
prefix = 's';

task_filter={'REST','EMOTION','SOCIAL','RELATIONAL','GAMBLING','WM'}; % TODO: add the rest % TODO: should it be REST1 (used by motion) or REST (used for data)?? 

% study metadata - used for filename and study info
study_info.dataset = 'hcp';
study_info.map = 'fc'; % "fc": functional connectivity, "act": activation map
study_info.test = 't'; % r, t, or t2

% outcome info
for i=1:length(task_filter)-1
    this_test=['test',num2str(i)];
    outcome.(this_test).sub_ids=NaN; % for mass univariate correlation; otherwise set NaN
    outcome.(this_test).score=NaN; % for mass univariate correlation; numerical for continuous or ordinal; string for factor; otherwise set NaN
    outcome.(this_test).score_label=NaN;
    outcome.(this_test).reference_condition=NaN;
    outcome.(this_test).contrast{1}=task_filter{1}; % for contrasts using brain_data fields (A, in A>B)
    outcome.(this_test).contrast{2}=task_filter{i+1}; % for contrasts using brain_data fields (B, in A>B)
    %outcome.(this_test).contrast(2,:)=repmat({'REST'},1,length(task_filter)-1); % for contrasts using brain_data fields (B, in A>B)
    outcome.(this_test).category='cognitive'; % cognitive, biometric, etc.
end

% 2. Load and organize data
% for brain data, last dimension is n_subs long

for task_idx=1:length(task_filter)
    
    this_task=task_filter{task_idx};
    fprintf('Running task %s\n',this_task);

    this_brain_data_dir=[data_dir,this_task,'/']; % TODO with task_filter
    if strcmp(this_task,'REST')
         this_motion_path=[data_dir,'/motion/rfMRI_',this_task,'1_allsubs_Movement_RelativeRMS_mean.txt'];
    else
        this_motion_path=[data_dir,'/motion/tfMRI_',this_task,'_allsubs_Movement_RelativeRMS_mean.txt'];
    end    
    % TODO with task filter
    this_motion_subids_path=[data_dir,'/motion/subids'];
    
    % 2.1. Load Matrices 

    % get matrix subject IDs
    flist=dir(this_brain_data_dir);
    flist(1:2)=[];
    n_subs=length(flist);

    % only the first time: make upper triangle mask (previously named triumask)
    if task_idx==1 
        tmp=importdata([this_brain_data_dir,flist(1).name]);
        brain_mask=triu(true(size(tmp,1)),1);
    end

    % preallocate matrices
    n_edges=sum(+(brain_mask(:)));
    this_brain_data=zeros(n_edges,n_subs);
    subids_matrices=zeros(n_subs,1);

    % load matrices
    for i=1:n_subs

        % append this subid
        subids_matrices(i)=str2double(flist(i).name(1:6));

        % load this matrix
        d=importdata([this_brain_data_dir,flist(i).name]);
        this_brain_data(:,i) = d(brain_mask);

        % print every 50 subs
        if mod(i,50)==0; fprintf('%d/%d\n',i,n_subs); end
    end

    % append to master brain_data 
    %brain_data{1,task_idx}=this_task;
    %brain_data{2,task_idx}=subids_matrices;
    %brain_data{3,task_idx}=this_brain_data;

    % 2.2. Load motion and align subids to brain data (remove subs without brain data)

    tmp=readmatrix(this_motion_path);
    this_motion_subids=readmatrix(this_motion_subids_path);
    this_motion=tmp(:,3);
    %this_motion=array2table([this_motion_subids,tmp(:,3)],'VariableNames', {'Subject','Motion'});

    % get motion for each subject with brain data, NaN is missing (prob some error if we have one without the other, esp motion without brain)
    this_motion_filtered=nan(length(subids_matrices), 1);
    [common_subids, brain_idx, this_motion_idx] = intersect(subids_matrices, this_motion_subids);
    this_motion_filtered(brain_idx) = this_motion(this_motion_idx);

    % append to master data
    %motion{task_idx}=this_motion_filtered;
    study_info.mask=brain_mask;
    brain_data.(this_task).sub_ids=subids_matrices;
    brain_data.(this_task).data=this_brain_data;
    brain_data.(this_task).motion=this_motion_filtered;

end

% 3. Save

out_path = [out_dir,strjoin({prefix,study_info.dataset,study_info.map,'noble','tasks'},'_'),'.mat'];
save(out_path,'study_info', 'brain_data', 'outcome')


