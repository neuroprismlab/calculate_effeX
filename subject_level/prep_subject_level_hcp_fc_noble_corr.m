%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Prepares study subject-level data and meta-data for Brain EffeX
% HCP rest matrices + traits
% Note 1: IMPORTANT: designed for data on mrrc server
% Note 2: only keeps subjects with both brain and trait data
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Setup

% set paths and strings
data_dir='~/project/HCP_S1200/';
out_dir = '/work/neuroprism/effect_size/data/subject_level/';
prefix = 's';

score_labels={'PMAT24_A_CR', 'Gender', 'Age', 'ProcSpeed_AgeAdj'};
traits_categories={'cognitive', 'socioeconomic', 'biometric', 'cognitive'};

this_brain_data_dir=[data_dir,'REST/'];
all_traits_path=[data_dir,'unrestricted_beh_2_4_2022_20_32_4.csv'];
motion_path=[data_dir,'motion/rfMRI_REST1_allsubs_Movement_RelativeRMS_mean.txt'];
motion_subids_path=[data_dir,'motion/subids'];

% 1. Study Info

study_info.dataset = 'hcp';
study_info.map = 'fc'; % 'fc': functional connectivity, 'act': activation map
study_info.test = 'r'; % 'r', 't', or 't2'

% brain data info
condition='REST'; % descriptive

% study metadata - used for filename and study info
% if r:
%study_info.this_brain_data_type = 'rest';
% if paired t (i.e., condition 1 vs. condition 2):
%study_info.condition_types = NaN;  % provide two conditions, e.g., study_info.study_info.condition_types = {'emotion', 'rest'}
% if unpaired t (i.e., group vs. 0):
%study_info.condition_type = NaN; % e.g., study_info.condition_type = 'emotion';
% if t2:
%study_info.group_types = NaN; %study_info.group_types = {'male', 'female'};
%study_info.this_brain_data_type = NaN; %study_info.this_brain_data_type = 'rest'; TODO: this will overwrite the above if left uncommented






% 2. Brain

% load matrices

% get matrix subject IDs
flist=dir(this_brain_data_dir);
flist(1:2)=[];
n_subs=length(flist);

% make upper triangle mask (previously named triumask)
tmp=importdata([this_brain_data_dir,flist(1).name]);
brain_mask=triu(true(size(tmp,1)),1);

% preallocate matrices
n_edges=sum(+(brain_mask(:)));
this_brain_data=zeros(n_edges,n_subs);

% load matrices
for i=1:n_subs

    % append this subid
    subids_matrices(i)=str2double(flist(i).name(1:6));

    % load this matrix
% TODO
    d=importdata([this_brain_data_dir,flist(i).name]);
    this_brain_data(:,i) = d(brain_mask);

    % print every 50 subs
    if mod(i,50)==0; fprintf('%d/%d\n',i,n_subs); end
end


% 3. Outcome

% Load (non-brain_data) trait

% load all traits + filter by specified subset
traits_all=readtable(all_traits_path);
traits_filtered = traits_all(:, [{'Subject'}, score_labels]);

% convert any strings to categorical vars
for col = 2:width(traits_filtered)  % skip first 'Subject' column
    if iscell(traits_filtered{1, col})
        this_trait=traits_filtered.Properties.VariableNames{col};
        traits_filtered.(this_trait) = categorical(traits_filtered.(this_trait));
    end
end

% save separately for each outcome

outcome = struct();
for i = 1:length(score_labels)
    this_score_label = score_labels{i};
    this_test=['test',num2str(i)];
    outcome.(this_test).sub_ids = traits_filtered.sub_ids;
    outcome.(this_test).score = traits_filtered.(this_score_label);
    outcome.(this_test).score_label = this_score_label;
    outcome.(this_test).reference_condition = condition;
    outcome.(this_test).contrast=NaN; % for contrasts using this_brain_data fields
    outcome.(this_test).category=traits_categories{i}; % cognitive, biometric, etc.
end


% 4. Motion

% Load motion and align subids to brain data (remove subs without brain data)

tmp=readmatrix(motion_path);
motion_subids=readmatrix(motion_subids_path);
motion=tmp(:,3);
%this_motion=array2table([this_motion_subids,tmp(:,3)],'VariableNames', {'Subject','Motion'});

% get motion for each subject with brain data, NaN is missing (prob some error if we have one without the other, esp motion without brain)
motion_filtered=nan(length(subids_matrices), 1);
[common_subids, brain_idx, motion_idx] = intersect(subids_matrices, motion_subids);
motion_filtered(brain_idx) = motion(motion_idx);





%{
% append to trait (if motion exists without trait, remove entry; if trait exists without motion, add NaN for motion)
traits_filtered.Motion = nan(height(traits_filtered), 1);
[common_subids, traits_idx, motion_idx] = intersect(traits_filtered.Subject, motion_subids);
traits_filtered.Motion(traits_idx) = motion(motion_idx);

% 3. Get data from subjects with both traits and matrices

subids = intersect(subids_matrices, traits_filtered.Subject);

% filter traits by common subs
traits_filtered = traits_filtered(ismember(traits_filtered.Subject, subids), :);

% filter matrices by common subs
% must first get the order dictated by subids
[~, indices] = ismember(subids, subids_matrices);
this_brain_data = this_brain_data(:, indices);
%}

% 6. Save: brain data, traits, motion, subids, brain mask, and study_info
% format of brain data and brain mask: matrix (for brain data, last dimension is n_subs long)
% format of traits, motion, and subids: either matrix or table

study_info.mask=brain_mask;
brain_data.(condition).sub_ids=subids; 
brain_data.(condition).data=this_brain_data;
brain_data.(condition).sub_ids_motion=subids; 
brain_data.(condition).motion=traits_filtered(:,size(traits_filtered,2));  
%outcome.(score_label).sub_ids=subids;
%outcome.(score_label).scores=traits_filtered(:,2:end-1);

out_path = [out_dir,strjoin({prefix,study_info.dataset,study_info.map,'noble','corr'},'_'),'.mat'];
save(out_path, 'brain_data', 'outcome','study_info')



