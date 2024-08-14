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
%           - category ---------- demographic, cognitive, biometric, psychiatric, etc.
%       - <name of outcome 2>
%           ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Setup

data_dir= '/work/neuroprism/data_shared/hbn/';
out_dir = '/work/neuroprism/effect_size/data/subject_level/';
prefix = 's';

score_labels = {'diag','nih','srs'};
outcome_categories = {'psychiatric','cognitive','psychiatric'};


% 1. Study Info

study_info.dataset = 'hbn';
study_info.map = 'fc'; % 'fc': functional connectivity, 'act': activation map
study_info.test = 'r'; % 'r', 't', or 't2'


% Load all data

outcome = struct();
for i = 1:length(score_labels)

    this_score_label = score_labels{i};
    S = load([data_dir,this_score_label,'_files.mat']);
    M = load([data_dir,this_score_label,'_motion.mat']);


    % 2. Brain

    % add to master variables
    condition = ['rest_',this_score_label]; % the task condition under which brain data was collected, e.g., 'REST', 'GAMBLING'
    brain_data.(condition).sub_ids = S.connectome_sublist; 

    % make mask - upper triangle
    msk = ones(size(S.mean_connectome(:,:,1)));
    msk = logical(triu(msk,1));
    study_info.mask = msk;

    % apply mask to all matrices
    n = size(S.mean_connectome,3);
    brain_data.(condition).data = zeros(nnz(msk), n);
    for j = 1:n
        tmp = S.mean_connectome(:,:,j);
        brain_data.(condition).data(:,j) = tmp(msk);
    end


    % 3. Motion

    % add to master variable
    brain_data.(condition).sub_ids_motion = S.connectome_sublist;
    brain_data.(condition).motion = M.([this_score_label,'_motion'])(:,3); % TODO: replace diag_motion


    % 4. Outcome

    % add to master variable

    this_test=['test',num2str(i)];
    outcome.(this_test).sub_ids = S.connectome_sublist;
    outcome.(this_test).score = S.([this_score_label,'_behav']);
    outcome.(this_test).score_label = this_score_label;
    outcome.(this_test).reference_condition = condition;
    outcome.(this_test).contrast = NaN; % for contrasts using this_brain_data fields
    outcome.(this_test).category = outcome_categories{i}; % cognitive, biometric, etc.
end


% 5. Save Data

out_path = [out_dir,strjoin({prefix,study_info.dataset,study_info.map,study_info.test},'_'),'.mat'];
save(out_path, 'study_info', 'brain_data', 'outcome')


