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
%           - contrast ---------- 2D cell array for t-test between brain data conditions
%           - category ---------- demographic, cognitive, biometric, psychiatric, etc.
%       - <name of outcome 2>
%           ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Setup

data_filename = '/work/neuroprism/data_shared/ukb/ukb_data_steph.mat';
out_dir = '/work/neuroprism/effect_size/data/subject_level/';
prefix = 's';

S = load(data_filename);


% 1. Study Info

study_info.dataset = 'ukb';
study_info.map = 'fc'; % 'fc': functional connectivity, 'act': activation map
study_info.test = 'r'; % 'r', 't', or 't2'


% 2. Brain

% add to master variables
condition = 'rest'; % the task condition under which brain data was collected, e.g., 'REST', 'GAMBLING'
brain_data.(condition).sub_ids = S.subject_id; 
brain_data.(condition).data = S.fc; % last dim is n_sub long

% make mask - upper triangle
n = (1 + sqrt(1 + 8 * size(S.fc,2))) / 2;
msk = ones(n,n);
msk = triu(msk,1);
study_info.mask = msk;


% 3. Motion

% add to master variable
brain_data.(condition).sub_ids_motion = S.subject_id;
brain_data.(condition).motion = S.framewise_displacement;


% 4. Outcome

% add to master variable
outcome_labels = {'age','fluid_intelligence','gender'};
outcome_categories = {'demographic','cognitive','demographic'};
outcome = struct();
for i = 1:length(outcome_labels)
    this_test=['test',num2str(i)];
    this_outcome_label = outcome_labels{i};
    outcome.(this_test).sub_ids = S.subject_id;
    outcome.(this_test).score = S.(this_outcome_label);
    outcome.(this_test).score_label = this_outcome_label;
    outcome.(this_test).reference_condition = condition;
    outcome.(this_test).contrast = NaN; % for contrasts using this_brain_data fields
    outcome.(this_test).category = outcome_categories{i}; % cognitive, biometric, etc.
end


% 5. Save Data

out_path = [out_dir,strjoin({prefix,study_info.dataset,study_info.map,'jiang'},'_'),'.mat'];
save(out_path, 'study_info', 'brain_data', 'outcome')


