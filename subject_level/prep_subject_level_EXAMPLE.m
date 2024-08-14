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

out_dir = '/work/neuroprism/effect_size/data/subject_level/';
prefix = 's'; % to mark the file type

% 1. Study Info

study_info.dataset = 'hcp';
study_info.map = 'fc'; % 'fc': functional connectivity, 'act': activation map
study_info.test = 'r'; % 'r', 't', or 't2'


% 2. Brain

% load brain data
% [subs_brain,d,msk] = load('tmp_brain_data_path');

% add to master variables
task = 'tmp'; % the task condition under which brain data was collected, e.g., 'REST', 'GAMBLING'
brain_data.(task).sub_ids = subs_brain; 
brain_data.(task).data = d; % last dim is n_sub long
study_info.mask = msk;


% 3. Motion
% TODO: consider separate variable for motion

% load motion data
% [subs_mot, mot] = load('tmp_motion_data_path');

% add to master variable
brain_data.(task).sub_ids_motion = subs_mot;
brain_data.(task).motion = mot;


% 4. Outcome

% load outcome data
% [subs_outcome, sc, con,cat] = load('tmp_outcome_data_path');

% add to master variable
this_test=['test',num2str(i)]; % save a separate field for each unique test
score_label = 'tmp'; % summarize the test to be conducted, e.g., 'male_vs_female', 'EMOTION_vs_REST', 'gF_corr'
outcome.(this_test).sub_ids = subs_outcome;
outcome.(this_test).score = sc; % continuous, integers (ordinal), binary (for t-test), or strings/factors (categorical)
outcome.(this_test).score_label = score_label;
outcome.(this_test).reference_condition = task;
outcome.(this_test).contrast = con; % for t-test between brain_data conditions: cell array, one contrast per column
outcome.(this_test).category = cat; % cognitive, biometric, etc.


% 5. Save Data

out_path = [out_dir,strjoin({prefix,study_info.dataset,study_info.map,'AUTHOR','1'},'_'),'.mat'];
save(out_path, 'study_info', 'brain_data', 'outcome')


