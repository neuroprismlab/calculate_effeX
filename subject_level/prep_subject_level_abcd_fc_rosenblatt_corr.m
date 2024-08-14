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

study_info.dataset = 'abcd';
study_info.map = 'fc'; % 'fc': functional connectivity, 'act': activation map
study_info.test = 'r'; % 'r', 't', or 't2'


% 2. Brain

% load brain data
% [subs_brain,d,msk] = load('tmp_brain_data_path');
load('/work/neuroprism/effect_size/data/subject_level/to_be_organized/s_abcd_fc_rosenblatt_corr/abcd_data_for_steph_july_30_2024.mat')

% add to master variables
task = 'rest'; % the task condition under which brain data was collected, e.g., 'REST', 'GAMBLING'
brain_data.(task).sub_ids = tb.src_subject_id; 
brain_data.(task).data = edges; % last dim is n_sub long
study_info.mask = triu(ones(268),1);


% 3. Motion
% TODO: consider separate variable for motion

% load motion data
% [subs_mot, mot] = load('tmp_motion_data_path');

% add to master variable
brain_data.(task).sub_ids_motion = tb.src_subject_id;
brain_data.(task).motion = motion_vals;


% 4. Outcome

% load outcome data
% [subs_outcome, sc, con,cat] = load('tmp_outcome_data_path');

% add to master variable

score_labels=tb_vars_of_interest;
outcome = struct();
for i = 1:length(score_labels)

    this_score_label = score_labels{i};

    % 4. Outcome

    % add to master variable

    this_test=['test',num2str(i)];
    outcome.(this_test).sub_ids = tb.src_subject_id;
    outcome.(this_test).score = tb.([this_score_label]);
    outcome.(this_test).score_label = this_score_label;
    outcome.(this_test).reference_condition = task;
    outcome.(this_test).contrast = NaN; % for contrasts using this_brain_data fields

    if i==1 | i==3
        outcome.(this_test).category = 'demographic'; % cognitive, biometric, etc.
    elseif i==4
        outcome.(this_test).category = 'biometric'; % cognitive, biometric, etc.
    elseif i==2
        outcome.(this_test).category = 'cognitive'; % cognitive, biometric, etc.
    else
        outcome.(this_test).category = 'psychiatric'; % cognitive, biometric, etc.
    end

end



% 5. Save Data

out_path = [out_dir,strjoin({prefix,study_info.dataset,study_info.map,'rosenblatt'},'_'),'.mat'];
save(out_path, 'study_info', 'brain_data', 'outcome')


