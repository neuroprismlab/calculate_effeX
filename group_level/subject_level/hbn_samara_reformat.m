% reformat HBN Rest and Movie data contributed by Ahmad Samara

% load data
data = load('/Users/shearer.h/Downloads/hbn_fc_1movie_2rest.mat');

% dataset information
dataset_name = 'hbn_samara';
map_type = 'fc'; % fc or act
n_nodes = 1000;

% ref_condition = 'movie';

% paths
output_dir = '/Users/shearer.h/Library/CloudStorage/GoogleDrive-halleerenate@gmail.com/My Drive/braineffex_data/HBN_samara/';
output_file = [output_dir, dataset_name, '_', map_type, '_samara.mat'];
bids_dir = '/Users/shearer.h/neurodesktop-storage/HBN_BIDS/';
do_group_level_out_path = '/Users/shearer.h/Library/CloudStorage/GoogleDrive-halleerenate@gmail.com/My Drive/braineffex_data/results/';


% set the outcome variables to use
outcome_vars(1).full   = 'age';
outcome_vars(1).short  = 'age';
outcome_vars(1).category = 'age (demographic)';
outcome_vars(2).full   = 'sex';
outcome_vars(2).short  = 'sex';
outcome_vars(2).category = 'sex (demographic)';
outcome_vars(3).full   = 'bmi';
outcome_vars(3).short  = 'bmi';
outcome_vars(3).category = 'biometric';
outcome_vars(4).full   = 'p_factor_mcelroy_harmonized_all_samples';
outcome_vars(4).short  = 'p_factor';
outcome_vars(4).category = 'psychiatric';
outcome_vars(5).full   = 'internalizing_mcelroy_harmonized_all_samples';
outcome_vars(5).short  = 'internalizing';
outcome_vars(5).category = 'psychiatric';
outcome_vars(6).full   = 'externalizing_mcelroy_harmonized_all_samples';
outcome_vars(6).short  = 'externalizing';
outcome_vars(6).category = 'psychiatric';
outcome_vars(7).full   = 'attention_mcelroy_harmonized_all_samples';
outcome_vars(7).short  = 'attention';
outcome_vars(7).category = 'psychiatric';

% fill in study_info
template.study_info.dataset = char(dataset_name);
template.study_info.map = char(map_type); % functional connectivity

% read in a subject ID list and clean up
subs = data.hbn_id;

n_edges = n_nodes*n_nodes;
n_edges_in_tri = n_nodes * (n_nodes - 1) / 2;
n_subs = length(subs);


%%%% BRAIN DATA

% Movie data (condition 1)
% initialize brain data
brain_data = zeros([n_edges_in_tri, n_subs]);
% create mask for upper triangle (will reuse for rest)
triu_mask = triu(true(n_nodes), 1);
% save mask to template
template.study_info.mask = triu_mask;

% identify edges involving known-bad parcels (533, 903 have no data for any subject)
% known problem with Schaefer 1000, see https://github.com/ThomasYeoLab/CBIG/issues/10
[row_idx, col_idx] = find(triu_mask);
bad_parcels = [533, 903];
bad_edge_idx = find(ismember(row_idx, bad_parcels) | ismember(col_idx, bad_parcels));

% loop through subjects and add their FC to brain data
for i = 1:n_subs
this_mat = data.hbn_fc(:,:,i,1);
brain_data(:,i) = this_mat(triu_mask);
end

% zero-fill known-bad parcels' edges
brain_data(bad_edge_idx, :) = 0;

% add Movie brain data to template
template.brain_data.movie.sub_ids = subs;
template.brain_data.movie.data = brain_data;


% Rest data (condition 2)
brain_data = zeros([n_edges_in_tri, n_subs]);
for i = 1:n_subs
this_mat = data.hbn_fc(:,:,i,2);
brain_data(:,i) = this_mat(triu_mask);
end

% zero-fill known-bad parcels' edges
brain_data(bad_edge_idx, :) = 0;

% add brain data to template
template.brain_data.rest.sub_ids = subs;
template.brain_data.rest.data = brain_data;


%%%% Outcome data

% read in demographic data
dem = readtable([bids_dir, 'study-', 'HBN', '_desc-participants.tsv'], 'FileType', 'text', 'Delimiter', '\t');

% loop through each outcome variable
for t = 1:length(outcome_vars)
    var = outcome_vars(t).full;
    this_category = outcome_vars(t).category;
    var_short = outcome_vars(t).short;
    this_outcome = [];
    this_outcome_subs = strings(0);

    for idx = 1:length(subs)
        this_sub = char(subs(idx));
        this_sub_short = this_sub(5:end);
        this_idx = string(dem.participant_id) == this_sub_short;
        this_sub_outcome = table2array(dem(this_idx,var));
        if any(this_idx)
            this_outcome = [this_outcome, this_sub_outcome];
            this_outcome_subs(end+1) = string(this_sub);
        end
    end

    template.outcome.(var_short).score = this_outcome';
    template.outcome.(var_short).score_label = var_short;
    template.outcome.(var_short).sub_ids = this_outcome_subs;
    template.outcome.(var_short).category = this_category;
    template.outcome.(var_short).contrast = NaN; 
    
    % add reference condition
    template.outcome.(var_short).reference_condition = 'rest';
end

% add another outcome test for movie vs. rest
template.outcome.movie_rest.score = NaN;
template.outcome.movie_rest.sub_ids = NaN;
template.outcome.movie_rest.score_label = NaN;
template.outcome.movie_rest.reference_condition = NaN;
template.outcome.movie_rest.contrast = {'rest', 'movie'};
template.outcome.movie_rest.category = 'task connectivity';

%%%% MOTION DATA
% getting motion data from RBC HBN data
cpac_dir = '/Users/shearer.h/data/HBN_XCP-D/';

% rest condition is FC from the two rest runs concatenated
% movie condition is FC from the two movie runs concatenated
% so I will get FD from each run per condition, and concatenate before
% obtaining the mean FD

% rest motion
cond = 'rest';

motion = [];
subs_motion = [];
count = 1;

for idx = 1:height(subs)
    sub_path = [cpac_dir, char(subs(idx, :))];
    ses_folders = dir(fullfile(sub_path, 'ses-*'));
    if ~isempty(ses_folders)

        % get file pattern for run 1
        file_pattern_r1 = fullfile(sub_path, 'ses-1', 'func', ...
            [char(subs(idx, :)), '_', 'ses-1', '_task-', cond, ...
            '_run-1_motion.tsv']);

        % find file path for run 1
        matched_file = dir(file_pattern_r1);

        % get FD from run 1
        if ~isempty(matched_file)
            this_file_name = fullfile(matched_file(1).folder, matched_file(1).name);
            %subs_motion = [subs_motion, subs(idx)];
            this_motion = readtable(this_file_name, 'FileType','text', 'Delimiter', '\t');
            this_fd = this_motion.framewise_displacement;
            % this_mean_fd = mean(this_fd.framewise_displacement);
            % motion = [motion, this_mean_fd];
            % count = count + 1;
            % if there's more than 1 NaN in the FD, print how many
            if sum(isnan(this_fd), 'all') > 1
                fprintf('number of NaNs in %s: %d\n', file_pattern_r1, sum(isnan(this_fd), 'all'));
            end
        else
            fprintf('No file found for run 1: %s\n', file_pattern_r1);
        end

        file_pattern_r2 = fullfile(sub_path, 'ses-1', 'func', ...
            [char(subs(idx, :)), '_', 'ses-1', '_task-', cond, ...
            '_run-2_motion.tsv']);

        matched_file_r2 = dir(file_pattern_r2);

        % get FD from run 2
        if ~isempty(matched_file_r2)
            this_file_name_r2 = fullfile(matched_file_r2(1).folder, matched_file_r2(1).name);
            %subs_motion_r2 = [subs_motion_r2, subs(idx)];
            this_motion_r2 = readtable(this_file_name_r2, 'FileType','text', 'Delimiter', '\t');
            this_fd_r2 = this_motion_r2.framewise_displacement;
            % this_mean_fd = mean(this_fd.framewise_displacement);
            % motion = [motion, this_mean_fd];
            % count = count + 1;
            % if there's more than 1 NaN in the FD, print how many
            if sum(isnan(this_fd_r2), 'all') > 1
                fprintf('number of NaNs in %s: %d\n', file_pattern_r2, sum(isnan(this_fd_r2), 'all'));
            end
        else
            fprintf('No file found for run 2: %s\n', file_pattern_r2);
        end

        % if FD available for both runs, add to motion
        if ~isempty(matched_file) & ~isempty(matched_file_r2)
            % concatenate this_fd and this_fd_r2
            this_fd_concat = [this_fd; this_fd_r2];

            % calculate mean FD for this sub
            this_mean_fd = mean(this_fd_concat, 'omitnan');

            motion = [motion, this_mean_fd];
            subs_motion = [subs_motion, subs(idx)];
            count = count + 1;
        end
    end
end

template.brain_data.(cond).motion = motion;
template.brain_data.(cond).sub_ids_motion = subs_motion;

% movie motion
cond = 'movie';

motion = [];
subs_motion = [];
count = 1;

for idx = 1:height(subs)
    sub_path = [cpac_dir, char(subs(idx, :))];
    ses_folders = dir(fullfile(sub_path, 'ses-*'));
    if ~isempty(ses_folders)

        % get file pattern for run 1
        file_pattern_r1 = fullfile(sub_path, 'ses-1', 'func', ...
            [char(subs(idx, :)), '_', 'ses-1', '_task-', cond, ...
            'TP_motion.tsv']);

        % find file path for run 1
        matched_file = dir(file_pattern_r1);

        % get FD from run 1
        if ~isempty(matched_file)
            this_file_name = fullfile(matched_file(1).folder, matched_file(1).name);
            %subs_motion = [subs_motion, subs(idx)];
            this_motion = readtable(this_file_name, 'FileType','text', 'Delimiter', '\t');
            this_fd = this_motion.framewise_displacement;
            % this_mean_fd = mean(this_fd.framewise_displacement);
            % motion = [motion, this_mean_fd];
            % count = count + 1;
            % if there's more than 1 NaN in the FD, print how many
            if sum(isnan(this_fd), 'all') > 1
                fprintf('number of NaNs in %s: %d\n', file_pattern_r1, sum(isnan(this_fd), 'all'));
            end
        else
            fprintf('No file found for run 1: %s\n', file_pattern_r1);
        end

        file_pattern_r2 = fullfile(sub_path, 'ses-1', 'func', ...
            [char(subs(idx, :)), '_', 'ses-1', '_task-', cond, ...
            'DM_motion.tsv']);

        matched_file_r2 = dir(file_pattern_r2);

        % get FD from run 2
        if ~isempty(matched_file_r2)
            this_file_name_r2 = fullfile(matched_file_r2(1).folder, matched_file_r2(1).name);
            %subs_motion_r2 = [subs_motion_r2, subs(idx)];
            this_motion_r2 = readtable(this_file_name_r2, 'FileType','text', 'Delimiter', '\t');
            this_fd_r2 = this_motion_r2.framewise_displacement;
            % this_mean_fd = mean(this_fd.framewise_displacement);
            % motion = [motion, this_mean_fd];
            % count = count + 1;
            % if there's more than 1 NaN in the FD, print how many
            if sum(isnan(this_fd_r2), 'all') > 1
                fprintf('number of NaNs in %s: %d\n', file_pattern_r2, sum(isnan(this_fd_r2), 'all'));
            end
        else
            fprintf('No file found for run 2: %s\n', file_pattern_r2);
        end

        % if FD available for both runs, add to motion
        if ~isempty(matched_file) & ~isempty(matched_file_r2)
            % concatenate this_fd and this_fd_r2
            this_fd_concat = [this_fd; this_fd_r2];

            % calculate mean FD for this sub
            this_mean_fd = mean(this_fd_concat, 'omitnan');

            motion = [motion, this_mean_fd];
            subs_motion = [subs_motion, subs(idx)];
            count = count + 1;
        end
    end
end

template.brain_data.(cond).motion = motion;
template.brain_data.(cond).sub_ids_motion = subs_motion;



% 
% for ref_condition = {'rest_run-1', 'movieDM'}
%     ref_condition = ref_condition{1};
%     motion = [];
%     subs_motion = [];
%     count = 1;
%     for idx = 1:height(subs)
%         sub_path = [cpac_dir, char(subs(idx, :))];
%         ses_folders = dir(fullfile(sub_path, 'ses-*'));
%         if ~isempty(ses_folders)
% 
%             % Use wildcard for acq- field
%             file_pattern = fullfile(sub_path, 'ses-1', 'func', ...
%                 [char(subs(idx, :)), '_', 'ses-1', '_task-', ref_condition, ...
%                 '_motion.tsv']);
% 
%             matched_files = dir(file_pattern);
% 
%             if ~isempty(matched_files)
%                 this_file_name = fullfile(matched_files(1).folder, matched_files(1).name);
%                 subs_motion = [subs_motion, subs(idx)];
%                 this_fd = readtable(this_file_name, 'FileType','text', 'Delimiter', '\t');
%                 this_mean_fd = mean(this_fd.framewise_displacement);
%                 motion = [motion, this_mean_fd];
%                 count = count + 1;
%             else
%                 fprintf('No file found for: %s\n', file_pattern);
%             end
%         end
%     end
%     if strcmp(ref_condition, 'rest_run-1')
%         % rename ref_condition to rest
%         ref_condition = 'rest';
%     end
%     if strcmp(ref_condition, 'movieDM')
%         ref_condition = 'movie';
%     end
%     template.brain_data.(ref_condition).motion = motion;
%     template.brain_data.(ref_condition).sub_ids_motion = subs_motion;
% end



% if subject IDs contain strings/characters, then fix that
sub_ids_1 = template.brain_data.rest.sub_ids(:);
%sub_ids_1 = template.brain_data.movie.sub_ids(:);
    
% Concatenate and get unique IDs (this returns a string array)
all_sub_ids = unique(sub_ids_1);

if ~isnumeric(all_sub_ids)
    
    % Create numeric IDs
    n_subjects = length(all_sub_ids);
    numeric_ids = (1000:(1000 + n_subjects - 1))';
    
    % Create the key struct (no table needed)
    subject_id_key = struct();
    subject_id_key.original = all_sub_ids;
    subject_id_key.numeric = numeric_ids;
    
    % Replace in brain_data
    [~, idx] = ismember(template.brain_data.rest.sub_ids, all_sub_ids);
    template.brain_data.rest.sub_ids = numeric_ids(idx);
    template.brain_data.movie.sub_ids = numeric_ids(idx);
    
    % Replace in motion data
    [~, idx] = ismember(template.brain_data.rest.sub_ids_motion, all_sub_ids);
    template.brain_data.rest.sub_ids_motion = numeric_ids(idx);
    [~, idx] = ismember(template.brain_data.movie.sub_ids_motion, all_sub_ids);
    template.brain_data.movie.sub_ids_motion = numeric_ids(idx);
    
    % Replace in outcome
    for out={outcome_vars.short}
        [~, idx] = ismember(template.outcome.(out{1}).sub_ids, all_sub_ids);
        template.outcome.(out{1}).sub_ids = numeric_ids(idx);
    end
end


%%%% SAVE
brain_data = template.brain_data;
outcome = template.outcome;
study_info = template.study_info;
save(output_file, 'brain_data', 'outcome', 'study_info', '-v7.3')


% to run group level:
do_group_level(do_group_level_out_path, output_dir, 'Testing', 0, 'NumNetworks', 7);
