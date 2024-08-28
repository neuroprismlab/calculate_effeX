%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate group-level statistics
% Authors: Hallee Shearer & Stephanie Noble
%
% Prereqs:
%   - [any unique parameters to set or scripts to run beforehand]
%   - [Note: you may need to create a mean motion file]
%
% Notes:
%   Motion note: For datasets processed with our pipeline, motion already regressed from timeseries data (regresses 24 motion parameters)
%   ** note on spreadsheet if not **
%
% Output:
%   - study_info
%       - dataset
%       - <test_components> (e.g., {'condition_label', 'score_label'}
%       - map
%       - test
%       - mask % TODO: there will be one for each brain condition
%       - category
%   - data
%       - <pooling strategy>
%            - <motion strategy>
%                 - b_standardized
%                 - p
%                 - std_brain
%                 - std_score
%                 - n         ------------ NaN if two-sample
%                 - n1        ------------ NaN if one-sample
%                 - n2        ------------ NaN if one-sample
%
% TODO: results_file_prefix naming convention - {dataset name}_{contributor name}_{date}
%       subject level: level1 + {date}
%       group level: study, effect map, level2, or group + {date}
%       combined: combined_effect_{date}
%
% TODO: consider adding number of subjects who are high motion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Setup

% script and data directories

%motion_dir = '/MOTIONPATH/'; % USER-DEFINED - REMOVED 051624
scripts_dir = '/home/h.shearer/hallee/calculate_effeX/group_level/helper_scripts/'; 
%scripts_dir = '/home/USERNAME/scripts/effect_size/';
results_dir = '/work/neuroprism/effect_size/data/group_level/';
%data_filename = '/work/neuroprism/data_shared/ukb/ukb_data_steph.mat'; % ADDED 051624
data_dir = '/work/neuroprism/effect_size/data/subject_level/';

% create a list of files in data_dir to loop through:
filenames = dir(data_dir);
filenames = filenames(~[filenames.isdir]);  % Remove directories from the list


% store the names of all datasets to loop through
datasets = {filenames.name}; 
% datasets(1:3) = [];
%TMP: testing
datasets= {'s_hcp_fc_noble_corr.mat'};

% results prefix for testing (will be appended to the start of each result
% file
res_prefix = 'testing';

% script paths

%data_loader_script_path = [scripts_dir,'data_loader']; % USER-DEFINED [this is only a placeholder and does not exist] - REMOVED 051624
atlas_tools_script_path = [scripts_dir,'atlas_tools']; % for pooling (structure_data, summarize_data)
regression_fast_script_path = [scripts_dir,'regression_fast']; % for fast multiple regression


% params

% if ukb data, pooling_params = [0] because we don't have a map
% TMP: change this once we have a ukb map
if contains(datasets, "ukb")
    pooling_params = [0];
else
    pooling_params = [0, 1];
end
motion_method_params = {'none', 'regression', 'threshold'};
low_motion_threshold = 0.1; % empirically, 5.7% of subjects are >0.1mm mFFD
n_network_groups = 10; % hard-coded for Shen atlas-based pooling

% set paths

%addpath(data_loader_script_path); % USER-DEFINED
addpath(genpath(atlas_tools_script_path));
addpath(regression_fast_script_path);


%% Load data and subject ID mapping

% Data structure for reference:
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
%       - test1
%           - sub_ids
%           - score ------------- continuous, integers (ordinal), binary (for t-test), or strings/factors (categorical)
%           - score_label
%           - reference_condition
%           - contrast ---------- 2D cell array, for t-test between brain data conditions
%           - category ---------- demographic, cognitive, biometric, psychiatric, etc.
%       - test2
%           ...
%
% USER-DEFINED - modify this section to load your data
% load the following:
%   - m (data matrices): should be n_sub x n_var 
%       - IMPORTANT: for summarization and pooling to work correctly, each row of m should be taken from the upper triangle of the original matrix consistent with Matlab indexing.
%         For example, this will be the case is m is created and saved in Matlab as follows:
%           triu_mask=logical(triu(ones(268),1));
%           m(this_subject,:) = mat(triu_mask);
%   - score (scoreavior / phenotype to predict): should be n_sub x 1
%   - subids_data (subject IDs for the matrices / images): should be n_sub x 1


%[m, score, subids_data] = data_loader;

% TODO: loop through outcomes, and get data according to each test / outcome name

%this_score='age';


for i = 1:length(datasets)
   
    dataset = datasets{i};
    fprintf(['Processing dataset: ',dataset,'\n'])
    
    data_path = [data_dir, dataset];

    S = load(data_path);

    % run checker on data
    S = checker(S);


    % create a struct to store the results 
    % results struct will have one struct for each contrast
    % each contrast struct will have one struct for each combination of pooling
    % and motion regression

    tests = fieldnames(S.outcome);

    for t = 1:length(tests)
        
        % for each outcome...
        test = tests{t}; % get the name of the outcome/score
        disp(['Running test "', test,'"'])
        
        % Setup

        results = []; % clear the struct from prior loops
        results.study_info.dataset = S.study_info.dataset;
        results.study_info.test = S.study_info.test;
        results.study_info.map = S.study_info.map;
        if isfield(S.study_info, 'mask') 
            results.study_info.mask = S.study_info.mask;
        end % TODO: make sure all input is mask not brain_mask
        if isfield(S.outcome.(test), 'level_map')
            results.study_info.level_map = S.outcome.(test).level_map;
        end


        test_type = infer_test_type(S, test);
        
        % TMP: skip the test if the test type is t
%         if strcmp(test_type, "t")
%             disp(['skipping test ', test, ' because test type is t'])
%             continue
%         end

        if strcmp(test_type, 'unknown')
            error(['could not infer test type for test ', test])
        end

        if test_type == 'r'
            % if study is r, then grab the reference_condition
            condition = S.outcome.(test).reference_condition;

            % use the reference condition as idx to grab ref brain data
            m = S.brain_data.(condition).data;

            % make sure the dims of m are subs x parcels
            if size(m,2) == length(S.brain_data.(condition).sub_ids)
                % if the second dimension is subjects,
                % flip dimensions
                m = m';
            end

            % get score data
            score = S.outcome.(test).score;
            score_label = S.outcome.(test).score_label;

            % get motion data
            motion = S.brain_data.(condition).motion;

            % save contrast to results
            results.study_info.test_components = {condition, score_label};
            results_file_prefix = [results_dir, res_prefix, S.study_info.dataset, '_', S.study_info.map, '_', S.study_info.test, '_', score_label, '_', condition];
            % TODO: remove 'hstest' after testing done

            % remove missing subject data
            [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);

        end

        %TODO: 
        if test_type == 't'
            
            if length(S.outcome.(test).contrast)==1
                % TODO: brain mask for t studies
                condition = S.outcome.(test).contrast{1};
                m = S.brain_data.(condition).data;
                
                % TODO: figure out pooling for activation maps
                % for now, skipping pooling
                pooling_params = [0];

                % save contrast to results
                results.study_info.test_components = {condition};

                score = ones(1,size(m,2));
                % score is the brain data for the second condition (could change
                % naming conventione here, but I was trying to keep things
                % consistent because after this step it should get treated the same
                % I think?)
                motion = S.brain_data.(condition).motion;
                % we have two motion variables here because motion is different for
                % each run! TODO: decide how to use this...

                % make sure the dims of m are subs x parcels
                if size(m,2) == length(S.brain_data.(condition).sub_ids)
                    % if the second dimension is subjects,
                    % flip dimensions
                    m = m';
                end
                
                % remove missing subject data
                [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);
                
                score = ones(1,size(m,1));

                results_file_prefix = [results_dir, res_prefix, S.study_info.dataset, '_', S.study_info.map, '_', S.study_info.test, '_', condition];

            else
                % TODO: brain mask for t studies
                % assuming that for t, contrast is an array of two strings
                % representing the two conditions in the contrast
                condition1 = S.outcome.(test).contrast{1};
                condition2 = S.outcome.(test).contrast{2};
                m = S.brain_data.(condition1).data;
                % m is the brain data for the first condition

                % save contrast to results
                results.study_info.test_components = {condition1, condition2};

                % make sure the dims of m are subs x parcels
                if size(m,2) == length(S.brain_data.(condition).sub_ids)
                    % if the second dimension is subjects,
                    % flip dimensions
                    m = m';
                end

                score = S.brain_data.(condition2).data;
                % score is the brain data for the second condition (could change
                % naming conventione here, but I was trying to keep things
                % consistent because after this step it should get treated the same
                % I think?)
                motion1 = S.brain_data.(condition1).motion;
                motion2 = S.brain_data.(condition2).motion;
                % we have two motion variables here because motion is different for
                % each run! TODO: decide how to use this...

                % remove missing subject data
                [m, score, motion1, motion2] = remove_missing_subs(m, score, S, test_type, test, condition1, condition2, motion1, motion2);
                motion=mean([motion1,motion2],2);
                    
                results_file_prefix = [results_dir, res_prefix, S.study_info.dataset, '_', S.study_info.map, '_', S.study_info.test, '_', condition1, '_', condition2];
        
            end
        end
        
        if strcmp(test_type, 't2')

            % If contrast is NaN, then the contributor has
            % combined the data and we will use outcome score as dummy variable
            contrast = S.outcome.(test).contrast;

            % Check if contrast has a length of 2, then that means the
            % contributor has not combined their data so we will combine for
            % them and create the dummy variable. 
            if iscell(contrast) && length(contrast) == 2
                condition1 = S.outcome.(test).contrast{1};
                condition2 = S.outcome.(test).contrast{2};
                m1 = S.brain_data.(condition1).data;
                m2 = S.brain_data.(condition2).data;
                m = cat(2, m1, m2); %TODO: build in a check for dimensions
                cond1_ids = S.brain_data.(condition1).sub_ids;
                cond2_ids = S.brain_data.(condition2).sub_ids;
                both_cond_ids = cat(1, cond1_ids, cond2_ids);

                % save contrast to results
                results.study_info.test_components = {condition1, condition2};


                % make sure the dims of m are subs x parcels
                if size(m,2) == length(S.brain_data.(condition1).sub_ids) + length(S.brain_data.(condition2).sub_ids)
                    % if the second dimension is subjects,
                    % flip dimensions
                    m = m';
                end

                % creating the dummy variable as 'score'
                score = categorical(cat(1, zeros(size(m1,2), 1), ones(size(m2,2), 1))); %TODO: test
                motion1 = S.brain_data.(condition1).motion;
                motion2 = S.brain_data.(condition2).motion;
                motion = cat(1, motion1, motion2); 

                results_file_prefix = [results_dir, res_prefix, S.study_info.dataset, '_', S.study_info.map, '_', test_type, '_', condition1, '_', condition2];

            elseif isnan(contrast) && size(S.outcome.(test).score_label,1) == 1
                condition = S.outcome.(test).reference_condition;
                m = S.brain_data.(condition).data;

                score_label = S.outcome.(test).score_label;

                % save contrast to results
                results.study_info.test_components = {condition, score_label};

                % make sure the dims of m are subs x parcels
                if size(m,2) == length(S.brain_data.(condition).sub_ids)
                    % if the second dimension is subjects,
                    % flip dimensions
                    m = m';
                end

                score = S.outcome.(test).score;
                motion = S.brain_data.(condition).motion;

                % remove missing subject data
                [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);

                results_file_prefix = [results_dir, res_prefix, S.study_info.dataset, '_', S.study_info.map, '_', test_type, '_', condition, '_', S.outcome.(test).score_label];
            end

        end


        for do_pooling = pooling_params

            %% Do large-scale pooling if specified

            if do_pooling
                m2 = []; 
                triumask=logical(triu(ones(n_network_groups)));  
                for s = 1:size(m,1) % over subjects
                    tr = structure_data(m(s,:),'triangleside','upper');
                    t2 = summarize_matrix_by_atlas(tr,'suppressimg',1)'; % transpose because it does tril by default
                    m2(s,:) = t2(triumask);
                end

                % update results file prefix
                % results_file_prefix2 = [results_file_prefix,'__by_net'];

                % set pooling_method for naming
                pooling_method = 'net';

            else
                % results_file_prefix2 = results_file_prefix;
                m2 = m;

                % set pooling_method for naming
                pooling_method = 'none';
            end

            disp(['   > pooling = ', pooling_method])


            for motion_method_it = 1:length(motion_method_params)

                %% Account/correct for motion as specified

                motion_method = motion_method_params{motion_method_it};
                disp(['     > motion method = ', motion_method])
                % set results name
                result_name = ['pooling_', pooling_method, '_motion_', motion_method];

                score2 = score; % changes if applying motion threshold

                if ~strcmp(motion_method,'none')

        %            % Align motion and data by subject - REMOVED 051624
        %            
        %            % align motion subIDs to data subIDs (removes subjects who do not have data or scoreavior)
        %            subids_motion = load(motion_subids_file);
        %            [subids,idx_matrices,idx_motion] = intersect(subids_data,subids_motion,'stable');
        %            if ~isequal(subids_data,subids_motion(idx_motion)); error('Matrix and motion intersected subject IDs don''t match. Likely subjects with matrices are missing motion.'); end;
        %            
        %            % load motion and index based on above alignment
        %            sub_motion = load(motion_file);
        %            sub_motion = sub_motion(idx_motion); % get subs with matrices
        %
        %            std_sub_motion = std(sub_motion); % TODO: save this or full motion vector
        %            mean_sub_motion = mean(sub_motion); % TODO: save this or full motion vector

                    % threshold if specified
                    if strcmp(motion_method,'threshold')
                        low_motion_idx = find(motion<low_motion_threshold); % TODO: consider saving
                        m2 = m2(low_motion_idx,:);
                        score2 = score2(low_motion_idx);
                        std_sub_motion = std(motion(low_motion_idx)); % TODO: save this or trimmed motion vector
                        mean_sub_motion = mean(motion(low_motion_idx)); % TODO: save this or trimmed motion vector

                    end

                    % update results file prefix
                    % results_file_prefix3=[results_file_prefix2,'__with_motion_',motion_method];

                end


                %% Run regression

                if strcmp(motion_method,'regression')
                    % include motion as a confound
                    [b_standardized,p,n,std_brain,std_score] = save_univariate_regression_results(test_type,m2,score2,motion);
                else
                    [b_standardized,p,n,std_brain,std_score] = save_univariate_regression_results(test_type,m2,score2);
                end

                results.data.(result_name).b_standardized = b_standardized;
                results.data.(result_name).p = p;
                results.data.(result_name).n = n;
                results.data.(result_name).std_brain = std_brain;
                results.data.(result_name).std_score = std_score;


            end
        end
        % save this contrast's results file as .mat
        save([results_file_prefix,'.mat'], 'results');

    end
end


% output path - USER-DEFINED

%results_file_prefix = [results_dir,'ukb_fc_t_rest_age__v2']; % see naming convention
%results_file_prefix = [results_dir,strrep(['ukb_fc_test_rest_age__v2'], ' ', '_')']; % see naming convention
% results_file_prefix = [results_dir, 'test_fc_r_rest']; 
% note: moved this into the for loop since the name will depend!

%% Calculate effects
% moved into for loop


%% Function definitions

% set up function for univariate y~x with optional confound as predictor

% NOTE: it may be tempting to use a stepwise procedure - estimate the residuals of m ~ motion, then use those residuals for subsequently estimating betas - but residualizing in this way biases effect size estimates towards zero (i.e., conservative), particularly with higher collinearity in predictors. Instead, including the confound in the actual model should be preferred as unbiased (though also note higher variance with collinearity) - https://besjournals.onlinelibrary.wiley.com/doi/10.1046/j.1365-2656.2002.00618.x
%   i.e., avoid: mdl = fitlm(brain, confound); brain = mdl.Residuals;


function [b_standardized,p,n,std_brain,std_score] = save_univariate_regression_results(test_type,brain,score,confounds)
    % brain: n_sub x n_var, score: n_sub x 1, Optional confounds: n_sub x n_var
    % brain is brain data, score is score, confounds is motion
   
    n = size(brain,1);
    std_brain = std(brain);
    std_score = std(score);
    
    if nargin==3
        confounds=[];
    elseif ~nargin==4
        error('%d arguments provided but only 3 or 4 allowed.',nargin)
    end
   
    % adjust design and which coefficient results to extract based on test type
    if strcmp(test_type, 'r')
        do_t_test=0;
        test_beta=2;
    elseif strcmp(test_type, 't2')
        do_t_test=1;
        test_beta=2;
    elseif strcmp(test_type, 't')
        do_t_test=1;
        test_beta=1;
        score=[]; % because score is just single group ID, but we already have this in the intercept
        std_score=1; % for one-sample t-test b->r conversion, in order to not affect the result
    end
        
    % Run regression across each brain variable

    if do_t_test
        
        % design matrix representing group ID or intercept is predictor - standard design for t-test 
        for c=1:size(brain,2)
            has_brain_data=~isnan(brain(:,c));
            n_has_brain_data = sum(has_brain_data);
             if ~isempty(score)
                 score2=score(has_brain_data);
%                 mdl(:,:,c) = Regression_fast([ones(1,n),score2,confounds], brain(has_brain_data,c), 1); % note: fitlm is built-in for this but too slow % TODO: check p-value calculation - some set to 0,  maybe singular for edge-wise
            else
                score2 = score;
            end
            if ~isempty(confounds)
                confounds2 = confounds(has_brain_data);
            else
                confounds2 = confounds;
            end
            mdl(:,:,c) = Regression_fast([ones(n_has_brain_data,1), score2, confounds2], brain(has_brain_data,c), 1); % note: fitlm is built-in for this but too slow % TODO: check p-value calculation - some set to 0,  maybe singular for edge-wise
            %end
        end
    
    else

        % score is outcome - standard design for correlation
        for c=1:size(brain,2)
            has_brain_data=~isnan(brain(:,c));
            if ~isempty(score)
                score2=score(has_brain_data);
            end
            mdl(:,:,c) = Regression_fast([ones(n,1),brain(has_brain_data,c),confounds], score2, 1); % note: fitlm is built-in for this but too slow % TODO: check p-value calculation - some set to 0,  maybe singular for edge-wise
            % same results without covariate: [r, p] = corr(brain,score);
        end
   
    end
    
    b = squeeze(mdl(test_beta,1,:))'; % unstandardized betas
    p = squeeze(mdl(test_beta,2,:))';
    if do_t_test
        b_standardized = b .* std_score ./ std_brain; % standardized betas - https://www3.nd.edu/~rwilliam/stats1/x92.pdf - TODO: this is not technically a "partial r" - decide what to do for subsequent R^2 or Cohen's d
    else
        b_standardized = b .* std_brain ./ std_score; % standardized betas - https://www3.nd.edu/~rwilliam/stats1/x92.pdf - TODO: this is not technically a "partial r" - decide what to do for subsequent R^2 or Cohen's d
    end
            
    % save([results_file_prefix,'.mat'],'b_standardized','p','std_brain','std_score','n')
end

