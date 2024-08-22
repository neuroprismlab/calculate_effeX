%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effect maps created by Rongtao + Stephanie
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
%       - brain_mask
%       - category
%   - data
%       - <pooling strategy>
%            - <motion strategy>
%                 - r
%                 - p
%                 - std_X
%                 - std_y
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
data_filename = '/work/neuroprism/effect_size/data/subject_level/s_abcd_fc_rosenblatt.mat';


% motion paths - USER-DEFINED - REMOVED 051624

%motion_file = [motion_dir,'EXAMPLEPATH_allsubs_Movement_RelativeRMS_mean.txt'];
%motion_subids_file = [motion_dir,'subids']; % subject IDs for motion measurements

% script paths

%data_loader_script_path = [scripts_dir,'data_loader']; % USER-DEFINED [this is only a placeholder and does not exist] - REMOVED 051624
atlas_tools_script_path = [scripts_dir,'atlas_tools']; % for pooling (structure_data, summarize_data)
regression_fast_script_path = [scripts_dir,'regression_fast']; % for fast multiple regression

% params

% pooling_params = [0,1]; %TODO: for now set pooling_params to 0 because
% UKB data
pooling_params = [0];
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

S = load(data_filename);

% run checker on data
S = checker(S);


% create a struct to store the results 
% results struct will have one struct for each contrast
% each contrast struct will have one struct for each combination of pooling
% and motion regression

tests = fieldnames(S.outcome);

for i = 1:length(tests)
    % for each outcome...
    test = tests{i}; % get the name of the outcome/score
    disp(['running: ', test])
    
    results = []; % clear the struct from prior loops
    results.study_info.dataset = S.study_info.dataset;
    results.study_info.test = S.study_info.test;
    results.study_info.map = S.study_info.map;
    results.study_info.mask = S.study_info.mask; % TODO: make sure all input is mask not brain_mask
    if isfield(S.outcome.(test), 'level_map')
        results.study_info.level_map = S.outcome.(test).level_map;
    end
        
    
    test_type = infer_test_type(S, test);
    
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
        results_file_prefix = [results_dir, 'hstest_', S.study_info.dataset, '_', S.study_info.map, '_', S.study_info.test, '_', score_label, '_', condition];
        % TODO: remove 'hstest' after testing done
        
        % remove missing subject data
        [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);
       
    end
    
    %TODO: 
    if test_type == 't'
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
        
        results_file_prefix = [results_dir, 'hstest_', S.study_info.dataset, '_', S.study_info.map, '_', S.study_info.test, '_', condition1, '_', condition2];
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
            
            results_file_prefix = [results_dir, 'hstest_', S.study_info.dataset, '_', S.study_info.map, '_', test_type, '_', condition1, '_', condition2];

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
            
            results_file_prefix = [results_dir, 'hstest_', S.study_info.dataset, '_', S.study_info.map, '_', test_type, '_', condition, '_', S.outcome.(test).score_label];
        end

    end
    

    for do_pooling = pooling_params

        %% Do large-scale pooling if specified

        if do_pooling

            m2 = []; 
            triumask=logical(triu(ones(n_network_groups)));  
            for i = 1:size(m,1) % over subjects
                t = structure_data(m(i,:),'triangleside','upper');
                t2 = summarize_matrix_by_atlas(t,'suppressimg',1)'; % transpose because it does tril by default
                m2(i,:) = t2(triumask);
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


        for motion_method_it = 1:length(motion_method_params)

            %% Account/correct for motion as specified
            
            motion_method = motion_method_params{motion_method_it};
            disp(['running motion methhod: ', motion_method])
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
                [r,p,n,std_X,std_y] = save_univariate_regression_results(m2,score2, motion);
            else
                [r,p,n,std_X,std_y] = save_univariate_regression_results(m2,score2);
            end
            results.data.(result_name).r = r;
            results.data.(result_name).p = p;
            results.data.(result_name).n = n;
            results.data.(result_name).std_X = std_X;
            results.data.(result_name).std_y = std_y;
            

        end
    end
    % save this contrast's results file as .mat
    save([results_file_prefix,'.mat'], 'results');
    
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
%   i.e., avoid: mdl = fitlm(X, confound); X = mdl.Residuals;

% TODO: either convert the above to point-biserial correlation coefficient or t-statistic or add a function to calculate these

% separate function first that 

function [r,p,n,std_X,std_y] = save_univariate_regression_results(X,y,X2)
    % X: n_sub x n_var, y: n_sub x 1, Optional X2: n_sub x n_var
    % X is brain data, y is score, X2 is motion
    n = length(y);
    std_X = std(X);
    std_y = std(y);
    
    if nargin==3 % regress score and motion
            
            % for each brain variable, perform regression
            for i=1:size(X,2)
                mdl(:,:,i) = Regression_fast([ones(n,1),X(:,i),X2], y, 1); % note: fitlm is built-in for this but too slow % TODO: check p-value calculation - some set to 0,  maybe singular for edge-wise
            end
            
            b = squeeze(mdl(2,1,:))'; % unstandardized betas
            r = b.*std_X / std_y; % standardized betas - https://www3.nd.edu/~rwilliam/stats1/x92.pdf - TODO: this is not technically a "partial r" - decide what to do for subsequent R^2 or Cohen's d
            p = squeeze(mdl(2,2,:))';
       
    elseif nargin==2 % regress brain data with score
            
            [r, p] = corr(X,y);
            
    else
        error('%d arguments provided but only 2 or 3 allowed.',nargin)
        
    end
    
    % TODO: the dimensions are flipped between the output from no motion
    % and from motion regression the way this is now
    
    % save([results_file_prefix,'.mat'],'r','p','std_X','std_y','n')
end

