%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate group-level statistics
% Authors: Hallee Shearer & Stephanie Noble
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT: User defines input and output directories.
% Each input filename within the input directory contains a data structure of the following form:
%
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
% NOTES:
%   Motion note: For datasets processed with our pipeline, motion already regressed from timeseries data (regresses 24 motion parameters)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT: One file per study in output directory
% Each output filename within the output directory will contain a data structure of the following form:
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TODO: results_file_prefix naming convention - {dataset name}_{contributor name}_{date}
%       subject level: level1 + {date}
%       group level: study, effect map, level2, or group + {date}
%       combined: combined_effect_{date}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Overarching setup

% script and data directories

scripts_dir = '/home/h.shearer/hallee/calculate_effeX/group_level/helper_scripts/';  % TODO: future - automatically get current folder
results_dir = '/work/neuroprism/effect_size/data/group_level/'; % USER-DEFINED
%data_filename = '/work/neuroprism/data_shared/ukb/ukb_data_steph.mat'; % ADDED 051624 - % TODO: remove when done testing
data_dir = '/work/neuroprism/effect_size/data/subject_level/'; % USER-DEFINED

% script paths
atlas_tools_script_path = [scripts_dir,'atlas_tools']; % for pooling (structure_data, summarize_data)
regression_fast_script_path = [scripts_dir,'regression_fast']; % for fast multiple regression

% params

motion_method_params = {'none', 'regression', 'threshold'};
low_motion_threshold = 0.1; % empirically, 5.7% of subjects are >0.1mm mFFD
n_network_groups = 10; % hard-coded for Shen atlas-based pooling
pooling_params = [0, 1];


% get list of input data filenames
filenames = dir(data_dir);
filenames = filenames(~[filenames.isdir]);  % Remove directories from the list
datasets = {filenames.name}; 

% set paths

%addpath(data_loader_script_path); % USER-DEFINED
addpath(genpath(atlas_tools_script_path));
addpath(regression_fast_script_path);



% testing stuff
testing=1;
if testing
    % run first test of every dataset
    
    datasets= {'s_hcp_fc_noble_corr.mat'};

    % if ukb data, pooling_params = [0] because we don't have a map
    % if activation, also skipping pooling- TODO: reinstate when done testing
    % TODO: TMP: change this once we have a ukb map
%     if contains(datasets, "ukb") ||  contains(datasets, "act")
%         pooling_params = [0];
%     end
    % TODO: for now we still want to not pool ukb and act even when not
    % testing
%this_score='age'; % TODO: remove when done testing

end
res_prefix = date; % appended to the start of each result file


%% Calculate effects for each test of each dataset

% TODO: loop through outcomes, and get data according to each test / outcome name

for i = 1:length(datasets)
   
    dataset = datasets{i};
    fprintf(['Processing dataset: ',dataset,'\n'])
    
    %TODO: TMP: remove this once we get ukb map and figure out act pooling
    if contains(dataset, "ukb") ||  contains(dataset, "act")
        pooling_params = [0];
    else
        pooling_params = [0,1];
    end
    
    % load & check data
    data_path = [data_dir, dataset];
    S = load(data_path);
    S = checker(S);

    % get results for each test

    tests = fieldnames(S.outcome);
    
    % TMP: for testing one specific test
    if testing
        tests = tests(1);
    end
    
    for t = 1:length(tests)
        
        test = tests{t}; 
        disp(['Running test "', test,'"'])
        
        % Create new results struct per test (contains all combinations of pooling + motion method)
        % TODO: later this is "results.data.(result_name)" - resolve this difference
        results = [];
        results.study_info.dataset = S.study_info.dataset;
        results.study_info.test = S.study_info.test;
        results.study_info.map = S.study_info.map;
        if isfield(S.study_info, 'mask')
            results.study_info.mask = S.study_info.mask;
        elseif isfield(S.study_info, 'brain_mask')% TODO: make sure all input is mask not brain_mask
            results.study_info.mask = S.study_info.brain_mask;
        end %TODO: move to checker (to change brain_mask to mask)
        if isfield(S.outcome.(test), 'level_map')
            results.study_info.level_map = S.outcome.(test).level_map;
        end

        % infer test type
        test_type = infer_test_type(S, test);
        
        % TMP: skip the test if the test type is t
%         if strcmp(test_type, "t")
%             disp(['skipping test ', test, ' because test type is t'])
%             continue
%         end


        % For each test type, extract and clean data

        if test_type == 'r'
            
            % extract relevant variables
            condition = S.outcome.(test).reference_condition; % reference condition specifies brain data for corr
            m = S.brain_data.(condition).data;
            motion = S.brain_data.(condition).motion;
            score = S.outcome.(test).score;
            score_label = S.outcome.(test).score_label;

            % check brain dims are n_subs x n_parcels, otherwise flip % TODO: dim checking should be added to checker, and then we just flip here
            if size(m,2) == length(S.brain_data.(condition).sub_ids)
                m = m';
            end
          
            % remove missing subject data (m, score, motion)
            % TODO: for this and elsewhere, should return and save new n after removing missing subs
            [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);
  
            % get test components and add to results
            results.study_info.test_components = {condition, score_label}; % TODO: previously reversed - does this need to be reversed for results_filename ?
            
            % TODO: remove 'hstest' after testing done

           
        %TODO: 
        elseif test_type == 't'
            
            if length(S.outcome.(test).contrast)==1 % unpaired t-test
                % TODO: should also check here if score is all 1's - contrast may not be specified
                
                % extract relevant variables
                % TODO: brain mask for t studies
                condition = S.outcome.(test).contrast{1};
                m = S.brain_data.(condition).data;
                motion = S.brain_data.(condition).motion;
                score_label = S.outcome.(test).score_label;
                
                % check brain dims are n_subs x n_parcels, otherwise flip % TODO: dim checking should be added to checker, and then we just flip here
                if size(m,2) == length(S.brain_data.(condition).sub_ids)
                    m = m';
                end
                
                % remove missing subject data
                score = ones(1,size(m,1)); % temporary
                [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);
                score = [];  % score is set below for unpaired one-sample t test; setting placeholder here - TODO: are these three lines really necessary? discuss
               
                % get test components and add to results
                results.study_info.test_components = {condition};
 
            else % paired t-test

                % extract relevant variables
                % TODO: brain mask for t studies
                % assuming that for t, contrast is 2D string array representing two contrast conditions
                % get brain and motion data for both conditions
                condition1 = S.outcome.(test).contrast{1};
                condition2 = S.outcome.(test).contrast{2};
                m1 = S.brain_data.(condition1).data; 
                m2 = S.brain_data.(condition2).data; 
                motion1 = S.brain_data.(condition1).motion;
                motion2 = S.brain_data.(condition2).motion;
                score = [];  % score is set below for paired one-sample t test; setting placeholder here 
                score_label = S.outcome.(test).score_label;
                
                % remove missing subject data
                % TODO: this will have to be updated for m (from m2-m1) and motion (from mean([motion1,motion2],2), checking subids
                [m1, m2, motion1, motion2] = remove_missing_subs(m1, m2, S, test_type, test, condition1, condition2, motion1, motion2);
                        
                % to facilitate paired test: take mean motion for a single confound covariate
                % TODO: check m1 and m2 are both ordered by the same subids before the following diff
                motion=mean([motion1,motion2],2);
                
                % to facilitate paired test: take CONDITION *2* - CONDITION *1* difference
                % TODO: check m1 and m2 are both ordered by the same subids before the following diff
                m = m2-m1;
                
                % get test components and add to results
                results.study_info.test_components = {condition1, condition2};
        
            end
        
        elseif strcmp(test_type, 't2')


            contrast = S.outcome.(test).contrast;

            if iscell(contrast) && length(contrast) == 2

                % if contrast is length 2, then combine data here and create dummy variable
                
                % extract relevant variables
                condition1 = S.outcome.(test).contrast{1};
                condition2 = S.outcome.(test).contrast{2};
                m1 = S.brain_data.(condition1).data;
                m2 = S.brain_data.(condition2).data;
                motion1 = S.brain_data.(condition1).motion;
                motion2 = S.brain_data.(condition2).motion;
                score_label = S.outcome.(test).score_label;

                % concatenate within brain1/2, motion1/2, condition_ids1/2
                m = cat(2, m1, m2); %TODO: build in a check for dimensions
                motion = cat(1, motion1, motion2); 
                cond1_ids = S.brain_data.(condition1).sub_ids;
                cond2_ids = S.brain_data.(condition2).sub_ids;
                both_cond_ids = cat(1, cond1_ids, cond2_ids);
                % TODO: consider the simpler:
                % m = [m1;m2];
                % motion = [motion1;motion2];
                % both_cond_ids = [cond1_ids;cond2_ids];

                % check brain dims are n_subs x n_parcels, otherwise flip % TODO: dim checking should be added to checker, and then we just flip here
                if size(m,2) == length(S.brain_data.(condition1).sub_ids) + length(S.brain_data.(condition2).sub_ids)
                    m = m';
                end
               
                % TODO: should also run remove_missing_subs to compare m vs motion here

                % creating the dummy variable as 'score'
                score = categorical(cat(1, zeros(size(m1,2), 1), ones(size(m2,2), 1))); %TODO: test
                % TODO: consider the simpler:
                % score = [zeros(size(m1,2)); ones(size(m2,2))];

                % get test components and add to results
                results.study_info.test_components = {condition1, condition2};

            elseif isnan(contrast) && size(S.outcome.(test).score_label,1) == 1
                
                % else if contrast is NaN, then data is already combined -> use outcome score as dummy variable
                
                % extract relevant variables
                condition = S.outcome.(test).reference_condition;
                m = S.brain_data.(condition).data;
                motion = S.brain_data.(condition).motion;
                score = S.outcome.(test).score;
                score_label = S.outcome.(test).score_label;

                % check brain dims are n_subs x n_parcels, otherwise flip % TODO: dim checking should be added to checker, and then we just flip here
                if size(m,2) == length(S.brain_data.(condition).sub_ids)
                    m = m';
                end
                
                % remove missing subject data
                [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);

                % get test components and add to results
                results.study_info.test_components = {condition, score_label};

            end
        end
        
        results_file_prefix = [results_dir, res_prefix, S.study_info.dataset, '_', S.study_info.map, '_', test_type, '_', strjoin(results.study_info.test_components, '_')];




        %% Run analysis for each pooling + motion method

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
                        if ~isempty(score2)
                            score2 = score2(low_motion_idx);
                        end
                        std_sub_motion = std(motion(low_motion_idx)); % TODO: save this or trimmed motion vector
                        mean_sub_motion = mean(motion(low_motion_idx)); % TODO: save this or trimmed motion vector
                        % TODO: consider saving number of subjects who are above motion threshold
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


                %% Append results data

                result_name = ['pooling_', pooling_method, '_motion_', motion_method];
                results.data.(result_name).b_standardized = b_standardized;
                results.data.(result_name).p = p;
                results.data.(result_name).n = n;
                results.data.(result_name).std_brain = std_brain;
                results.data.(result_name).std_score = std_score;
                results.data.(result_name).pooling_method = pooling_method;
                results.data.(result_name).motion_method = motion_method;
                
                

            end
        end
        
        % if category provided for the test, save
        if isfield(S.outcome.(test), 'category')
            results.study_info.category = S.outcome.(test).category;
        else
            results.study_info.category = NaN;
        end

        % Save all results for this test
        save([results_file_prefix,'.mat'], 'results');
       
        % previous results filenames for reference - TODO: delete when done
        %results_file_prefix = [results_dir,'ukb_fc_t_rest_age__v2']; % see naming convention
        %results_file_prefix = [results_dir,strrep(['ukb_fc_test_rest_age__v2'], ' ', '_')']; % see naming convention
        % results_file_prefix = [results_dir, 'test_fc_r_rest']; 


    end
end







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
    % TODO: clean up. There's a little more logic here than I think is needed
    if strcmp(test_type, 't2')
        do_t_test=1;
    elseif strcmp(test_type, 't')
        do_t_test=1;
        score=[]; % because score is just single group ID, but we already have this in the intercept
        std_score=1; % for one-sample t-test b->r conversion, in order to not affect the result
    else
        do_t_test=0;
    end
        
    % Run regression across each brain variable

    if do_t_test
        % design matrix representing group ID or intercept is predictor - standard design for t-test 
        
       if strcmp(test_type, 't') % can't use corr, need intercept

            mdl = Regression_fast_mass_univ_y([ones(n,1), score, confounds], brain, 1); % note: fitlm is built-in for this but too slow for this purpose
            
            b_standardized = mdl(:,1,1);
            p = mdl(:,1,2);

        else %t2
            if isempty(confounds)
                [b_standardized,p]=corr(brain,score);
            else
                [b_standardized,p]=partialcorr(brain,score,confounds);
            end
     
            % TODO: revisit whether addl info needed for subsequent R^2 or d - https://www3.nd.edu/~rwilliam/stats1/x92.pdf 
            % TODO: check p-value calculation - some previously set to 0, maybe singular for edge-wise 
        end
    
    else

        % score is outcome - standard design for correlation

        if isempty(confounds)
            [b_standardized,p]=corr(score,brain);
        else
            [b_standardized,p]=partialcorr(score,brain,confounds);
        end
   
    end

end
            

