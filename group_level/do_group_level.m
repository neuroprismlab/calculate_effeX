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
%       - mask
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



%% Overarching setup

% script and data directories

results_dir = '/work/neuroprism/effect_size/data/group_level/'; % USER-DEFINED
data_dir = '/work/neuroprism/effect_size/data/subject_level/'; % USER-DEFINED

current_file = mfilename('fullpath');
scripts_dir = [fileparts(current_file),'/helper_scripts/'];

% params

motion_method_params = {'none', 'regression', 'threshold'};
low_motion_threshold = 0.1; % empirically, 5.7% of subjects are >0.1mm mFFD
n_network_groups = 10; % hard-coded for Shen atlas-based pooling
pooling_params = [0, 1];
testing=1; % USER-DEFINED

% get list of input data filenames
filenames = dir(data_dir);
filenames = filenames(~[filenames.isdir]);  % Remove directories from the list
datasets = {filenames.name}; 

% set paths
addpath(genpath(scripts_dir));

% setup for tests

if testing
    datasets = {'s_hcp_fc_noble_corr.mat'};
    testing_str ='test';
else
    testing_str = [];
end


%% Calculate effects for each test of each dataset

disp(upper(testing_str))

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
    
    % TODO: TMP: for testing one specific test, figure out and remove
    if dataset == "s_pnc_fc_ye.mat"
        tests = tests(1:10);
    end
    
    for t = 1:length(tests)
        
        test = tests{t}; 
        disp(['Running test "', test,'"'])
        
        % Create new results struct per test (contains all combinations of pooling + motion method)
        % TODO: later this is "results.data.(result_name)" - resolve this difference
        results = [];
        results.study_info.dataset = S.study_info.dataset;
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
        
        results.study_info.test = test_type;
        

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
                score = cat(1, zeros(size(m1,2), 1), ones(size(m2,2), 1)); %TODO: test
                % TODO: consider the simpler:
                % score = [zeros(size(m1,2)); ones(size(m2,2))];

                % get test components and add to results
                results.study_info.test_components = {condition1, condition2};
                
%                 % calculate n1 and n2
%                 n1 = sum(S.outcome.(test).score==1);
%                 n2 = sum(S.outcome.(test).score==0);
%                 
%                 % save n1 and n2
%                 results.study_info.n1 = n1;
%                 results.study_info.n2 = n2;

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
                
%                 % calculate n1 and n2
%                 n1 = sum(S.outcome.(test).score==1);
%                 n2 = sum(S.outcome.(test).score==0);
%                 
%                 % save n1 and n2
%                 results.study_info.n1 = n1;
%                 results.study_info.n2 = n2;

            end
        end


        results_file_prefix = [results_dir, strjoin([S.study_info.dataset, S.study_info.map, test_type, results.study_info.test_components, date, testing_str], '_')];



        %% Run analysis for each do_pooling + motion + do_multivariate method

        %for do_multivariate = multivariate_params
        
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

                        % Align motion and data by subject - REMOVED 051624  (removed the commented lines lines that align subjects) - TODO: check + confirm the checker uses similar logic as the removed lines
            
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
                    %if do_multivariate
                    %    test_type=['multi_',test_type];
                    %    % TODO: consider whether safe to overwrite
                    %end

                    if strcmp(motion_method,'regression')
                        % include motion as a confound
                        [b_standardized,p,n,std_brain,std_score] = run_test(test_type,m2,score2,motion);
                    else
                        [b_standardized,p,n,std_brain,std_score] = run_test(test_type,m2,score2);
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
                    
                    

                end % motion
            end % pooling
        %end % multivariate

        % if category provided for the test, save
        if isfield(S.outcome.(test), 'category')
            results.study_info.category = S.outcome.(test).category;
        else
            results.study_info.category = NaN;
        end

        % Save all results for this test
        save([results_file_prefix,'.mat'], 'results');
       
    end % tests
end % datasets







%% Function definitions

% NOTE: Deconfounding via "residualizing" (i.e., fit brain ~ motion, then use brain residuals for subsequently estimating betas - equivalently mdl=fitlm(brain,confound); brain=mdl.Residuals) is known to bias univariate effect size estimates towards zero (i.e., conservative), particularly in the case of higher collinearity between predictors. Including the confound directly in the model should be preferred as unbiased (though there will also be higher variance in estimates with collinearity) - https://besjournals.onlinelibrary.wiley.com/doi/10.1046/j.1365-2656.2002.00618.x . We avoid deconfounding via residualizing for univariate estimates, but this appears to be the standard for multivariate estimates and to our knowledge the extent to which bias may be introduced is unclear.

function [b_standardized,p,n,std_brain,std_score] = run_test(test_type,brain,score,confounds)
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
    switch test_type
        % Run mass univariate tests (for each brain region) or multivariate tests
        
        
        
        case 't'
            % Standard 1-Sample t-Test (Mass Univariate): brain is outcome
            
            score=[]; % because score is just single group ID, but we already have this in the intercept
            std_score=1; % for one-sample t-test b->r conversion, in order to not affect the result
    
            % need intercept so can't use corr
            mdl = Regression_fast_mass_univ_y([ones(n,1), score, confounds], brain, 1); % note: fitlm is built-in for this but too slow for this purpose
            b_standardized = mdl(1,:,1);
            p = mdl(1,:,2);

        
        
        case 't2'
            % Standard 2-Sample t-Testi (Mass Univariate): group ID is predictor, brain is outcome
        
            if isempty(confounds)
                [b_standardized,p]=corr(brain,score);
            else
                [b_standardized,p]=partialcorr(brain,score,confounds);
            end
     
            % TODO: revisit whether addl info needed for subsequent R^2 or d - https://www3.nd.edu/~rwilliam/stats1/x92.pdf 
            % TODO: check p-value calculation - some previously set to 0, maybe singular for edge-wise



        case 'r'
            % Standard Correlation (Mass Univariate): score is outcome (standard correlation)
        
            if isempty(confounds)
                [b_standardized,p]=corr(score,brain);
            else
                [b_standardized,p]=partialcorr(score,brain,confounds);
            end


        %{
        case 'multi_t'
            %  Hotelling t-Test (Multivariate): brain is outcome

            % 1. Dimensionality reduction - slow (~10 sec)
            [~,brain_reduced]=pca(brain,'NumComponents',n_components,'Centered','off'); % make sure not to center - we're measuring dist from 0!

            % 2. Optional: regress confounds from brain
            if isempty(confounds)
                brain2 = brain_reduced;
            else
                confound_centered = confound - mean(confound);
                P_confound = confound_centered / (confound_centered' * confound_centered) * confound_centered'; % confound projection matrix
                brain2 = brain_reduced - P_confound * brain_reduced;
            end

            % 3. Hotelling t-test 
            [t_sq,p] = Hotelling_T2_Test(brain2); 
            b_standardized = sqrt(t_sq); % this may also be the biased unbiased estimate of mahalanobis d
            
            % d=t --> same result as from a direct estimate of d (sample):
            % d = sqrt(mahal(zeros(1,size(brain2,2)),brain2));
            % and same p as from manova1:
            % [~, p, ~] = manova1(brain2, ones(size(brain2, 1), 1));



        case {'multi_t2', 'multi_r'}
            % Canonical Correlation (Multivariate): brain is predictor, score is outcome (equivalent to the opposite for t-test analogue)

            % 1. Dimensionality reduction - slow (~10 sec)
            [~,brain_reduced]=pca(brain,'NumComponents',n_components);

            % 2. Optional: regress confounds from brain and score
            if isempty(confounds)
                brain2=brain_reduced;
                score2=score;
            end
                confound_centered = confound - mean(confound);
                P_confound = confound_centered / (confound_centered' * confound_centered) * confound_centered'; % confound projection matrix
                brain2 = brain_reduced - P_confound * brain_reduced;
                score2 = score - P_confound * score;
            end

            % 3. Canonical Correlation (top component)
            [brain_comp,score_comp,b_standardized,~,~,stats] = canoncorr(brain2,score2);
            p = stats.pChisq; % TODO: compare with look up from Winkler et al table
            %d = 2*r/sqrt(1-r^2); % TODO: confirm
        %}


    end
end
            

