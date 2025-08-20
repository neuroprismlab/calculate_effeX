%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate group-level statistics
% Authors: Hallee Shearer & Stephanie Noble
%
% PURPOSE:
% This script calculates group-level statistical maps from subject-level statistical maps.
% It processes subject-level brain data (FC or activation) and computes stats for various
% statistical tests with different combinations of:
%     - Dimensionality: univariate or multivariate
%     - Spatial Pooling: voxel/parcel level or network-level
%     - Motion correction: none, regression, or thresholding of mean FD per subject
%
% OVERVIEW:
% 1. Setup paths and user inputs
% 2. For each dataset and each statistical test available to run within that dataset:
%     - extract and format brain data based on test type to be conducted
%     - run all combinations of analysis parameters (pooling x motion x dimensionality)
%     - compute simultaneous confidence intervals
% 3. Save results
% 
% SUPPORTED TEST TYPES:
% - 'r' : correlation between brain data and measures
% - 't' : one-sample t-test
% - 't2' : two-sample t-test between groups
% - 'multi_*': multivariate versions of each test using canonical correlation or Hotelling's T^2
%
% KEY VARIABLES:
% - m/m2: brain data matrices (subjects x variables)
% - score: behavioral/outcome measures (e.g., age, test score, ...)
% - motion: mean FD per subject for confound control
% - results: output structure containing all analysis combinations
%
%
% INPUT FORMAT: 
% Each input filename within the input directory contains a data structure of the following form:
%
%   - study_info
%       - dataset
%       - map
%       - test
%       - [mask] ---------------- if mask is the same across all conditions; otherwise set below
%       - date
%   - brain_data
%       - <name of condition 1>
%           - sub_ids
%           - data -------------- last dim is n_sub long
%           - sub_ids_motion
%           - motion
%           - [mask] ------------ if mask differs by condition; otherwise set above
%           - [mask_hdr] -------- for nifti files
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
%           - [level_map] ------- number-to-string map for binary or categorical score (including ref level)
%       - test2
%           ...
%
% NOTES:
%   Motion note: For datasets processed with our Yale pipeline (ABCD, HBN, PNC, SLIM), motion is already regressed from timeseries data (regresses 24 motion parameters)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT: One file per study in output directory
% Each output filename within the output directory will contain a data structure of the following form:
%
%   - study_info
%       - dataset
%       - map
%       - category
%       - mask
%       - [mask_hdr]
%       - test
%       - test_components ------- (e.g., {'condition_label', 'score_label'}
%       - [level_map]
%   - data
%       - <pooling strategy>
%            - <motion strategy>
%                 - stat
%                 - p
%                 - std_brain
%                 - std_score
%                 - [n]         ------------ NaN if two-sample
%                 - [n1]        ------------ NaN if one-sample
%                 - [n2]        ------------ NaN if one-sample
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OVERARCHING SETUP

% DIRECTORY CONFIGURATION
% set input and output directories for data processing
results_dir = '/work/neuroprism/effect_size/data/group_level/'; % USER-DEFINED
data_dir = '/work/neuroprism/effect_size/data/subject_level/'; % USER-DEFINED

% PARAMETER CONFIGURATION
% user can change these according to which parameters to run
motion_method_params = {'none', 'regression', 'threshold'}; 
low_motion_threshold = 0.1; % mean Framewise Displacement threshold (in mm)
% set to other threshold if wanted, then name the motion method accordingly, e.g., 'threshold2' for 0.2mm

% analysis choices to run
n_network_groups = 10; % hard-coded for Shen atlas-based pooling
pooling_params = [0, 1]; % 0 = no pooling, 1 = network-level pooling
multivariate_params = [0, 1]; % 0 = univariate tests, 1 = multivariate tests

testing=0; % USER-DEFINED - change to 1 if testing
if testing
    datasets = datasets(1); % run only certain dataset(s) while testing
    testing_str = '_test'; % add a string to all filenames when testing
else
    testing_str = [''];
end

% ----------------------------------------------------------------------------

% flags
save_info.save = 1; % save results
save_info.overwrite = 1; % allow overwriting existing files (initialized)
save_info.use_same = 0; % apply same overwrite decision to all datasets
save_info.asked = 0; % track if user has been prompted for overwrite decisiion

% get current script directory and set helper paths
[current_dir,~,~] = fileparts(mfilename('fullpath'));
scripts_dir = [current_dir,'/helper_scripts/'];
reference_dir = [current_dir, '/reference/'];

% get list of input data filenames
filenames = dir(data_dir);
filenames = filenames(~[filenames.isdir]);  % Remove directories from the list
datasets = {filenames.name}; 

% set paths
addpath(genpath(scripts_dir));
addpath(genpath(reference_dir));

%% ================ MAIN PROCESSING LOOP: DATASETS & TESTS ================

disp(upper(testing_str))

for i = 1:length(datasets) % loop through all available datasets
   
    dataset = datasets{i};
    fprintf(['Processing dataset: ',dataset,'\n'])
    
    % ===================== FILE HANDLING AND OVERWRITE LOGIC =====================
    % check if results already exist and prompt decisions to overwrite or not
    data_path = [data_dir, dataset];

    S = load(data_path,'study_info');
    results_file_pre_prefix = [results_dir, strjoin({S.study_info.dataset, S.study_info.map, testing_str}, '_'),'*'];

    % ==================== OVERWRITE CONFIRMATION LOGIC ====================
     if ~isempty(dir(results_file_pre_prefix))
         if ~save_info.asked || ~save_info.use_same
 
             user_response = input(sprintf(['Results for ',results_file_pre_prefix,' already exist. Overwrite? [yes/no]\n> ']),'s');
             if strcmp(user_response,'yes')
                 fprintf(['Replacing results ',results_file_pre_prefix,'.\n']);
                 save_info.overwrite = 1;
             else
                 fprintf(['Keeping existing results ',results_file_pre_prefix,'.\n']);
                 save_info.overwrite = 0;
             end
 
             if ~save_info.asked
                 user_response = input(sprintf(['Repeat for all datasets? [yes/no]\n> ']),'s');
                 if strcmp(user_response,'yes')
                     fprintf('Using this setting for all.\n');
                     save_info.use_same = 1;
                 else
                     fprintf('Okay, will ask each time.\n');
                    save_info.use_same = 0;
                end
                save_info.asked = 1;
             end
         end
         save_info.save = save_info.overwrite;
     else
         save_info.save = 1;
     end
       
    if save_info.save

    
    % ==================== DATA LOADING AND VALIDATION ====================
    % Load subject-level data structure and validate its format

    S = load(data_path);
    S = checker(S);

    % get results for each test


    % ================== ITERATE THROUGH ALL TESTS IN DATASET ================
    % Each dataset may contain multiple statistical tests to perform
    
    tests = fieldnames(S.outcome);

    for t = 1:length(tests)
        
        test = tests{t}; 
        disp(['Running test "', test,'"'])
        
        % ================ INITIALIZE RESULTS STRUCTURE FOR TEST ===============
        % Clear results from previous iteration and set up new structure
        results = [];
        
        % Create new results struct per test (contains all combinations of pooling + motion method)
        % if a level map exists, add to results study info 
        if isfield(S.outcome.(test), 'level_map')
            results.study_info.level_map = S.outcome.(test).level_map;
        end

        % infer test type to conduct based on the data available for this test
        test_type = infer_test_type(S, test);

        % store test type, map type, dataset name, and date in results struct
        results.study_info.test = test_type;
        results.study_info.map = S.study_info.map;
        results.study_info.dataset = S.study_info.dataset;
        results.study_info.date = date;
        
        % save mask if stored in study_info, add to results
        % otherwise, mask will be saved to results by test type below
        if isfield(S.study_info, 'mask')
            results.study_info.mask = S.study_info.mask;
        end
       
        
        % ============== TEST TYPE-SPECIFIC DATA EXTRACTION & PROCESSING ==============
        % Extract and format data based on statistical test type (correlation, t-test, etc.)

        % ------------------- CORRELATION TEST (r) -------------------
        % Extract brain data, behavioral scores, and motion parameters
        if test_type == 'r'
            
            % extract relevant variables
            condition = S.outcome.(test).reference_condition; % reference condition specifies brain data for corr
            m = S.brain_data.(condition).data;
            motion = S.brain_data.(condition).motion;
            score = S.outcome.(test).score;
            score_label = S.outcome.(test).score_label;

            % check brain dims are n_subs x n_parcels, otherwise flip 
            if size(m,2) == length(S.brain_data.(condition).sub_ids)
                m = m';
            end
          
            % remove missing subject data (m, score, motion)
            % TODO: for this and elsewhere, should return and save new n after removing missing subs
            [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);
  
            % get test components and add to results
            results.study_info.test_components = {condition, score_label}; 


        % ------------------- T-TEST (one-sample) -------------------
        elseif test_type == 't'
            
            if length(S.outcome.(test).contrast)==1 % unpaired t-test
                % Single condition one-sample t-test
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
                
                % if brain mask is here, save it
                if isfield(S.brain_data.(condition), 'mask')
                    results.study_info.mask = S.brain_data.(condition).mask;
                    results.study_info.mask_hdr = S.brain_data.(condition).mask_hdr;
                end
                 
 
            else % two-condition paired t-test

                % extract relevant variables
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
                motion=mean([motion1,motion2],2);
                
                % to facilitate paired test: take CONDITION *2* - CONDITION *1* difference
                m = m2-m1;
                
                % get test components and add to results
                results.study_info.test_components = {condition1, condition2};
        
            end

        % ------------------- TWO-SAMPLE T-TEST (t2) -------------------
        % Compare brain activity between two groups
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
                
            elseif isnan(contrast) && size(S.outcome.(test).score_label,1) == 1
                
                % else if contrast is NaN, then data is already combined -> use outcome score as dummy variable
                
                % extract relevant variables
                condition = S.outcome.(test).reference_condition;
                m = S.brain_data.(condition).data;
                motion = S.brain_data.(condition).motion;
                score = S.outcome.(test).score;
                score_label = S.outcome.(test).score_label;

                % map score to {0,1} for regression
                score = +(score == max(score));
                
                % TODO: need level_map for interpretation
                %if isfield(results.study_info ,'level_map')
                %    results.study_info.level_map;
                %end 

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


        % ================= RESULTS FILE MANAGEMENT =================
        % Set up results file paths and load existing results if they exist
        
        results_file_prefix = [results_dir, strjoin([S.study_info.dataset, S.study_info.map, test_type, results.study_info.test_components, testing_str], '_')];

        results_file_path = [results_file_prefix, '.mat'];
        if exist(results_file_path, 'file')
            tmp = load(results_file_path);
            if isfield(tmp, 'results')
                results.data = tmp.results.data;
            end
        end
        if ~isfield(results, 'data')
            results.data = struct();
        end
        

        %% ========== NESTED ANALYSIS LOOPS: MULTIVARIATE x POOLING x MOTION ==========
        % Run analysis for each combination of multivariate/univariate, pooling, and motion correction methods

        % ----------- MULTIVARIATE/UNIVARIATE CONFIGURATION -----------
        % Determine whether to run multivariate or univariate statistical tests
        for do_multivariate = multivariate_params
            if do_multivariate
                test_type=['multi_',test_type]; % prepend 'multi_' for multivariate analyses
                % TODO: consider whether safe to overwrite
                mv_test_type = test_type;
            else
                mv_test_type = 'none'; % run univariate analysis
            end
            disp(['   > statistical test = ', test_type])

            %% ============= ATLAS-BASED POOLING CONFIGURATION =============
            % Apply network-level pooling using atlas parcellations if specified
            for do_pooling = pooling_params

                %% Do large-scale pooling if specified
                if do_pooling

                    % ----------- FUNCTIONAL CONNECTIVITY POOLING -----------
                    % Pool connectivity matrix by networks using upper triangular mask
                    if strcmp(results.study_info.map,'fc')
                        
                        triumask = logical(triu(ones(n_network_groups)));  
                        m2 = []; 
                        for s = 1:size(m,1) % over subjects
                            tr = structure_data(m(s,:),'triangleside','upper');
                            t2 = summarize_matrix_by_atlas(tr,'suppressimg',1)'; % transpose because it does tril by default
                            m2(s,:) = t2(triumask);
                        end
                        
                        % TODO: above triumask should be read from the provided mask to avoid accidents

                    % ----------- TASK ACTIVATION POOLING -----------
                    % Pool activation data using dilated Shen network atlas
                    elseif strcmp(results.study_info.map,'act')

                        % load, dilate, and mask shen network atlas
                        % for now, we're using an appromiximate mapping without fully registering since we're working with large-scale networks
                        shen_nets_filename = [current_dir,'/helper_scripts/atlas_tools/atlas_mappings/shen_1mm_268_parcellation__in_subnetworks.nii.gz'];
                        shen_nets = niftiread(shen_nets_filename);
                        shen_nets = imdilate(shen_nets, strel('cube', 3)); % dilate for approximate registration
                        shen_nets = imresize3(shen_nets, 0.5,'Method','nearest'); % downsample to match data resolution
                        shen_nets = shen_nets(find(results.study_info.mask)); % work with data in the reduced dims of mask, not full 3D
                        % niftiwrite(shen_nets,[results_dir,'dilated_shen_nets.nii']) % run before the previous line
                        % niftiwrite(double(results.study_info.mask),[results_dir,'tmp_msk.nii'])  

                        % average activation within each network for each subject
                        m2 = [];
                        for s = 1:size(m,1) % loop over subjects
                             m2(s,:) = average_within_atlas(m(s,:)', shen_nets, 0);
                        end

                    end

                    % set pooling_method for result naming
                    pooling_method = 'net';

                % ----------- NO POOLING: USE ORIGINAL DATA -----------
                else
                    m2 = m;

                    % set pooling_method for naming
                    pooling_method = 'none';
                end

                disp(['     > pooling = ', pooling_method])

                % preserve copies of variables that could be altered in the loops
                % specifically this is for thresholding since subjects are removed in the loop
                base_m2   = m2;      % brain matrix after pooling choice
                baseScore = score;   % score vector aligned with base_m2
                baseMotion = motion; % motion vector aligned with base_m2

                % ------------ LOOP THROUGH MOTION CORRECTION METHODS ------------
                for motion_method_it = 1:length(motion_method_params)

                    %% Account/correct for motion as specified

                    motion_method = motion_method_params{motion_method_it};
                    disp(['       > motion method = ', motion_method])

                    raw_name = ['pooling_', pooling_method, '_motion_', motion_method, '_mv_', mv_test_type];
                    result_name = matlab.lang.makeValidName(raw_name);
                    
                    % if already present, skip this method
                    if isfield(results.data, result_name)
                        disp(['Skipping (already exists): ', result_name]);
                        continue
                    end
                    
                    % reset working copies for each motion method
                    m2_work   = base_m2;
                    score2    = baseScore;
                    motion2   = baseMotion;   % keep a local to make intent clear

                    % ------------ MOTION CORRECTION --------------
                    if ~strcmp(motion_method,'none')

                        % threshold if specified
                        % remove subjects with mean FD above low_motion_threshold
                        if strcmp(motion_method,'threshold')
                            low_motion_idx = find(motion2<low_motion_threshold); % TODO: consider saving
                            m2_work = m2_work(low_motion_idx,:);
                            if ~isempty(score2)
                                score2 = score2(low_motion_idx);
                            end
                            std_sub_motion = std(motion2(low_motion_idx)); % TODO: save this or trimmed motion vector
                            mean_sub_motion = mean(motion2(low_motion_idx)); % TODO: save this or trimmed motion vector
                            % TODO: consider saving number of subjects who are above motion threshold
                        end
                    end


                    %--------------- MOTION REGRESSION ----------------
                    % if motion method is regression, include motion as a confound regressor
                    if strcmp(motion_method,'regression')
                        % include motion as a confound
                        [stat,p,n,n1,n2,std_brain,std_score, stat_fullres, p_fullres] = run_test(test_type,m2_work,score2,motion2);
                    else % otherwise run test without motion as a regressor
                        [stat,p,n,n1,n2,std_brain,std_score] = run_test(test_type,m2_work,score2);
                    end


                    % ------------ APPEND RESULTS -----------------

                    result_name = ['pooling_', pooling_method, '_motion_', motion_method,'_mv_', mv_test_type];
                    results.data.(result_name).stat = stat;
                    results.data.(result_name).p = p;
                    results.data.(result_name).n = n;
                    results.data.(result_name).n1 = n1;
                    results.data.(result_name).n2 = n2;
                    results.data.(result_name).std_brain = std_brain;
                    results.data.(result_name).std_score = std_score;
                    results.data.(result_name).pooling_method = pooling_method;
                    results.data.(result_name).motion_method = motion_method;
                    results.data.(result_name).mv_method = mv_test_type;
                    
                    if strcmp(motion_method,'regression')
                        results.data.(result_name).stat_fullres = stat_fullres;
                        results.data.(result_name).p_fullres = p_fullres;
                    end
                end % motion
            end % pooling
        end % multivariate

        % if category provided for the test, save
        if isfield(S.outcome.(test), 'category')
            results.study_info.category = S.outcome.(test).category;
        else
            results.study_info.category = NaN;
        end

        % Save all results for this test
        save([results_file_prefix,'.mat'], 'results');
       
    end % tests
    end % overwrite
end % datasets


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: Deconfounding via "residualizing" (i.e., fit brain ~ motion, then use brain residuals for subsequently estimating betas - equivalently mdl=fitlm(brain,confound); brain=mdl.Residuals) is known to bias univariate effect size estimates towards zero (i.e., conservative), particularly in the case of higher collinearity between predictors. Including the confound directly in the model should be preferred as unbiased (though there will also be higher variance in estimates with collinearity) - https://besjournals.onlinelibrary.wiley.com/doi/10.1046/j.1365-2656.2002.00618.x . We avoid deconfounding via residualizing for univariate estimates, but this appears to be the standard for multivariate estimates and to our knowledge the extent to which bias may be introduced is unclear.
% TODO: this is also the standard analogue to partial corr estimates, as validated by our effect_size_ref_validation. Revisit.

% ------------- run_test -------------
% this function performs the statistical test on brain data with optional motion confound regression. 
% handles both univariate and multivariate analyses.

function [stat,p,n,n1,n2,std_brain,std_score, varargout] = run_test(test_type,brain,score,confounds)
    % brain: n_sub x n_var, score: n_sub x 1, Optional confounds: n_sub x n_var
    % brain is brain data, score is score, confounds is motion

    % Check arguments
    if nargin==3
        confounds=[];
    elseif ~nargin==4
        error('%d arguments provided but only 3 or 4 allowed.',nargin)
    end

    %  ------------ DATA CLEANING -----------
    % Remove incomplete cases
    if contains(test_type, 'r') || contains(test_type, 't2')
        if isempty(confounds)
            complete_cases = all(~isnan(brain), 2) & all(~isnan(score), 2);
            brain = brain(complete_cases,:);
            score = score(complete_cases,:);
        else % include confounds
            complete_cases = all(~isnan(brain), 2) & all(~isnan(score), 2) & all(~isnan(confounds), 2);
            brain = brain(complete_cases,:);
            score = score(complete_cases,:);
            confounds = confounds(complete_cases,:); 
        end 
    else % 't' - ignore score
       if isempty(confounds)
            complete_cases = all(~isnan(brain), 2);
            brain = brain(complete_cases,:);
        else % include confounds
            complete_cases = all(~isnan(brain), 2) & all(~isnan(confounds), 2);
            brain = brain(complete_cases,:);
            confounds = confounds(complete_cases,:); 
        end 
    end

    % ---------- BASIC DESCRIPTIVE STATISTICS ---------------
    % calculate sample sizes and standard deviations
    n = size(brain,1); % total sample size
    if contains(test_type,'t2') % Two-sample test group sizes
        n1 = sum(score==0);
        n2 = sum(score==1);
    else
        n1 = NaN; % not applicable for other tests
        n2 = NaN;
    end
    std_brain = std(brain); % standard deviation of brain data across subjects
    % standard deviation of scores for tests that include scores (r and t2)
    if contains(test_type,'t2') || contains(test_type,'r') % score only exists for t2 & r
        std_score = std(score);
    else
        std_score = NaN; 
    end
     
    %  ------------- MULTIVARAITE ANALYSIS SETUP --------------
    % configure dimensionality reduction for multivariate tests
    if contains(test_type,'multi')
        n_components = floor(n/50); % set the number of components
        n_vars = size(brain,2);
        if n_components > n_vars
            n_components = n_vars; % can't have more components than variables
        end
    end
 

    % "stat" is exactly the statistic specified by "test_type" (e.g., t-statistic for stat="t")
    % Note: when dealing with confounds, "stat/p" is an estimate of the statistic in a multiple regression framework
    %       and stat_fullres / p_fullres is an estimate of the statistic in the case of zero confounds

    % -------------  STATISTICAL TEST EXECUTION -----------
    % perform the specified statistical test
    
    switch test_type
        % Run mass univariate tests (for each brain region) or multivariate tests

        % one-sample t test
        case 't'
            % Standard 1-Sample t-Test (Mass Univariate): brain is outcome
            % note re Regression_fast: fitlm is built-in for this but too slow for this purpose; need intercept so can't use corr
           
            [stat,p,B] = Regression_faster_mass_univ_y([ones(n,1), confounds], brain, 1);
            
            if ~isempty(confounds)
                brain_res_plus_intercept = brain - [ones(n,1), confounds] * B + B(1);
                [stat_fullres, p_fullres] = Regression_faster_mass_univ_y(ones(n,1), brain_res_plus_intercept, 1);
                varargout{1} = stat_fullres;
                varargout{2} = p_fullres;
            end

        % two-sample t-test
        case 't2'
            % Standard 2-Sample t-Test (Mass Univariate): group ID is predictor, brain is outcome
            
            if isempty(confounds)
                [stat,p] = Regression_faster_mass_univ_y([ones(n,1), score], brain, 2);
            else
                [~, stat_fullres, p_fullres, ~, stat, p] = partial_and_semipartial_corr(brain, score, confounds, n);
                varargout{1} = stat_fullres;
                varargout{2} = p_fullres;
            end
     
            % TODO: check p-value calculation - some previously set to 0, maybe singular for edge-wise

        % brain-behavior correlation
        case 'r'
            % Standard Correlation (Mass Univariate): score is outcome (standard correlation)
           
            if isempty(confounds)
                [stat,p]=corr(score,brain);
            else
                [stat_fullres, ~, p_fullres, stat, ~, p] = partial_and_semipartial_corr(score, brain, confounds, n);
                varargout{1} = stat_fullres;
                varargout{2} = p_fullres;
            end

        % multivariate one-sample t test
        case 'multi_t'
            %  Hotelling t-Test (Multivariate): brain is outcome
           
            % 1. If confounds: regress confounds from brain
            if isempty(confounds)
                brain2 = brain;
            else
                confound_centered = confounds - mean(confounds);
                P_confound = confound_centered / (confound_centered' * confound_centered) * confound_centered'; % confound projection matrix
                brain2 = brain - P_confound * brain;
            end
 
            % 2. Dimensionality reduction - slow (~10 sec)
            [~,brain_reduced] = pca(brain2, 'NumComponents', n_components, 'Centered', 'off'); % make sure not to center - we're measuring dist from 0! % aiming for 50 samples/feature for stable results a la Helmer et al.
            
            % 3. Hotelling t-test 
            [t_sq,p] = Hotelling_T2_Test(brain_reduced); 
            stat = sqrt(t_sq); % this may also be the biased unbiased estimate of mahalanobis d

            if ~isempty(confounds)
                stat_fullres = stat;
                p_fullres = p;
                stat = [];
                p = [];

                varargout{1} = stat_fullres;
                varargout{2} = p_fullres;
            end

        % multivariate two-sample t-test and brain behavior correlation
        case {'multi_t2', 'multi_r'}
            % Canonical Correlation (Multivariate): brain is predictor, score is outcome (equivalent to the opposite for t-test analogue)

            % 1. If confounds: regress confounds from brain and score
            if isempty(confounds)
                brain2=brain;
                score2=score;
            else
                % confounds = confounds(complete_cases,:); % filter for only complete cases
                confound_centered = confounds - mean(confounds);
                P_confound = confound_centered / (confound_centered' * confound_centered) * confound_centered'; % confound projection matrix
                brain2 = brain - P_confound * brain;
                score2 = score - P_confound * score;
            end
            
            % 2. Dimensionality reduction - slow (~10 sec)
            [~,brain_reduced] = pca(brain2, 'NumComponents', n_components); % aiming for 50 samples/feature for stable results a la Helmer et al.

            
            % 3. Canonical Correlation (top component)
            [brain_comp,score_comp,stat,~,~,stats] = canoncorr(brain_reduced,score2);
            p = stats.pChisq; % TODO: compare with look up from Winkler et al table

            % If confounds: repeat Canonical Correlation with semipartial analogue
            if ~isempty(confounds)
                
                stat_fullres = stat;
                p_fullres = p;
                
                varargout{1} = stat_fullres;
                varargout{2} = p_fullres;

                % relative to total variance in y, without regression (akin to semipartial r)
                switch test_type
                case 'multi_t2'
                    [~,brain_reduced_w_confound] = pca(brain, 'NumComponents', n_components); 
                    % stat_full_y is the correlation obtained without regressing confounds from y
                    [brain_comp,score_comp,stat,~,~,stats] = canoncorr(brain_reduced_w_confound,score2);
                case 'multi_r'
                    [brain_comp,score_comp,stat,~,~,stats] = canoncorr(brain_reduced,score);
                end
                %p = stats.pChisq;
            end

            % If t-test: convert to t-stat
            if strcmp(test_type,'multi_t2')
                if ~isempty(confounds)
                    % use conversion from partial_and_semipartial_corr.m function
                    stat = r_to_test_stats(stat, n, 1, 2);
                    varargout{1} = r_to_test_stats(stat_fullres, n, 1, 2);
                else
                    stat = r_to_test_stats(stat, n, 0, 2);
                end
            end
    end
end
            
