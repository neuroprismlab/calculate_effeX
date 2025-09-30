function calculate_group_stats(results_dir, data_dir, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate group-level statistics
% Authors: Hallee Shearer & Stephanie Noble
%
% USAGE:
%   calculate_group_stats(results_dir, data_dir)
%   calculate_group_stats(results_dir, data_dir, 'MotionMethods', {'none', 'regression'})
%   calculate_group_stats(results_dir, data_dir, 'LowMotionThreshold', 0.2)
%
% INPUTS:
%   results_dir - Path to output directory for results
%   data_dir    - Path to input directory containing subject-level data
%
% OPTIONAL NAME-VALUE PAIRS:
%   'MotionMethods'       - Cell array of motion correction methods
%                          Default: {'none', 'regression', 'threshold'}
%   'LowMotionThreshold'  - Mean FD threshold in mm for motion thresholding
%                          Default: 0.1
%   'Testing'            - Flag for testing mode (0 or 1)
%                          Default: 0
%
% PURPOSE:
% This function calculates group-level statistical maps from subject-level 
% statistical maps. It processes subject-level brain data (FC or activation) 
% and computes stats for various statistical tests with different combinations of:
%     - Dimensionality: univariate or multivariate
%     - Spatial Pooling: voxel/parcel level or network-level
%     - Motion correction: none, regression, or thresholding of mean FD per subject
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARSE INPUT ARGUMENTS
p = inputParser;
addRequired(p, 'results_dir', @ischar);
addRequired(p, 'data_dir', @ischar);
addParameter(p, 'MotionMethods', {'none', 'regression', 'threshold'}, @iscell);
addParameter(p, 'LowMotionThreshold', 0.1, @isnumeric);
addParameter(p, 'Testing', 0, @(x) isnumeric(x) && (x==0 || x==1));

parse(p, results_dir, data_dir, varargin{:});

motion_method_params = p.Results.MotionMethods;
low_motion_threshold = p.Results.LowMotionThreshold;
testing = p.Results.Testing;

%% PARAMETER CONFIGURATION
n_network_groups = 10; % hard-coded for Shen atlas-based pooling
pooling_params = [0, 1]; % 0 = no pooling, 1 = network-level pooling
multivariate_params = [0, 1]; % 0 = univariate tests, 1 = multivariate tests

if testing
    testing_str = '_test';
else
    testing_str = '';
end

%% SETUP PATHS
[current_dir,~,~] = fileparts(mfilename('fullpath'));
scripts_dir = [current_dir,'/helper_scripts/'];
reference_dir = [current_dir, '/reference/'];

% get list of input data filenames
filenames = dir(data_dir);
filenames = filenames(~[filenames.isdir]);
datasets = {filenames.name};

if testing && ~isempty(datasets)
    datasets = datasets(1);
end

% set paths
addpath(genpath(scripts_dir));
addpath(genpath(reference_dir));

%% INITIALIZE FLAGS
save_info.save = 1;
save_info.overwrite = 1;
save_info.use_same = 0;
save_info.asked = 0;

%% MAIN PROCESSING LOOP
disp(upper(testing_str))

for i = 1:length(datasets)
   
    dataset = datasets{i};
    fprintf(['Processing dataset: ',dataset,'\n'])
    
    % FILE HANDLING AND OVERWRITE LOGIC
    data_path = [data_dir, dataset];

    S = load(data_path,'study_info');
    results_file_pre_prefix = [results_dir, strjoin({S.study_info.dataset, S.study_info.map, testing_str}, '_'),'*'];

    % OVERWRITE CONFIRMATION LOGIC
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
    
    % DATA LOADING AND VALIDATION
    S = load(data_path);
    S = checker(S);

    % ITERATE THROUGH ALL TESTS IN DATASET
    tests = fieldnames(S.outcome);

    for t = 1:length(tests)
        
        test = tests{t}; 
        disp(['Running test "', test,'"'])
        
        % INITIALIZE RESULTS STRUCTURE FOR TEST
        results = [];
        
        if isfield(S.outcome.(test), 'level_map')
            results.study_info.level_map = S.outcome.(test).level_map;
        end

        test_type = infer_test_type(S, test);

        results.study_info.test = test_type;
        results.study_info.map = S.study_info.map;
        results.study_info.dataset = S.study_info.dataset;
        results.study_info.date = date;
        
        if isfield(S.study_info, 'mask')
            results.study_info.mask = S.study_info.mask;
        end
       
        
        % TEST TYPE-SPECIFIC DATA EXTRACTION & PROCESSING
        if test_type == 'r'
            
            condition = S.outcome.(test).reference_condition;
            m = S.brain_data.(condition).data;
            motion = S.brain_data.(condition).motion;
            score = S.outcome.(test).score;
            score_label = S.outcome.(test).score_label;

            if size(m,2) == length(S.brain_data.(condition).sub_ids)
                m = m';
            end
          
            [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);
  
            results.study_info.test_components = {condition, score_label}; 

        elseif test_type == 't'
            
            if length(S.outcome.(test).contrast)==1
                
                condition = S.outcome.(test).contrast{1};
                m = S.brain_data.(condition).data;
                motion = S.brain_data.(condition).motion;
                score_label = S.outcome.(test).score_label;
                
                if size(m,2) == length(S.brain_data.(condition).sub_ids)
                    m = m';
                end
                
                score = ones(1,size(m,1));
                [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);
                score = [];
               
                results.study_info.test_components = {condition};
                
                if isfield(S.brain_data.(condition), 'mask')
                    results.study_info.mask = S.brain_data.(condition).mask;
                    results.study_info.mask_hdr = S.brain_data.(condition).mask_hdr;
                end
                 
            else
                
                condition1 = S.outcome.(test).contrast{1};
                condition2 = S.outcome.(test).contrast{2};
                m1 = S.brain_data.(condition1).data; 
                m2 = S.brain_data.(condition2).data; 
                motion1 = S.brain_data.(condition1).motion;
                motion2 = S.brain_data.(condition2).motion;
                score = [];
                score_label = S.outcome.(test).score_label;
                
                [m1, m2, motion1, motion2] = remove_missing_subs(m1, m2, S, test_type, test, condition1, condition2, motion1, motion2);
                        
                motion=mean([motion1,motion2],2);
                m = m2-m1;
                
                results.study_info.test_components = {condition1, condition2};
        
            end

        elseif strcmp(test_type, 't2')

            contrast = S.outcome.(test).contrast;

            if iscell(contrast) && length(contrast) == 2
                
                condition1 = S.outcome.(test).contrast{1};
                condition2 = S.outcome.(test).contrast{2};
                m1 = S.brain_data.(condition1).data;
                m2 = S.brain_data.(condition2).data;
                motion1 = S.brain_data.(condition1).motion;
                motion2 = S.brain_data.(condition2).motion;
                score_label = S.outcome.(test).score_label;

                m = cat(2, m1, m2);
                motion = cat(1, motion1, motion2); 
                cond1_ids = S.brain_data.(condition1).sub_ids;
                cond2_ids = S.brain_data.(condition2).sub_ids;
                both_cond_ids = cat(1, cond1_ids, cond2_ids);

                if size(m,2) == length(S.brain_data.(condition1).sub_ids) + length(S.brain_data.(condition2).sub_ids)
                    m = m';
                end

                score = cat(1, zeros(size(m1,2), 1), ones(size(m2,2), 1));

                results.study_info.test_components = {condition1, condition2};
                
            elseif isnan(contrast) && size(S.outcome.(test).score_label,1) == 1
                
                condition = S.outcome.(test).reference_condition;
                m = S.brain_data.(condition).data;
                motion = S.brain_data.(condition).motion;
                score = S.outcome.(test).score;
                score_label = S.outcome.(test).score_label;

                score = +(score == max(score));

                if size(m,2) == length(S.brain_data.(condition).sub_ids)
                    m = m';
                end
                
                [m, score, motion] = remove_missing_subs(m, score, S, test_type, test, condition, motion);

                results.study_info.test_components = {condition, score_label};
                
            end
        end

        % RESULTS FILE MANAGEMENT
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
        

        % NESTED ANALYSIS LOOPS
        for do_multivariate = multivariate_params
            if do_multivariate
                test_type=['multi_',test_type];
                mv_test_type = test_type;
            else
                mv_test_type = 'none';
            end
            disp(['   > statistical test = ', test_type])

            for do_pooling = pooling_params

                if do_pooling

                    if strcmp(results.study_info.map,'fc')
                        
                        triumask = logical(triu(ones(n_network_groups)));  
                        m2 = []; 
                        for s = 1:size(m,1)
                            tr = structure_data(m(s,:),'triangleside','upper');
                            t2 = summarize_matrix_by_atlas(tr,'suppressimg',1)';
                            m2(s,:) = t2(triumask);
                        end

                    elseif strcmp(results.study_info.map,'act')

                        shen_nets_filename = [current_dir,'/helper_scripts/atlas_tools/atlas_mappings/shen_1mm_268_parcellation__in_subnetworks.nii.gz'];
                        shen_nets = niftiread(shen_nets_filename);
                        shen_nets = imdilate(shen_nets, strel('cube', 3));
                        shen_nets = imresize3(shen_nets, 0.5,'Method','nearest');
                        shen_nets = shen_nets(find(results.study_info.mask));

                        m2 = [];
                        for s = 1:size(m,1)
                             m2(s,:) = average_within_atlas(m(s,:)', shen_nets, 0);
                        end

                    end

                    pooling_method = 'net';

                else
                    m2 = m;
                    pooling_method = 'none';
                end

                disp(['     > pooling = ', pooling_method])

                base_m2   = m2;
                baseScore = score;
                baseMotion = motion;

                for motion_method_it = 1:length(motion_method_params)

                    motion_method = motion_method_params{motion_method_it};
                    disp(['       > motion method = ', motion_method])

                    raw_name = ['pooling_', pooling_method, '_motion_', motion_method, '_mv_', mv_test_type];
                    result_name = matlab.lang.makeValidName(raw_name);
                    
                    if isfield(results.data, result_name)
                        disp(['Skipping (already exists): ', result_name]);
                        continue
                    end
                    
                    m2_work   = base_m2;
                    score2    = baseScore;
                    motion2   = baseMotion;

                    if ~strcmp(motion_method,'none')

                        if strcmp(motion_method,'threshold')
                            low_motion_idx = find(motion2<low_motion_threshold);
                            m2_work = m2_work(low_motion_idx,:);
                            if ~isempty(score2)
                                score2 = score2(low_motion_idx);
                            end
                            std_sub_motion = std(motion2(low_motion_idx));
                            mean_sub_motion = mean(motion2(low_motion_idx));
                        end
                    end

                    if strcmp(motion_method,'regression')
                        [stat,p,n,n1,n2,std_brain,std_score, stat_fullres, p_fullres] = run_test(test_type,m2_work,score2,motion2);
                    else
                        [stat,p,n,n1,n2,std_brain,std_score] = run_test(test_type,m2_work,score2);
                    end

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
                end
            end
        end

        if isfield(S.outcome.(test), 'category')
            results.study_info.category = S.outcome.(test).category;
        else
            results.study_info.category = NaN;
        end

        save([results_file_prefix,'.mat'], 'results');
       
    end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stat,p,n,n1,n2,std_brain,std_score, varargout] = run_test(test_type,brain,score,confounds)
    
    if nargin==3
        confounds=[];
    elseif ~nargin==4
        error('%d arguments provided but only 3 or 4 allowed.',nargin)
    end

    if contains(test_type, 'r') || contains(test_type, 't2')
        if isempty(confounds)
            complete_cases = all(~isnan(brain), 2) & all(~isnan(score), 2);
            brain = brain(complete_cases,:);
            score = score(complete_cases,:);
        else
            complete_cases = all(~isnan(brain), 2) & all(~isnan(score), 2) & all(~isnan(confounds), 2);
            brain = brain(complete_cases,:);
            score = score(complete_cases,:);
            confounds = confounds(complete_cases,:); 
        end 
    else
       if isempty(confounds)
            complete_cases = all(~isnan(brain), 2);
            brain = brain(complete_cases,:);
        else
            complete_cases = all(~isnan(brain), 2) & all(~isnan(confounds), 2);
            brain = brain(complete_cases,:);
            confounds = confounds(complete_cases,:); 
        end 
    end

    n = size(brain,1);
    if contains(test_type,'t2')
        n1 = sum(score==0);
        n2 = sum(score==1);
    else
        n1 = NaN;
        n2 = NaN;
    end
    std_brain = std(brain);
    if contains(test_type,'t2') || contains(test_type,'r')
        std_score = std(score);
    else
        std_score = NaN; 
    end
     
    if contains(test_type,'multi')
        n_components = floor(n/50);
        n_vars = size(brain,2);
        if n_components > n_vars
            n_components = n_vars;
        end
    end
 
    switch test_type
        case 't'
            [stat,p,B] = Regression_faster_mass_univ_y([ones(n,1), confounds], brain, 1);
            
            if ~isempty(confounds)
                brain_res_plus_intercept = brain - [ones(n,1), confounds] * B + B(1);
                [stat_fullres, p_fullres] = Regression_faster_mass_univ_y(ones(n,1), brain_res_plus_intercept, 1);
                varargout{1} = stat_fullres;
                varargout{2} = p_fullres;
            end

        case 't2'
            if isempty(confounds)
                [stat,p] = Regression_faster_mass_univ_y([ones(n,1), score], brain, 2);
            else
                [~, stat_fullres, p_fullres, ~, stat, p] = partial_and_semipartial_corr(brain, score, confounds, n);
                varargout{1} = stat_fullres;
                varargout{2} = p_fullres;
            end

        case 'r'
            if isempty(confounds)
                [stat,p]=corr(score,brain);
            else
                [stat_fullres, ~, p_fullres, stat, ~, p] = partial_and_semipartial_corr(score, brain, confounds, n);
                varargout{1} = stat_fullres;
                varargout{2} = p_fullres;
            end

        case 'multi_t'
            if isempty(confounds)
                brain2 = brain;
            else
                confound_centered = confounds - mean(confounds);
                P_confound = confound_centered / (confound_centered' * confound_centered) * confound_centered';
                brain2 = brain - P_confound * brain;
            end
 
            [~,brain_reduced] = pca(brain2, 'NumComponents', n_components, 'Centered', 'off');
            
            [t_sq,p] = Hotelling_T2_Test(brain_reduced); 
            stat = sqrt(t_sq);

            if ~isempty(confounds)
                stat_fullres = stat;
                p_fullres = p;
                stat = [];
                p = [];

                varargout{1} = stat_fullres;
                varargout{2} = p_fullres;
            end

        case {'multi_t2', 'multi_r'}
            if isempty(confounds)
                brain2=brain;
                score2=score;
            else
                confound_centered = confounds - mean(confounds);
                P_confound = confound_centered / (confound_centered' * confound_centered) * confound_centered';
                brain2 = brain - P_confound * brain;
                score2 = score - P_confound * score;
            end
            
            [~,brain_reduced] = pca(brain2, 'NumComponents', n_components);
            
            [brain_comp,score_comp,stat,~,~,stats] = canoncorr(brain_reduced,score2);
            p = stats.pChisq;

            if ~isempty(confounds)
                
                stat_fullres = stat;
                p_fullres = p;
                
                varargout{1} = stat_fullres;
                varargout{2} = p_fullres;

                switch test_type
                case 'multi_t2'
                    [~,brain_reduced_w_confound] = pca(brain, 'NumComponents', n_components); 
                    [brain_comp,score_comp,stat,~,~,stats] = canoncorr(brain_reduced_w_confound,score2);
                case 'multi_r'
                    [brain_comp,score_comp,stat,~,~,stats] = canoncorr(brain_reduced,score);
                end
            end

            if strcmp(test_type,'multi_t2')
                if ~isempty(confounds)
                    stat = r_to_test_stats(stat, n, 1, 2);
                    varargout{1} = r_to_test_stats(stat_fullres, n, 1, 2);
                else
                    stat = r_to_test_stats(stat, n, 0, 2);
                end
            end
    end
end