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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Setup

% script and data directories

%motion_dir = '/MOTIONPATH/'; % USER-DEFINED - REMOVED 051624
scripts_dir = '/work/neuroprism/effect_size/scripts/generate_studies/helper_scripts/'; 
%scripts_dir = '/home/USERNAME/scripts/effect_size/';
results_dir = '/work/neuroprism/effect_size/data/individual_studies/';
data_filename = '/work/neuroprism/data_shared/ukb/ukb_data_steph.mat'; % ADDED 051624

% motion paths - USER-DEFINED - REMOVED 051624

%motion_file = [motion_dir,'EXAMPLEPATH_allsubs_Movement_RelativeRMS_mean.txt'];
%motion_subids_file = [motion_dir,'subids']; % subject IDs for motion measurements

% script paths

%data_loader_script_path = [scripts_dir,'data_loader']; % USER-DEFINED [this is only a placeholder and does not exist] - REMOVED 051624
atlas_tools_script_path = [scripts_dir,'atlas_tools']; % for pooling (structure_data, summarize_data)
regression_fast_script_path = [scripts_dir,'regression_fast']; % for fast multiple regression

% params

pooling_params = 0;
motion_method_params = {'none', 'regression', 'threshold'};
low_motion_threshold = 0.1; % empirically, 5.7% of subjects are >0.1mm mFFD
n_network_groups = 10; % hard-coded for Shen atlas-based pooling

% set paths

%addpath(data_loader_script_path); % USER-DEFINED
addpath(genpath(atlas_tools_script_path));
addpath(regression_fast_script_path);


%% Load data and subject ID mapping

% USER-DEFINED - modify this section to load your data
% load the following:
%   - m (data matrices): should be n_sub x n_var 
%       - IMPORTANT: for summarization and pooling to work correctly, each row of m should be taken from the upper triangle of the original matrix consistent with Matlab indexing.
%         For example, this will be the case is m is created and saved in Matlab as follows:
%           triu_mask=logical(triu(ones(268),1));
%           m(this_subject,:) = mat(triu_mask);
%   - beh (behavior / phenotype to predict): should be n_sub x 1
%   - subids_data (subject IDs for the matrices / images): should be n_sub x 1


%[m, beh, subids_data] = data_loader;

this_beh='age';

data_matfile = load(data_filename);
m = data_matfile.fc;
beh = data_matfile.(this_beh); % age, fluid intelligence, gender
mean_motion = data_matfile.('framewise_displacement');

%return;

% output path - USER-DEFINED

results_file_prefix = [results_dir,'ukb_fc_t_rest_age__v2']; % see naming convention
%results_file_prefix = [results_dir,strrep(['ukb_fc_test_rest_age__v2'], ' ', '_')']; % see naming convention


%% Calculate effects

for do_pooling = pooling_params

    %% Do large-scale pooling if specified

    if do_pooling

        m2 = []; 
        triumask=logical(triu(ones(n_network_groups)));  
        for i = 1:size(m,1) % over subjects
            t = structure_data(m(i,:),'triangleside','upper');
            t2 = t; % HALLEE removed this line because we don't have 55 node UKB map - summarize_matrix_by_atlas(t,'suppressimg',1)'; % transpose because it does tril by default
            m2(i,:) = t2(triumask);
        end

        % update results file prefix
        results_file_prefix2 = [results_file_prefix,'__by_net'];

    else
        results_file_prefix2 = results_file_prefix;
        m2 = m;
    end


    for motion_method_it = 1:length(motion_method_params)

        %% Account/correct for motion as specified

        motion_method = motion_method_params{motion_method_it};

        beh2 = beh; % changes if applying motion threshold
        
        if ~strcmp(motion_method,'none')

%            % Align motion and data by subject - REMOVED 051624
%            
%            % align motion subIDs to data subIDs (removes subjects who do not have data or behavior)
%            subids_motion = load(motion_subids_file);
%            [subids,idx_matrices,idx_motion] = intersect(subids_data,subids_motion,'stable');
%            if ~isequal(subids_data,subids_motion(idx_motion)); error('Matrix and motion intersected subject IDs don''t match. Likely subjects with matrices are missing motion.'); end;
%            
%            % load motion and index based on above alignment
%            mean_motion = load(motion_file);
%            mean_motion = mean_motion(idx_motion); % get subs with matrices
%
%            std_mean_motion = std(mean_motion); % TODO: save this or full motion vector
%            mean_mean_motion = mean(mean_motion); % TODO: save this or full motion vector
 
            % threshold if specified
            if strcmp(motion_method,'threshold')
                low_motion_idx = find(mean_motion<low_motion_threshold); % TODO: consider saving
                m2 = m2(low_motion_idx,:);
                beh2 = beh2(low_motion_idx);
                std_mean_motion = std(mean_motion(low_motion_idx)); % TODO: save this or trimmed motion vector
                mean_mean_motion = mean(mean_motion(low_motion_idx)); % TODO: save this or trimmed motion vector
 
            end

            % update results file prefix
            results_file_prefix3=[results_file_prefix2,'__with_motion_',motion_method];
        else
            results_file_prefix3 = results_file_prefix2;
        end


        %% Run regression

        if strcmp(motion_method,'regression')
            % include motion as a confound
            [r,p,n,std_X,std_y] = save_univariate_regression_results(m2,beh2,results_file_prefix3, mean_motion);
        else
            [r,p,n,std_X,std_y] = save_univariate_regression_results(m2,beh2,results_file_prefix3);
        end
    
    end

end


%% Function definitions

% set up function for univariate y~x with optional confound as predictor

% NOTE: it may be tempting to use a stepwise procedure - estimate the residuals of m ~ motion, then use those residuals for subsequently estimating betas - but residualizing in this way biases effect size estimates towards zero (i.e., conservative), particularly with higher collinearity in predictors. Instead, including the confound in the actual model should be preferred as unbiased (though also note higher variance with collinearity) - https://besjournals.onlinelibrary.wiley.com/doi/10.1046/j.1365-2656.2002.00618.x
%   i.e., avoid: mdl = fitlm(X, confound); X = mdl.Residuals;

% TODO: either convert the above to point-biserial correlation coefficient or t-statistic or add a function to calculate these

function [r,p,n,std_X,std_y] = save_univariate_regression_results(X,y,results_file_prefix,X2)
    % X: n_sub x n_var, y: n_sub x 1, results_file_prefix: string, Optional X2: n_sub x n_var
    n = length(y);
    std_X = std(X);
    std_y = std(y);
    
    if nargin==4
        for i=1:size(X,2)
            mdl(:,:,i) = Regression_fast([ones(n,1),X(:,i),X2], y, 1); % note: fitlm is built-in for this but too slow % TODO: check p-value calculation - some set to 0,  maybe singular for edge-wise
        end
        b = squeeze(mdl(2,1,:))'; % unstandardized betas
        r = b.*std_X / std_y; % standardized betas - https://www3.nd.edu/~rwilliam/stats1/x92.pdf - TODO: this is not technically a "partial r" - decide what to do for subsequent R^2 or Cohen's d
        p = squeeze(mdl(2,2,:))';
    elseif nargin==3
        [r,p] = corr(X,y);
    else
        error('%d arguments provided but only 3 or 4 allowed.',nargin)
    end
    
    save([results_file_prefix,'.mat'],'r','p','std_X','std_y','n')
end

