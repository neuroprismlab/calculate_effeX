%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference for statistic and effect size calculation
%
% Calculates:
%   stats (r/t/r_semipart, p)
%   effect sizes + CIs (d, r^2/eta^2)
%
% For validation, see effect_size_conversion_ref_validation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[current_dir,~,~] = fileparts(mfilename('fullpath'));
scripts_dir = [current_dir,'/../helper_scripts/'];
scripts_dir2 = [current_dir,'/helper/'];
addpath(genpath(scripts_dir))
addpath(genpath(scripts_dir2))

clear all

% user-defined
% stat_list = {'t'};
% stat_list = {'t','t2','r'};
stat_list = {'multi_t'};
use_confound_opts = [0,1];
use_confound = 1;
validate_conversion = 1;
add_effect = 1;
add_effect_confound = 1;
use_diff_group_sizes = 0;

% setup
num_sdx_r_to_d = 2;
n = 200;
m = 10000; % m = 1;

z = randn([n,1]); % confound

if add_effect
    w_effect=1;
    intercept = 0.5;
else
    w_effect=0;
    intercept = 0;
end

if add_effect_confound
    w_confound=1;
else
    w_confound=0;
end


for stat_idx = 1:length(stat_list)

    stat = stat_list{stat_idx};

    % finish setting up data

    if strcmp(stat,'t2') || strcmp(stat,'multi_t2')
    if use_diff_group_sizes
        n_diff_btw_groups = 10;
    else
        n_diff_btw_groups = 0;
    end
    n1 = (n+n_diff_btw_groups) / 2;
    n2 = (n-n_diff_btw_groups) / 2;
    x = [zeros(n1,1); ones(n2,1)]; % score
    else
        x = randn([n,1]); % score
    end


    y = intercept + w_effect*x/10 + w_confound*z + randn([n,m]); % brain

    for use_confound = use_confound_opts

        switch stat
        
            case 't'
        
                % stat & test
                if use_confound

                    % 1. multiple regression logic: effect size represents add'l variance explained over null (H0: intercept=0)
                    [t,p,B] = Regression_faster_mass_univ_y([ones(n,1), z], y, 1);
                    n_coefficients = 2; % TODO: r combine script: will have to hard code (or pass this)

                    % TODO
                    % 2. zero-confound logic: effect size represents
                    % using a variation on partial corr logic (pre-regress out of x and y) - instead, regress out of y (y_resid ~ intercept) - effect size from t-stat implicitly measures the decrease in SS with/without intercept
                    y_residuals_with_intercept = y - [ones(n,1), z] * B + B(1);
                    [t_fullres,p_fullres] = Regression_faster_mass_univ_y(ones(n,1), y_residuals_with_intercept);

                    stat2 = t_fullres;
                    p2 = p_fullres;

                    % TODO: validate
                    d_fullres = t1_to_d(t_fullres,n);
                    % d = B(1) / std(y_residuals); % TODO: move to validation
                    omega_sq_fullres = t1_to_omega_sq(t_fullres,n,1); % TODO: check n_coefficients

                else
                    [t, p] = Regression_faster_mass_univ_y(ones(n,1), y);
                    n_coefficients = 1; % TODO: r combine script: same as above
                end
                
                stat = t;

                % effect sizes
                % if using confound, this is the expected effect size when that confound
                % (e.g., motion) is 0
                d = t1_to_d(t,n);
                omega_sq = t1_to_omega_sq(t,n,n_coefficients); % proportion of the total variance in y that is attributed to the fact that y_mean is not null (0)
        
                % CI - IN R
                % d_ci <- sapply(d, function(x) d_ci(x, n1 = n, alpha = 0.05))
                % omega_sq_ci <- ci.R2(R2 = omega_sq, N = n, conf.level = 0.95*n_var) % via MBESS
        
            case 't2'
        
                if use_confound
        
                    % stat & test

                    [r_part, t_part, p_part, r_semipart, t_semipart, p_semipart] = partial_and_semipartial_corr(y,x,z,n);

                    stat = t_semipart;
                    stat2 = t_part;
                    p = p_semipart;
                    p2 = p_part;

                    % [r_part2, p_part2] = partialcorr(x,y,z); % Note: calculating *both* semipartial and partial above is faster than calculating just partial using this built-in

                    % effect size
                    % let's convert to t to get more precise d for uneven groups
                    % t_semipart = r_to_t(r_semipart,n,1);
                    d = t2_to_d(r_semipart,n1,n2);
                    r_sq_semipart = r_semipart.^2;

                    % t_part = r_to_t(r_part,n,1);
                    d_part = t2_to_d(r_part,n1,n2);
                    r_sq_part = r_part.^2;
        
                    % CI - IN R
                    % TODO
                    % r_sq_full = r_sq_from_pairwise(r_yx, r_yz, r_xz);
                    % r_sq_semipart_CI = ci.spcor(0.05/n_var, r_semipart, r_sq_full, n); % https://rdrr.io/github/dgbonett/statpsych/man/ci.spcor.html
                    % d_semipart_CI <- r_to_d(sqrt(r_sq_semipart_CI));
                    % TODO: be good to check that you get same CI result from t_semipart or straight from d
        
                else
        
                    % stat & test
                    [t, p] = Regression_faster_mass_univ_y([ones(n,1), x], y, 2);
        
                    stat = t;

                    % effect size
                    d = t2_to_d(t, n1, n2); % for unequal group sizes
                    r = t_to_r(t, n, 0);
                    r_sq = r.^2;
        
                    % CI - IN R
                    % d_ci <- sapply(d, function(x) d_ci(x, n1 = n/2, n2 = n/2, alpha = 0.05/n_var))
                    % r_sq_ci <- ci.r_sq(r_sq = r^2, N = n, conf.level = 0.95*n_var)
        
                end
        
            case 'r'
        
                if use_confound
        
                    % stat & test

                    [r_part, t_part, p_part, r_semipart, t_semipart, p_semipart] = partial_and_semipartial_corr(y,x,z,n);
                    
                    stat = r_semipart;
                    stat2 = r_part;
                    p = p_semipart;
                    p2 = p_part;

                    % r_pairwise = corr_generalized(y, x, z);
                    % r_yx = r_pairwise{1,2};
                    % r_yz = r_pairwise{1,3};
                    % r_xz = r_pairwise{2,3};
                    % 
                    % [r_semipart, p_semipart] = r_semipartial_from_pairwise(r_yx, r_xz, r_yz, n); % now the *outcome* is x (score)
                    % [r_part, p_part] = r_partial_from_pairwise(r_yx, r_yz, r_xz, n);

                    % effect size
                    d_semipart = r_to_d(r_semipart, num_sdx_r_to_d);
                    r_sq_semipart = r_semipart .^2;
                    
                    d_part = r_to_d(r_part, num_sdx_r_to_d);
                    r_sq_part = r_part.^2;
        
                    % CI - IN R
                    % r_sq_full = r_sq_from_pairwise(r_yx, r_xz, r_yz); % now the *outcome* is x (score)
                    % r_sq_semipart_CI = ci.spcor(alpha = 0.05*n_var, r_semipart, r_sq_full, n); % https://rdrr.io/github/dgbonett/statpsych/man/ci.spcor.html
                    % d_semipart_CI = r_to_d(sqrt(r_sq_semipart_CI));
        
                else
        
                    % stat & test
                    [t, p] = Regression_faster_mass_univ_y([ones(n,1), x], y, 2); % using y as outcome, but result identical to mass univariate with x as outcome
                    r = t_to_r(t, n, 0);

                    stat = r;
        
                    % effect size
                    d = r_to_d(r, num_sdx_r_to_d);
                    r_sq = r.^2;
        
                    % CI
                    %d_ci <- sapply(r, function(x) d_ci__from_r(x, n = d_maps[[i]][[t]]$n[1], num_sdx_r_to_d = num_sdx_r_to_d, alpha = 0.05/n_var))
                    %library(MBESS)
                    %ci_r_sq <- ci.r_sq(r_sq = r^2, N = n, conf.level = 0.95*n_var)
                end
                
            case 'multi_t'
        
                    % 1. Dimensionality reduction - slow (~10 sec)
                    y_reduced = dimensionality_reduction(y, n, 1);
                    
                    % 2. Optional: regress z from y
                    if use_confound
                        % z = z(complete_cases,:); % filter for only complete cases
                        confound_centered = z - mean(z);
                        P_confound = confound_centered / (confound_centered' * confound_centered) * confound_centered'; % confound projection matrix
                        y2 = y_reduced - P_confound * y_reduced;

                        % TODO - follow t above
                        % using a variation on partial corr logic (pre-regress out of x and y): instead, regress out of y (y_resid ~ intercept) - effect size from t-stat implicitly measures the decrease in SS with/without intercept

                    else
                        y2 = y_reduced;
                    end
        
                    % 3. Hotelling t-test
                    % stats & effect size
                    [t_sq, p] = Hotelling_T2_Test(y2); % TODO: consider epsilon_sq - pros: should be more unbiased - cons: other measures are using omega_sq
                    t = sqrt(t_sq);

                    stat = t;
        
                    % more effect sizes
                    d = t1_to_d(t, n);
                    if use_confound; n_confounds=1; else n_confounds=0; end
                    r = t_to_r(t, n, n_confounds); % treating it like the null is a reference population
                    r_sq = r.^2;

                    % TODO: add hierarchical - IN PROGRESS
        
                    % PLSR with/without intercept
                    [~, ~, ~, ~, ~, ~, MSE_null] = plsregress(z, y_reduced, 1, 'Intercept', false);
                    [~, ~, ~, ~, ~, ~, MSE_alt] = plsregress(z, y_reduced, 1, 'Intercept', true);
                    
                    % t-test
                    MSE_diff = MSE_null(2,2) - MSE_alt(2,2); % rows are x,y; col are no components (intercept only), 1st component etc
                    n_coeff = 1;
                    SE_diff = sqrt((2 * MSE_alt(2,2)^2) / (n - n_coeff - 1));
                    t_stat = MSE_diff / SE_diff;

                    % eta_sq = 1 - MSE_alt / MSE_null;
                    
                    % first 2 terms are the loadings, second 2 are the
                    % scores, 5 is regression coefficients, 6 is percentage
                    % of variance explained, 7 is MSE, 8 is addl stats

                    % CI - IN R
                    % d_ci <- sapply(d, function(x) d_ci(x, n1 = n, alpha = 0.05))
                    % omega_sq_ci <- ci.R2(R2 = omega_sq, N = n, conf.level = 0.95*n_var) % TODO: validate
                    % TODO: VALIDATE this is appropriate for multivariate
        
        
            case {'multi_t2', 'multi_r'}
        
                    % Canonical Correlation (Multivariate): x is predictor, y is outcome (equivalent to the opposite for t-test analogue)
        
                    % 1. Dimensionality reduction - slow (~10 sec)
                    y_reduced = dimensionality_reduction(y, n, 0);
        
                    % 2. Optional: regress z from x and y
                    if use_confound
                        % z = z(complete_cases,:); % filter for only complete cases
                        confound_centered = z - mean(z);
                        P_confound = confound_centered / (confound_centered' * confound_centered) * confound_centered'; % confound projection matrix
                        y2 = y_reduced - P_confound * y_reduced;
                        x2 = x - P_confound * x;

                        % 3. Canonical Correlation (top component)
                        [y_comp_part, x_comp_part, r_part, ~, ~, stats_part] = canoncorr(y2,x2);
                        p_part = stats_part.pChisq;

                        [y_comp_semipart, x_comp_semipart, r_semipart, ~, ~, stats_semipart] = canoncorr(y,x2);
                        p_semipart = stats_semipart.pChisq;

                        p = p_semipart;
                        p2 = p_part;
                        if strcmp(stat,'multi_t2')
                            stat = r_to_t(r_semipart);
                            stat2 = r_to_t(r_part);
                        else
                            stat = r_semipart;
                            stat2 = r_part;
                        end


                    else
                        % 3. Canonical Correlation (top component)
                        [y_comp, x_comp, r, ~, ~, stats] = canoncorr(y_reduced,x);
                        p = stats.pChisq;

                        if strcmp(stat,'multi_t2')
                            stat = r_to_t(r);
                        else
                            stat = r;
                        end
                    end

        
                    % effect size
                    % TODO: currently using conversion for univariate - would it be preferable to use partial_eta_sq? or redundancy index?
                    % d_part = r_to_d(r_part, num_sdx_r_to_d);
                    % r_sq_part = r_part.^2;
                    %
                    % d_semipart = r_to_d(r_semipart, num_sdx_r_to_d);
                    % r_sq_semipart = r_semipart.^2;
        
                    % CI IN R
                    %d_ci <- sapply(r, function(x) d_ci__from_r(x, n = d_maps[[i]][[t]]$n[1], num_sdx_r_to_d = num_sdx_r_to_d, alpha = 0.05/n_var))
                    %library(MBESS)
                    %ci_r_sq <- ci.r_sq(r_sq = r^2, N = n, conf.level = 0.95*n_var)
                
        end
        
        if validate_conversion
            effect_size_ref_validation
        end
        
        

    end

    % clearvars n m num_sdx_r_to_d x x y z mdl  use_confound

end

%% Conversion functions
% here, the test statistic is always a vector

% Test stat

function r = t_to_r(t, n, n_confounds)
    % only makes sense for t from 2-sample t-test
    % n_confounds > 0 only for semipartial % TODO: and partial - correct other references
    r = t ./ (t .^2 + (n - 2 - n_confounds)) .^0.5;
end

function t = r_to_t(r, n, n_confounds)
    % only makes sense for t from 2-sample t-test
    % n_confounds > 0 only for semipartial
    t = r .* ((n - 2 - n_confounds) ./ (1 - r .^2)) .^0.5;
end

function p = t_to_p(t, n, n_confounds, n_groups)
    % n_confounds > 0 only for semipartial
    % two-sided
    p = 2 * (1 - tcdf(abs(t), (n - n_groups - n_confounds) )); % note: tcdf is vectorized
end

function r = corr_generalized(varargin)
    % generalized pairwise correlation for an arbitrary number of input variables
    n = nargin;
    r = cell(n, n);
    for i = 1:n
        for j = i+1:n
            r{i, j} = corr(varargin{i}, varargin{j});
        end
    end
end

function [r_xy, r_xz, r_yz] = mass_univariate_correlation(x, y, z)
    % pairwise correlation between y, x, and z
    r_xy = corr(x,y);
    r_xz = corr(x,z);
    r_yz = corr(y,z);
end
function [r_yx, r_yz, r_xz] = mass_univariate_correlation__from_regfast(y, x, z)
    % pairwise correlation between y, x, and z, where y contains multiple columns (note:
    % avoids pairwise corr btw columns of y)

    n = length(z);
    [t_yx, ~] = Regression_faster_mass_univ_y([ones(n,1), x], y, 2);
    [t_yz, ~] = Regression_faster_mass_univ_y([ones(n,1), z], y, 2);
    [t_xz, ~] = Regression_faster_mass_univ_y([ones(n,1), x], z, 2);
    
    df = n - 2; % df = n-2
    r_yx = t_yx ./ (t_yx .^2 + df) .^0.5;
    r_yz = t_yz ./ (t_yz .^2 + df) .^0.5;
    r_xz = t_xz ./ (t_xz .^2 + df) .^0.5;
end

function [r_semipart, p_semipart] = r_semipartial_from_pairwise(r_yx, r_yz, r_xz, n)
    % y is outcome, x is predictor, z is confound (note: denominator does not include outcome)
    % https://www.listendata.com/2017/03/partial-correlation.html
    
    r_semipart = (r_yx - r_yz .* r_xz) ./ (1 - r_xz .^2) .^0.5;

    % for test, following Kim 2015 (creator of ppcor) - validated in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9230004/ 
    t_semipart = r_to_t(r_semipart, n, 1);
    p_semipart = t_to_p(t_semipart, n, 1, 2);
end

function [r_part, p_part] = r_partial_from_pairwise(r_yx, r_yz, r_xz, n)
    % y is outcome, x is predictor, z is confound (note: denominator does not include outcome)
    % https://www.listendata.com/2017/03/partial-correlation.html
    
    r_part = (r_yx - r_yz .* r_xz) ./ ((1 - r_xz .^2)*(1 - r_yz .^2)) .^0.5;

    % for test, following Kim 2015 (creator of ppcor) - validated in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9230004/ 
    t_part = r_to_t(r_part, n, 1);
    p_part = t_to_p(t_part, n, 1, 2);
end

function r_sq = r_sq_from_pairwise(r_yx, r_yz, r_xz)
    % y is outcome, x is predictor, z is confound (note: denominator does not include outcome)
     % https://stats.stackexchange.com/questions/32294/regression-r_sq-and-correlations
    r_sq = (r_yx.^2 + r_yz.^2 - 2 * r_yx .* r_yz .* r_xz) ./ (1 - r_xz.^2);
end

% Effect size

function d = r_to_d(r, num_sdx)
    d = num_sdx * r ./ (1 - r.^2) .^ (1/2);
end

function d = t1_to_d(t, n)
    % for 1-sample t-test
    d = t ./ sqrt(n);
end

function d = t2_to_d(t, n1, n2)
     % for 2-sample t-test
     % allows unequal group sizes
    d = t .* sqrt(1 / n1 + 1 / n2);
end

function omega_sq = t1_to_omega_sq(t,n,n_coefficients)
    % for 1-sample t-test
    % omega_sq is less biased estimator of eta-squared
    omega_sq = (t.^2 - 1) ./ (t.^2 + (n-n_coefficients));
    % eta_sq = t.^2 ./ (t.^2 + (n-1));
end

% Multivariate

function y_reduced = dimensionality_reduction(y, n, center)
    % aiming for 50 samples/feature for stable results a la Helmer et al.

    n_components = floor(n/50); 
    n_vars = size(y,2);
    if n_components > n_vars
        n_components = n_vars;
    end

    if center
        % make sure not to center - we're measuring dist from 0!
        [~,y_reduced] = pca(y, 'NumComponents', n_components, 'Centered', 'off');
    else
        [~,y_reduced] = pca(y, 'NumComponents', n_components);
    end

end

