%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validate statistic and effect size calculation as defined by reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allowable error rates - user-defined
allowable_error_t = 0.005;
allowable_error_r = 0.003;
allowable_error_r_sq = 0.001;
allowable_error_p = 0.001;
allowable_error_d = 0.001;
allowable_error_omega_sq = 0.004;

if use_confound
    confound_str = 'with confound';
else
    confound_str = 'without confound';
end

fprintf(['\n----------------------------------------\n' ...
    'Checking test %s (%s)\n' ...
    '----------------------------------------\n'],stat,confound_str)

switch stat

    case 't'

        if use_confound

            % VALIDATION: confirm omega_sq consistent with our logic
            % Note: logic = reflects difference in SS vs. null models, following logic in https://www.graphpad.com/support/faqid/918/
            
            X = z;
            B = (X' * X) \ (X' * y);
            y_fit = X * B;
            df = size(X, 1) - size(X, 2);
            residuals = y - y_fit;
            SS_null = sum(residuals .^ 2) / df;

            X = [ones(size(x,1),1), z];
            B = (X' * X) \ (X' * y);
            y_fit = X * B;
            df = size(X, 1) - size(X, 2);
            residuals = y - y_fit;
            SS_alt = sum(residuals .^ 2) / df;

            omega_sq_ref = 1 - SS_alt ./ SS_null; % TODO: really this is eta_sq
            % omega_sq = (eta_sq - (df_effect / df_total) * (1 - eta_sq)) / (1 + (df_effect / df_total));

            compare_stats(omega_sq, omega_sq_ref, allowable_error_omega_sq, 1) % critical
            

            % motivating example: https://www.graphpad.com/support/faqid/918/
            % y = [4.5,5.6;3.7,6.4;5.3,6.4;5.4,6.0;3.9,5.7];
            % y_diff = y(:,2) - y(:,1);
            % SS_null = sum(y_diff.^2); % y ~ 0
            % SS_alt = sum((y_diff-mean(y_diff)).^2); % y ~ 1
            % eta_sq = 1-SS_alt/SS_null;
            
            
            % TODO
            % SIMILARITY: using a variation on semipartial corr logic (pre-regress out of x, equivalent to hierarchical MR) - instead, estimate effect size from decrease in SS with/without intercept
            % SIMILARITY: using a variation on partial corr logic (pre-regress out of x and y): instead, regress out of y (y_resid ~ intercept) - effect size from t-stat implicitly measures the decrease in SS with/without intercept

            % SIMILARITY: How different from pre-regressing out confound (leaving intercept)
            % TODO: different strategies so expect these to differ, but want to understand what exactly accounts for that difference

            y_deconfounded = residuals + repmat(B(1,:),n,1);
            [~, p_ref, ~, stats] = ttest(y_deconfounded);
            t_ref = stats.tstat;
            d_ref = t_ref ./ sqrt(n);

            compare_stats(t,t_ref,allowable_error_t)
            compare_stats(p,p_ref,allowable_error_p)
            compare_stats(d,d_ref,allowable_error_d)


            % VALIDATION: Compare with omega_sq <- d

            omega_sq_ref2 = (d.^2 * n - 1) ./ (d.^2 * n + n-1);
            compare_stats(omega_sq, omega_sq_ref2, allowable_error_omega_sq, 1)


        else


            % VALIDATION: compare with standard function

            [~, p_ref, ~, stats] = ttest(y);
            t_ref = stats.tstat;
            d_ref = t_ref ./ sqrt(n);

            compare_stats(t,t_ref,allowable_error_t, 1)
            compare_stats(p,p_ref,allowable_error_p, 1)
            compare_stats(d,d_ref,allowable_error_d, 1)


            % VALIDATION: compare with base logic for d

            d_ref2 = mean(y) ./ std(y); % b = mean(y)
            compare_stats(d,d_ref2,allowable_error_d, 1) % NO MATCH - this might be allowable to be different


            % VALIDATION: compare with base logic for R^2
            % https://www.graphpad.com/support/faqid/918/
            
            SS_alt = sum((y-mean(y)).^2);
            SS_null = sum(y.^2);
            eta_sq_ref = 1 - SS_alt./SS_null;

            %% TODO: convert to omega_sq
            % compare_stats(omega_sq, omega_sq_ref, allowable_error_omega_sq, 1)


            % VALIDATION: omega_sq from t
            % logic: diff between t->omega_sq and t->eta_sq is -1 in numerator; adjust the same way for d->omega_sq
            % eta_sq_from_t = t.^2 ./ (t.^2 + n-1);
            % eta_sq_from_d = (d.^2 * n) ./ (d.^2 * n + n-1);
            % https://haiyangjin.github.io/2020/05/eta2d/

            omega_sq_ref2 = (d.^2 * n - 1) ./ (d.^2 * n + n-1);
            compare_stats(omega_sq, omega_sq_ref2, allowable_error_omega_sq, 1)

        end


    case 't2'

        if use_confound

            % reduce dimensions for speedup
            m = 5;
            y = y(:,1:m);


            % VALIDATION: compare with standard function
            % note: only runs for 1-D vectors
            % note: this answers the question: how much more does group contribute that motion contributes?

            [r_semipart_ref,~] = msemipartialcorr(y(:,1),x,z);

            compare_stats(r_semipart, r_semipart_ref, allowable_error_r, 1)


            % VALIDATION: compare base logic - starting with corr
            % note: unnecessarily computes corr between all columns of y, slowing things down by more than 100x for m=10,000
            
            [r,p_ref] = corr([y,x,z]);
            r_yx = r(1:m, m+1);
            r_yz = r(1:m, m+2);
            r_xz = r(m+1, m+2);
            r_semipart_ref2 = (r_yx - r_yz .* r_xz) ./ (1 - r_xz.^2) .^0.5;
            %r_semipart_ref = (r(1,2) - r(1,3) * r(2,3)) / sqrt(1 - r(2,3)^2); % where 2 is the predictor of interest and 3 is the confound

            compare_stats(r_semipart, r_semipart_ref2, allowable_error_r, 1)


            % VALIDATION: compare with logic for partial corr

            X = [ones(size(x,1),1), z];
            B = (X' * X) \ (X' * y);
            y_fit = X * B;
            y_residuals = y - y_fit;
            B = (X' * X) \ (X' * x);
            x_fit = X * B;
            x_residuals = x - x_fit;
            mdl = fitlm(x_residuals,y_residuals(:,1));
            t_part_ref = mdl.Coefficients.tStat(2);
            compare_stats(t_part, t_part_ref, allowable_error_t, 1)

            % playing with comparing regressing z out of only y
            % [~, ~, ~, stats] = ttest2(y_residuals(x==1), y_residuals(x==0));
            % t_ref = stats.tstat;
            % d_ref = t_ref * sqrt(1 / n1 + 1 / n2);
            % compare_stats(d, d_ref, allowable_error_d)


            % SIMILARITY: compare the above with logic for multiple regression
            % note: this is for x -> y and z -> y
            % note: this also answers the question: what would the effect be in the absence of any confounds? But this further assumes no relationship between x, z

            [~,~,pval,~,stats] = stepwisefit([x,z],y(:,1),'Display','off');
            t_part_ref2 = stats.TSTAT(1);
            % mdl = fitlm([x,z],y(:,1));
            % t_ref2 = mdl.Coefficients.tStat(2);

            compare_stats(t_part_ref(1), t_part_ref2, allowable_error_t)


        else


            % VALIDATION: compare with standard function
            % note: only does first test (no mass univariate)

            [~, ~, ~, stats] = ttest2(y(x==1), y(x==0));
            t_ref = stats.tstat;
            d_ref = t_ref * sqrt(1 / n1 + 1 / n2);

            compare_stats(t, t_ref, allowable_error_t, 1)
            compare_stats(d, d_ref, allowable_error_d, 1)


            %  VALIDATION: compare with standard function + conversion logic

            [r_ref2, p_ref2] = corr(x,y);
            t_ref2 = r_ref2 .* sqrt((n - 2) ./ (1 - r_ref2.^2));
            compare_stats(r, r_ref2, allowable_error_r, 1)
            compare_stats(p, p_ref2, allowable_error_p, 1)
            compare_stats(t, t_ref2, allowable_error_t, 1)


            % SIMILARITY: compare with d <- r conversion assuming equal group sizes -  see how big difference equal group sizes makes
            
            d_ref2 = num_sdx_r_to_d .* r ./ ((1 - r.^2) .^ (1/2));
            compare_stats(d, d_ref2, allowable_error_d)

        end

    case 'r'

        if use_confound

            % reduce dimensions for speedup
            m = 5;


            % VALIDATION: compare with standard function (from exchange)
            % note: only runs for univariate y

            [r_semipart_ref,~] = msemipartialcorr(x,y(:,1),z);
            compare_stats(r_semipart, r_semipart_ref, allowable_error_r, 1)


            % VALIDATION: compare base logic - starting with corr
            % note: unnecessarily computes corr between all columns of y, slowing things down by more than 100x for m=10,000

            [r,~] = corr([y(:,1:m),x,z]);
            r_yx2 = r(1:m, m+1);
            r_yz2 = r(1:m, m+2);
            r_xz2 = r(m+1, m+2);
            r_semipart_ref2 = (r_yx2 - r_yz2 * r_xz2) ./ (1 - r_yz2 .^2) .^0.5;

            compare_stats(r_semipart, r_semipart_ref2, allowable_error_r, 1)


            % VALIDATION: compare base logic - starting with R^2's

            model1 = fitlm([y(:,1),z], x); % recall again that x (score) is outcome here
            r_sq_x_yz = model1.Rsquared.Ordinary;
            model2 = fitlm(z, x);
            r_sq_x_z = model2.Rsquared.Ordinary;
            r_sq_semipart_ref2 = r_sq_x_yz - r_sq_x_z;

            compare_stats(r_sq_semipart(1), r_sq_semipart_ref2, allowable_error_r_sq) % MAIN - confirm that it matches the expected definition from r^2s
            

            % SIMILARITY: compare with base logic for removing direct effect of confound on predictor (brain)
            % note: expecting x -> y and z -> x -> y OR z -> x, so just remove effect of z -> x (here, x = brain)
            % note: not expecting z -> y, though y may -> z - no theoretical reason to expect z to affect y in a non-trivial way (e.g., motion to influence IQ)
            % note: this answers the question: what would the effect be in the absence of any confounds?

            X = [ones(size(x,1),1), z];
            B = (X' * X) \ (X' * y(:,1:m));
            y_fit = X * B;
            y_residuals = y(:,1:m) - y_fit;
            [r_ref,p_ref] = corr(x,y_residuals);
            compare_stats(r_semipart(1:m), r_ref, allowable_error_r)

        else


            % VALIDATION: compare with standard function
            
            % reduce dimensions for speedup
            m = 5;
            y = y(:,1:m);

            [r_ref,p_ref] = corr(x,y);
            compare_stats(r(1:m), r_ref, allowable_error_r, 1)
            compare_stats(p, p_ref, allowable_error_p, 1)

        end

    case 'multi_t'

        % VALIDATION: compare with standard functions
        % TODO: d=t --> same result as from a direct estimate of d (sample) ?
        % TODO: should d change with df ?

        % note: y2 is already confound-regressed
        d_ref = sqrt(mahal(zeros(1,size(y2,2)),y2)); % mahal -> d
        [~, p_ref, ~] = manova1(y2, ones(size(y2, 1), 1)); % manova1 -> p

        compare_stats(d, d_ref, allowable_error_d, 1)
        compare_stats(p, p_ref, allowable_error_p, 1)



        % TODO: FOR COMPARISON: center then project orthogonal then re-add
        % intercept
        
        % Step 1: Center the brain data
        brain_mean = mean(brain);
        brain_centered = brain - brain_mean;
        
        % Step 2: Center the confounds
        confound_centered = confounds - mean(confounds);
        
        % Step 3: Compute the projection matrix for the centered confounds
        P_confound = confound_centered / (confound_centered' * confound_centered) * confound_centered';
        
        % Step 4: Project the centered brain data to be orthogonal to the centered confounds
        brain_centered_projected = brain_centered - P_confound * brain_centered;
        
        % Step 5: Re-add the mean of the original brain data
        brain2 = brain_centered_projected + brain_mean;

    case {'multi_t2', 'multi_r'}
        fprintf('No tests available for multi_t2 or multi_r \n')

         % TODO: compare p-value with look up from Winkler et al table
        
end

% TODO: check multi

function compare_stats(stat, stat_ref, allowable_error, varargin)

    if nargin > 3
        critical = varargin{1};
    else
        critical = 0;
    end

    if ~critical
        allowable_error = allowable_error * 2;
    end

    if length(stat_ref) < length(stat)
        stat = squeeze(stat(1:length(stat_ref)));
        stat_ref = squeeze(stat_ref);
    end

    
    mean_error = mean(abs(stat(:)-stat_ref(:)));

    fprintf('%s', inputname(2))

    if mean_error > allowable_error

        if critical
            fprintf('\t<strong>**** CRITICAL ERROR **** </strong> ')
        else 
            fprintf('\t** Minor error ** ')
        end

        stack = dbstack;
        if length(stack) > 1
            callerInfo = stack(2);
            callerLineNumber = callerInfo.line;
        else
            callerLineNumber = [];
        end
        fprintf('error = %0.4f (line no = %d) \n', mean_error, callerLineNumber)
        
    else
        fprintf('\n')
    end

end