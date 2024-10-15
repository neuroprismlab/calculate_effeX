%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Validate statistic and effect size calculation as defined by reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allowable error rates - user-defined
allowable_error_t = 0.002;
allowable_error_r = 0.002;
allowable_error_r_sq = 0.004;
allowable_error_p = 0.002;
allowable_error_d = 0.015;
allowable_error_omega_sq = 0.015;

if use_confound
    confound_str = 'with confound';
else
    confound_str = 'without confound';
end
fprintf('** Test %s (%s) **\n',stat,confound_str)

switch stat

    case 't'

        if use_confound

            % stat & test from older script
            % confirms results are consistent with previous version
            [~, p_ref, ~, stats] = ttest(y);
            t_ref = stats.tstat;
            d_ref = t_ref * sqrt(n);

            % manual test 1: difference in SS between models
            % following logic in https://www.graphpad.com/support/faqid/918/
            SS_alt = sum((y-mean(y)).^2);
            SS_null = sum(y.^2);
            omega_sq_no_deconfounding = 1 - SS_alt/SS_null; % this shouldn't be different
            
            x = z;
            B = (x' * x) \ (x' * y);
            y_fit = x * B;
            df = size(x, 1) - size(x, 2);
            residuals = y - y_fit;
            SS_null = sum(residuals .^ 2) / df;

            x = [ones(size(x,1),1), z];
            B = (x' * x) \ (x' * y);
            y_fit = x * B;
            df = size(x, 1) - size(x, 2);
            residuals = y - y_fit;
            SS_alt = sum(residuals .^ 2) / df;

            omega_sq_ref = 1 - SS_alt/SS_null;




            % manual test 2: regressing out confound

            % omega_sq from d
            omega_sq_ref2 = d.^2 / (d.^2 + n);
            % eta_sq_ref2 = (d.^2 - 1) / (d.^2 + (n / (n-1)));

            % compare
            compare_stats(t,t_ref,allowable_error_t)
            compare_stats(p,p_ref,allowable_error_p)
            compare_stats(omega_sq, omega_sq_ref, allowable_error_omega_sq) % MAIN TEST: from raw data
            compare_stats(omega_sq, omega_sq_ref2, allowable_error_omega_sq)  % NO MATCH: not sure why, but

        else

            % stat & test from built in
            [~, p_ref, ~, stats] = ttest(y);
            t_ref = stats.tstat;
            d_ref = t_ref * sqrt(n);

            % d from raw
            d_ref = mean(y) / std(y); % b = mean(y)

            % equivalent definition as R^2
            % https://www.graphpad.com/support/faqid/918/
            SS_alt = sum((y-mean(y)).^2);
            SS_null = sum(y.^2);
            omega_sq_ref = 1 - SS_alt/SS_null;

            % omega_sq from d
            omega_sq_ref2 = d.^2 / (d.^2 + n);
            % eta_sq_ref2 = (d.^2 - 1) / (d.^2 + (n / (n-1)));

            % compare
            compare_stats(t,t_ref,allowable_error_t)
            compare_stats(p,p_ref,allowable_error_p)
            compare_stats(omega_sq, omega_sq_ref, allowable_error_omega_sq) % MAIN TEST
            compare_stats(omega_sq, omega_sq_ref2, allowable_error_omega_sq) % NO MATCH - allowed to be somewhat different
            compare_stats(d,d_ref,allowable_error_d) % NO MATCH - this might be allowable to be different

        end


    case 't2'

        if use_confound

            % reduce dimensions for speedup
            m = 5;
            y = y(:,1:m);

            % r_semipart from dedicated function
             % only runs for 1-D vectors
            [r_semipart_ref,~] = msemipartialcorr(y(:,1),x,z);

            % r_semipart from corr - VALIDATION: confirmed it matches the above but unnecessarily computes corr between all columns of y, slowing things down by more than 100x for m=10,000
            [r,p_ref] = corr([y(:,1:m),x,z]);
            r_yx = r(1:m, m+1);
            r_yz = r(1:m, m+2);
            r_xz = r(m+1, m+2);
            r_semipart_ref2 = (r_yx - r_yz .* r_xz) ./ (1 - r_xz.^2) .^0.5;
            %r_semipart_ref = (r(1,2) - r(1,3) * r(2,3)) / sqrt(1 - r(2,3)^2); % where 2 is the predictor of interest and 3 is the confound

            % compare
            compare_stats(r_semipart, r_semipart_ref, allowable_error_r) % MAIN
            compare_stats(r_semipart, r_semipart_ref2, allowable_error_r) % MAIN

        else

            % d from t-test
            % only does first test, though (no mass univariate)
            [~, ~, ~, stats] = ttest2(y(x==1), y(x==0));
            t_ref = stats.tstat;
            d_ref = t_ref * sqrt(1 / n1 + 1 / n2);

            % r, p, and then t from built-in corr
            [r_ref2, p_ref2] = corr(x,y);
            t_ref2 = r_ref2 .* sqrt((n - 2) ./ (1 - r_ref2.^2));

            % d from r AND assuming equal group sizes:
            d_ref2 = num_sdx_r_to_d .* r ./ ((1 - r.^2) .^ (1/2));



            % compare
            compare_stats(t, t_ref, allowable_error_t) % MAIN
            compare_stats(t, t_ref2, allowable_error_t)
            compare_stats(d, d_ref, allowable_error_d)
            compare_stats(d, d_ref2, allowable_error_d) % for evaluating how big a difference using equal group sizes makes
            compare_stats(r, r_ref2, allowable_error_r) % MAIN
            compare_stats(p, p_ref2, allowable_error_p) % MAIN

        end

    case 'r'

        if use_confound

            % reduce dimensions for speedup
            m = 5;
            y = y(:,1:m);

            % r_semipart from dedicated function - VALIDATION: confirmed it matches the above (but only runs for univariate y):
            [r_semipart_ref,~] = msemipartialcorr(x,y(:,1),z);

            % r_semipart from corr - VALIDATION: confirmed it matches the above but unnecessarily computes corr between all columns of y, slowing things down by more than 100x for m=10,000
            [r,~] = corr([y(:,1:m),x,z]);
            r_yx2 = r(1:m, m+1);
            r_yz2 = r(1:m, m+2);
            r_xz2 = r(m+1, m+2);
            r_semipart_ref2 = (r_yx2 - r_yz2 * r_xz2) ./ (1 - r_yz2 .^2) .^0.5;

            % r_sq_semipart from model - VALIDATION: confirmed it matches the above:
            model1 = fitlm([x,z], y(:,1));
            r_sq_y_xz = model1.Rsquared.Ordinary;
            model2 = fitlm(z, y(:,1));
            r_sq_y_z = model2.Rsquared.Ordinary;
            r_sq_semipart_ref2 = r_sq_y_xz - r_sq_y_z;

            % compare
            compare_stats(r_semipart, r_semipart_ref, allowable_error_r) % MAIN
            compare_stats(r_semipart, r_semipart_ref2, allowable_error_r) % MAIN
            compare_stats(r_sq_semipart, r_sq_semipart_ref2, allowable_error_r_sq) % MAIN - confirm that it matches the expected definition from r^2s

        else

            % reduce dimensions for speedup
            m = 5;
            y = y(:,1:m);

            % stat & test
            [r_ref,p_ref] = corr(x,y);

            % compare
            compare_stats(r, r_ref, allowable_error_r) % MAIN
            compare_stats(p, p_ref, allowable_error_p) % MAIN

        end

    case 'multi_t'

        % d=t --> same result as from a direct estimate of d (sample):
        d_ref = sqrt(mahal(zeros(1,size(y2,2)),y2));
        % and same p as from manova1:
        [~, p_ref, ~] = manova1(y2, ones(size(y2, 1), 1));

        % TODO: should d alter df with 

        % compare
        compare_stats(d, d_ref, allowable_error_d) % MAIN
        compare_stats(p, p_ref, allowable_error_p) % MAIN

    case {'multi_t2', 'multi_r'}
        fprintf('No tests available for multi_t2 or multi_r \n')

         % TODO: compare p-value with look up from Winkler et al table
        
end

% TODO: check multi

function compare_stats(stat, stat_ref, allowable_error)

    if length(stat_ref) < length(stat)
        stat = stat(1:length(stat_ref));
    end

    fprintf('Checking stat %s ... ', inputname(2))
    mean_error = mean(abs(stat(:)-stat_ref(:)));
    if mean_error > allowable_error
        fprintf('** Error = %0.2f ** \n', mean_error)
    else
    fprintf('Matches \n', mean_error)
        
    end
end