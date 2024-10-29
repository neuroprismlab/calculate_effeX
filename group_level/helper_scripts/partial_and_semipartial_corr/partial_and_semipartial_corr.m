function [r_part, t_part, p_part, r_semipart, t_semipart, p_semipart] = partial_and_semipartial_corr(y,x,z,n)
    % y is outcome, x is predictor, and z is confound, n is sample size

    r_pairwise = corr_generalized(y, x, z);
    r_yx = r_pairwise{1,2};
    r_yz = r_pairwise{1,3};
    r_xz = r_pairwise{2,3};

    r_part = r_partial_from_pairwise(r_yx, r_yz, r_xz);
    r_semipart = r_semipartial_from_pairwise(r_yx, r_yz, r_xz);

    [t_part, p_part] = r_to_test_stats(r_part, n, 1, 2);
    [t_semipart, p_semipart] = r_to_test_stats(r_semipart, n, 1, 2);

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

function [r_semipart, t_semipart, p_semipart] = r_semipartial_from_pairwise(r_yx, r_yz, r_xz)
    % y is outcome, x is predictor, z is confound (note: denominator does not include outcome)
    % https://www.listendata.com/2017/03/partial-correlation.html
    r_semipart = (r_yx - r_yz .* r_xz) ./ (1 - r_xz .^2) .^0.5;
end

function [r_part, t_part, p_part] = r_partial_from_pairwise(r_yx, r_yz, r_xz)
    % y is outcome, x is predictor, z is confound (note: denominator does not include predictor x outcome)
    % https://www.listendata.com/2017/03/partial-correlation.html
    r_part = (r_yx - r_yz .* r_xz) ./ ((1 - r_xz .^2) * (1 - r_yz .^2)) .^0.5;
end

function [t,p] = r_to_test_stats(r, n, n_confounds, n_groups)
    % for test, following Kim 2015 (creator of ppcor) - validated in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9230004/ 
    % two-sided

    % TODO: test whether same as directly from r -> p
    t = r .* ((n - 2 - n_confounds) ./ (1 - r .^2)) .^0.5;
    p = 2 * (1 - tcdf(abs(t), (n - n_groups - n_confounds) )); % note: tcdf is vectorized

end


