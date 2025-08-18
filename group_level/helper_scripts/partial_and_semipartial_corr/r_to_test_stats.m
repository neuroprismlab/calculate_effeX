function [t,p] = r_to_test_stats(r, n, n_confounds, n_groups)
    % for test, following Kim 2015 (creator of ppcor) - validated in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9230004/ 
    % two-sided

    % TODO: test whether same as directly from r -> p
    t = r .* ((n - 2 - n_confounds) ./ (1 - r .^2)) .^0.5;
    p = 2 * (1 - tcdf(abs(t), (n - n_groups - n_confounds) )); % note: tcdf is vectorized

end


