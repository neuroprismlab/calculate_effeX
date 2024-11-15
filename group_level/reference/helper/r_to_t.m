function t = r_to_t(r, n, n_confounds)
    % only makes sense for t from 2-sample t-test
    % n_confounds > 0 only for semipartial
    t = r .* ((n - 2 - n_confounds) ./ (1 - r .^2)) .^0.5;
end
