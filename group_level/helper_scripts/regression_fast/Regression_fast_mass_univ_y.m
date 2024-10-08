function[T,P] = Regression_fast_mass_univ_y(x,y)

% INPUTS
% x:    nxp regressor matrix of n observations (rows) across
%       p regressors (columns)
% y:    Comun vector of responses 
% Pval: Optional, computes P-values of the estimates

% OUTPUT
% Val:  Column vector of the estimated coefficients (Pval == false/zero/missing), or
%       nx3 matrix of the estimated coefficients (first columsn) with p-values (second
%       column) and t-statistics (third column)

% USE
% Val = Regression_fast(x,y)

% Lightly edited by Stephanie Noble to return mass univariate estimates for multiple y simultaneously, to return t-stats, and small speedups (matrix left division, no logic, etc)

B   = (x'*x)\x'*y;

    y_fit = x*B;
    df  = -diff(size(x));

    s2 = sum((y-y_fit).^2)/df;
    %s   = (sum((y-y_fit).^2)/df)^0.5;
    se = sqrt(diag((x' * x) \ eye(size(x, 2))) * s2);
    T   = B./se;
   
    P = 2 * (1 - tcdf(abs(T), df));

end
