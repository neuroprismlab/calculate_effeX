function [T,P,B] = Regression_faster_mass_univ_y(x,y,coeff_id)

% INPUTS
% x:    nxp regressor matrix of n observations (rows) across
%       p regressors (columns)
% y:    nxm matrix of responses, for running m unique mass univariate tests
% coeff_id: Optional, returns specific coefficient

% OUTPUT
% T: pxm estimated t-statistics (or <p x m if coeff_id specified)
% P: pxm estimated p-values (or <p x m if coeff_id specified)

% USE
% [T,P] = Regression_fast(x,y)

% Edited by Stephanie Noble (Aug 2024) to return mass univariate estimates
% for multiple y simultaneously and to only return select estimates

B   = (x'*x)\x'*y;
y_fit = x*B;
df  = -diff(size(x));

residuals = y - y_fit;
s2 = sum(residuals .^ 2) / df;
se = sqrt(diag((x' * x) \ eye(size(x, 2))) * s2);
T = B./se;

P = 2 * (1 - tcdf(abs(T), df));

if nargin == 3
    T = T(coeff_id,:);
    P = P(coeff_id,:);
end

end