function[Val] = Regression_fast_mass_univ_y(x,y,Pval)

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

% Lightly edited by Stephanie Noble (Aug 2024) to return mass univariate estimates for multiple y simultaneously, and also to return t-stats

if nargin == 2
    Pval = false;
end

B   = (x'*x)\x'*y;
Val = B;

if Pval    
    y_fit = x*B;
    df  = -diff(size(x));

    s2 = sum((y-y_fit).^2)/df;
    %s   = (sum((y-y_fit).^2)/df)^0.5;
    se  = (diag(inv(x'*x))*s2).^0.5;
    T   = B./se;
    
    P   = (T>=0).*(1 - tcdf(T,df))*2 + (T<0).*(tcdf(T,df))*2;
    
    Val = cat(3,B,P,T);
end

end
