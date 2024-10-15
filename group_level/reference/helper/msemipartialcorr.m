function [sr1,sr2] = msemipartialcorr(y,x1,x2)
% function [sr1,sr2] = msemipartialcorr(y,x1,x2)
% 
% computes the semipartial correlation between y and x1 (sr1) when x2 has been partialed out from x1
% in terms of a Venn diagram, if a is the area of exclusive overlap of x1 and y, b is the area
% of exclusive overlap of x2 and y, and c is the area where x1, x2, and y all overlap, the semi-
% partial correlation computed by this function corresponds to a
% 
% put differently, the effects of x2 have been removed from x1 but not from y (hence, "semi-" partial); 
% in this system, "removing the effect" is equivalent to subtracting from x1 the x1 values estimated from x2
% 
% the identical, but inverted argument, holds for sr2
% 
% validated with SPSS V20
% 
% May 2014, Maik C. Stuettgen
%% input check
if any(size(y)~=size(x1)) || any(size(y)~=size(x2))
  error('input vectors don''t match')
end
if size(y,1)~=length(y),y=y';end
if size(x1,1)~=length(x1),x1=x1';end
if size(x2,1)~=length(x2),x2=x2';end
%% computations
ry_x1  = corr(x1,y);
ry_x2  = corr(x2,y);
rx1_x2 = corr(x1,x2);
sr1 = (ry_x1-ry_x2*rx1_x2) / sqrt(1-rx1_x2^2);
sr2 = (ry_x2-ry_x1*rx1_x2) / sqrt(1-rx1_x2^2);