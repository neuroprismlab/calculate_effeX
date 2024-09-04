

% Example and performance comparison

y = randn(100,1);
x = randn(100,4);

Val = Regression_fast(x,y);
Val = Regression_fast(x,y,false);


Val = Regression_fast(x,y,1);
Val = Regression_fast(x,y,true);

%% Comparison against fitlm, which however provides much more information
% about the regression, While e.g. only coefficients estimates are of interest.

N = 100;
y = randn(100,1);
x = randn(100,4);
xy = [y,x];

tic
for i = 1:N
    Regression_fast(x,y,1);
end
toc

tic
for i = 1:N
    fitlm(x,y);
end
toc

%% Comparsion against polyfit for the linear case

tic
for i = 1:N
    Regression_fast(x,y(:,1),0);
end
toc

tic
for i = 1:N
    polyfit(y,x(:,1),1);
end
toc

%% Comparison against regress, which however does not provide information on SE, and
% P-Values for the coefficients

tic
for i = 1:N
    Regression_fast(x,y,0);
end
toc

tic
for i = 1:N
    regress(y,x);
end
toc

%%


