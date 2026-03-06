% MATLAB code to generate data for testing PLS class methods
%

clear
clc

data = xlsread("People.xlsx", 1, "B2:M33");
X = data;
X(:, 4) = [];
y = data(:, 4);

X = bsxfun(@minus, X, mean(X));
X = bsxfun(@rdivide, X, std(X));
y = bsxfun(@minus, y, mean(y));
y = bsxfun(@rdivide, y, std(y));

[XL, YL, XS, YS, beta, pctvar, mse, stats] = plsregress(X, y, 4);

csvwrite('pls-xloadings.csv', XL)
csvwrite('pls-yloadings.csv', YL)
csvwrite('pls-coeffs.csv', beta(2:end))
csvwrite('pls-weights.csv', stats.W)
csvwrite('pls-xscores.csv', XS)
csvwrite('pls-yscores.csv', YS)
csvwrite('pls-xres.csv', stats.Xresiduals)
csvwrite('pls-yres.csv', stats.Yresiduals)
csvwrite('pls-expvar.csv', pctvar)

%yscrs =
%
%   11.1420    0.0915    0.0668   -0.0015
%    5.6355   -0.1266   -0.0208    0.0050
%    5.6355   -0.0683    0.0099    0.0308
%   -5.3774   -0.5884   -0.1461   -0.0654
%   -2.6242   -0.0651    0.0344   -0.0154
%   -1.2476    0.0568    0.0901    0.0140
%    2.8823   -0.4560   -0.1198   -0.0433

%% compute VIP scores
% code is taken from here: https://se.mathworks.com/matlabcentral/answers/443239-how-to-calculate-the-variable-importance-in-projection-from-outputs-of-plsregress

W0 = bsxfun(@rdivide,stats.W,sqrt(sum(stats.W.^2,1)));
sumSq = sum(XS.^2,1).*sum(YL.^2,1);
vipScores = sqrt(size(XL,1) * sum(bsxfun(@times,sumSq,W0.^2),2) ./ sum(sumSq,2));
disp(vipScores)

csvwrite('pls-vipscores.csv', vipScores)
