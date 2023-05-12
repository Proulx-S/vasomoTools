function funTs = dtrnd4psd(funTs)
%Detrend time series for multitaper estimation of psd
%   Simple detrending using order-2 polynomials, the highest order that
%   can't fitting a sinwave.

%% Define polynomial regressors
tr = funTs.tr/1000;
t = 0:tr:(funTs.nframes-1)*tr;
X = [];
X(:,end+1) = ones(size(t)); X(:,end) = X(:,end)./norm(X(:,end));
X(:,end+1) = t; X(:,end) = X(:,end) - mean(X(:,end)); X(:,end) = X(:,end)./norm(X(:,end));
X(:,end+1) = t.^2; X(:,end) = X(:,end) - mean(X(:,end)); X(:,end) = X(:,end)./norm(X(:,end));

%% Fit
beta = funTs.vec'/X';

%% Remove from data
funTs.vec = funTs.vec - X*beta';

