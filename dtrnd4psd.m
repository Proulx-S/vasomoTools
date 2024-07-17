function [funTs,funTsMean] = dtrnd4psd(funTs,t,prcBOLDflag)
%Detrend time series for multitaper estimation of psd
%   Simple detrending using order-2 polynomials, the highest order that
%   can't fitting a sinwave.
if ~exist('prcBOLDflag','var') || isempty(prcBOLDflag)
    prcBOLDflag = 0;
end
if isstruct(funTs)
    tr = funTs.tr/1000;
    if isfield(funTs,'t')
        t = funTs.t;
    else
        t = 0:tr:(funTs.nframes-1)*tr;
    end
    if ~isfield(funTs,'vec') || isempty(funTs.vec)
        funTs = vol2vec(funTs);
        wasVolFlag = 1;
    else
        wasVolFlag = 0;
    end
    vec = funTs.vec;
else
    vec = funTs;
end

%% Define polynomial regressors
X = [];
X(:,end+1) = ones(size(t)); X(:,end) = X(:,end)./norm(X(:,end));
X(:,end+1) = t; X(:,end) = X(:,end) - mean(X(:,end)); X(:,end) = X(:,end)./norm(X(:,end));
X(:,end+1) = t.^2; X(:,end) = X(:,end) - mean(X(:,end)); X(:,end) = X(:,end)./norm(X(:,end));

%% Fit
beta = vec'/X';
% figure('WindowStyle','docked');
% plot(squeeze(t),mean(vec,2)); hold on
% plot(squeeze(t),mean(X*beta',2)); hold off
% ind = randperm(size(vec,2),1);
% plot(squeeze(t),vec(:,ind)); hold on
% plot(squeeze(t),X*beta(ind,:)'); hold off


%% Remove from data
vec = vec - X*beta';

%% Percent BOLD
if prcBOLDflag
    vec = vec ./ (X(:,1)*beta(:,1)');
end

%% Output
if isstruct(funTs)
    funTs.vec = vec;
    % funTs.vecMean = X(1,1)*beta(:,1)';
    if wasVolFlag
        funTs = vec2vol(funTs);
    end
else
    funTs = vec;
    % funTsMean = X(1,1)*beta(:,1)';
end
funTsMean = nan(size(funTs.vol2vec));
funTsMean(funTs.vol2vec) = X(1,1)*beta(:,1)';
