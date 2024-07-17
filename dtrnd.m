function [funTs,funTsMean] = dtrnd(funTs,t,prcBOLDflag)
%Detrend time series for multitaper estimation of psd
%   Simple detrending using order-2 polynomials, the highest order that
%   can't fitting a sinwave.
if ~exist('prcBOLDflag','var') || isempty(prcBOLDflag)
    prcBOLDflag = 0;
end
if isstruct(funTs)
    if length(funTs)==1
        tr = funTs.tr/1000;
        if isfield(funTs,'t')
            T = funTs.t;
        else
            T = 0:tr:(funTs.nframes-1)*tr;
        end
        if ~isfield(funTs,'vec') || isempty(funTs.vec)
            funTs = vol2vec(funTs);
            wasVolFlag = 1;
        else
            wasVolFlag = 0;
        end
        vec = funTs.vec;

        if isfield(funTs,'nruns')
            nruns = funTs.nruns;
        elseif isfield(funTs,'vec') && ~isempty(funTs.vec)
            nruns = size(funTs.vec,4);
            funTs.nruns = nruns;
        else
            dbstack; error('code that');
        end
        imMean = nan(1,nnz(funTs.vol2vec),1,funTs.nruns);
    elseif length(funTs)>1
        for i = 1:length(funTs)
            [funTsX(i),~] = dtrnd(funTs(i)); funTs(i).vol = []; funTs(i).vec = [];
        end
        funTs = funTsX; clear funTsX
        funTsMean = [];
        return
    end
else
    vec = squeeze(funTs);
    vec = permute(vec(:,:),[1 3 4 2]);
    nruns = size(vec,4);
    sz = size(vec); sz(1) = 1;
    funTsMean = nan(sz);
    if exist('t','var') && ~isempty(t)
        T = t;
    else
        T = repmat((0:size(vec,1)-1)',[1 1 1 1 1 nruns]);
    end
end
if isduration(T)
    T = seconds(T);
end


for runInd = 1:nruns
    %% Define polynomial regressors
    t = T(:,:,:,:,:,runInd);
    X = [];
    X(:,end+1) = ones(size(t)); X(:,end) = X(:,end)./norm(X(:,end));
    X(:,end+1) = t; X(:,end) = X(:,end) - mean(X(:,end)); X(:,end) = X(:,end)./norm(X(:,end));
    X(:,end+1) = t.^2; X(:,end) = X(:,end) - mean(X(:,end)); X(:,end) = X(:,end)./norm(X(:,end));

    %% Fit
    beta = permute(vec(:,:,:,runInd),[2 1])/X';
    % figure('WindowStyle','docked');
    % plot(squeeze(t),mean(vec,2)); hold on
    % plot(squeeze(t),mean(X*beta',2)); hold off
    % ind = randperm(size(vec,2),1);
    % plot(squeeze(t),vec(:,ind)); hold on
    % plot(squeeze(t),X*beta(ind,:)'); hold off

    %% Remove from data
    vec(:,:,:,runInd) = vec(:,:,:,runInd) - X*permute(beta,[2 1]);

    %% Percent BOLD
    if prcBOLDflag
        vec(:,:,:,runInd) = vec(:,:,:,runInd) ./ (X(:,1)*permute(beta(:,1),[2 1]));
    end

    %% Image mean
    imMean(:,:,:,runInd) = X(1,1)*beta(:,1)';
end

%% Output
if isstruct(funTs)
    funTs.vec = vec;
    % funTs.vecMean = X(1,1)*beta(:,1)';
    if wasVolFlag
        funTs = vec2vol(funTs);
    end
    funTsMean = nan([1 1 funTs.nruns size(funTs.vol2vec,1:3)]);
    funTsMean(:,:,:,funTs.vol2vec) = permute(permute(imMean,[5 6 4 1 2 3]),[4 5 6 1 2 3]);
    funTs.imDtrndMean = funTsMean;
    % funTs.imDtrndMean = nan(size(funTs.vol2vec));
    % funTs.imDtrndMean(funTs.vol2vec) = X(1,1)*beta(:,1)';
else
    funTs = vec;
    for runInd = 1:nruns
        funTsMean(:,:,:,runInd) = X(1,1)*beta(:,1)';
    end
    % funTsMean = X(1,1)*beta(:,1)';
end
% funTsMean = nan([funTs.nruns size(funTs.vol2vec,1:3) 1]);
% funTsMean(:,funTs.vol2vec) = permute(imMean,[4 2 1 3]); %X(1,1)*beta(:,1)';
% funTsMean = permute(funTsMean,[2 3 4 5 6 1]);
