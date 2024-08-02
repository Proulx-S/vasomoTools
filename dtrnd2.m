function [funTs,poly,order] = dtrnd2(funTs,t,prcBOLDflag,order)
%Detrend time series for multitaper estimation of psd
%   Simple detrending using order-2 polynomials, the highest order that
%   can't fitting a sinwave.
%   Runs in different struct are processed seperatly
%   Runs in the same struct are processed together, as a single timeseries
%   with missing data
%   When funTs is a scalar, t must be provided
if ~exist('prcBOLDflag','var'); prcBOLDflag = []; end
if isempty(prcBOLDflag);        prcBOLDflag = 0; end
if ~exist('t','var'); t = []; end
if ~exist('order','var'); order = []; end

if iscell(funTs)
    for I = 1:numel(funTs)
        disp(['detrending cell ' num2str(I) '/' num2str(numel(funTs))])
        if ~strcmp(mfilename,'dtrnd2'); dbstack; error('double-check that'); end
        funTs{I} = dtrnd2(funTs{I},t,prcBOLDflag,order);
    end
    poly = [];


elseif isstruct(funTs)
    funTs = vol2vec(funTs);
    for I = 1:numel(funTs)
        if ~strcmp(mfilename,'dtrnd2'); dbstack; error('change function name (dtrnd2)'); end
        disp(['detrending struct ' num2str(I) '/' num2str(numel(funTs))])
        if ~isfield(funTs(I),'t') || isempty(funTs(I).t)
            tr = funTs(I).tr/1000;
            funTs(I).t = (0:tr:(funTs(I).nframes-1)*tr)';
        end
        [funTs(I).vec,funTs(I).poly] = dtrnd2(funTs(I).vec,funTs(I).t,prcBOLDflag,order);
        funTs(I).poly.vol2vec = funTs(I).vol2vec;
    end
    order = size(funTs(I).poly.beta,2)-1;
    funTs = setNiceFieldOrder(funTs,{'vol' 'vol2vec' 'vec' 't' 'poly' 'volInfo' 'vecInfo'});
    poly = [];



elseif isnumeric(funTs)
    %% Define polynomial regressors
    if length(t)==1
        tr = t;
        t = 0:tr:(size(funTs,1)-1)*tr;
    end
    X = [];
    if isempty(order)
        X(:,end+1) = ones(size(t(:))); X(:,end) = X(:,end)./norm(X(:,end));
        X(:,end+1) = t(:); X(:,end) = X(:,end) - mean(X(:,end)); X(:,end) = X(:,end)./norm(X(:,end));
        X(:,end+1) = t(:).^2; X(:,end) = X(:,end) - mean(X(:,end)); X(:,end) = X(:,end)./norm(X(:,end));
    elseif order~=-1
        X(:,end+1) = ones(size(t(:))); X(:,end) = X(:,end)./norm(X(:,end));
        for i = 1:order
            X(:,end+1) = t(:).^i; X(:,end) = X(:,end) - mean(X(:,end)); X(:,end) = X(:,end)./norm(X(:,end));
        end
    else
        poly = [];
        return
    end

    %% Fit
    sz = size(funTs,[1 2 3 4]);
    beta = reshape(permute(funTs,[2 1 4 3]),[sz(2) prod(sz([1 3 4]))]) / X';

    %% Remove from data
    funTs = funTs - permute(reshape(permute(X*permute(beta,[2 1]),[2 1]),sz([2 1 3 4])),[2 1 3 4]);
    
    beta0 = X(1,1)*permute(beta(:,1),[2 1]);
    %% Percent BOLD
    if prcBOLDflag
        funTs = funTs./beta0.*100;
    end

    %% Outputs
    poly.beta = beta; clear beta
    poly.t = permute(t(:),[2 3 1]); clear t
    poly.beta0 = permute(beta0,[2 1 3 4]); clear beta0
    poly.info = 'vox x polyOrder x time';
    poly.info2 = 'X(:,end+1) = t(:).^2; X(:,end) = X(:,end) - mean(X(:,end)); X(:,end) = X(:,end)./norm(X(:,end));';
    
       

else
    dbstack; error('this should not happen')
end
