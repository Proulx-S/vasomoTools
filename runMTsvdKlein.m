function svdStruct = runMTsvdKlein(funTs,fpass,W,K,mask,vecNorm,scaleFlag,M)
if ~exist('scaleFlag','var') || isempty(scaleFlag)
    scaleFlag = 0;
end
if ~exist('M','var') || isempty(M); M=[]; end

% %% Store image mean before processing
% tsMean = nan(size(funTs.vol2vec));
% tsMean(:) = mean(funTs.vec,1);

%% Apply timeseries normalization
% normFact (vector) based on psd so we need to use its square root here
if exist('normFact','var') && exist('vecNorm','var') && ~isempty(vecNorm)
    for i = 1:length(funTs)
        funTs(i).vec = funTs(i).vec ./ sqrt(vecNorm(logical(mask(funTs(i).vol2vec))));
    end
end

%% Detrend time series (detrend up to order-2 polynomial, since this is the highest order not fitting a sinwave)
tic
funTs2 = cell(size(funTs));
parfor i = 1:length(funTs)
    funTs2{i} = dtrnd4psd(funTs(i));
end
funTs = [funTs2{:}];
toc

%% Put runs in same structure
if length(funTs)>1
    funTs = vec2vol(funTs);
    funTs2 = funTs(1);
    funTs2.vol = cat(6,funTs.vol);
    funTs2.volMean = cat(6,funTs.volMean);
    funTs = vol2vec(funTs2); clear funTs2
    funTs.nruns = length(funTs);
end

if exist('mask','var') && ~isempty(mask)
    funTs = vol2vec(vec2vol(funTs),mask,1);
else
    mask = [];
end

%% Compute SVD
svdStruct = doIt(funTs,fpass,W,K,mask,vecNorm,scaleFlag,M);


function svdStruct = doIt(funTs,fpass,W,K,mask,vecNorm,scaleFlag,M)
funTs = vol2vec(funTs);
%% Set multitaper parameters
tr = funTs.tr/1000;
nFrame = funTs.nframes;
param.tapers = [];
param.Fs = 1/tr;
param.augFac = size(funTs.vec,4);

Wflag = exist('W','var') && ~isempty(W);
Kflag = exist('K','var') && ~isempty(K);
if Wflag && Kflag
    error('Cannot specify both W and K');
elseif Wflag
    display(['w  (halfwidth) requested  : ' num2str(W,'%0.5f ')])
    T = tr.*funTs.nframes;
    TW = T*W;
    K = round(TW*2-1);
    if K<=1; warning('K must be minimum 2, changing it'); K = 2; end
    TW = (K+1)/2;
    W = TW/T; Wreal = W;
    param.tapers = [TW K];
    display(['w  (halfwidth) used       : ' num2str(Wreal,'%0.5f ')])
    display(['tw (time-halfwidth) used  : ' num2str(TW)])
    display(['k  (number of tapers) used: ' num2str(K)])
elseif Kflag
    display(['k  (number of tapers) requested  : ' num2str(K)])
    if K<=1; warning('K must be minimum 2, changing it'); K = 2; end
    TW = (K+1)/2;
    T = tr.*funTs.nframes;
    W = TW/T; Wreal = W;
    param.tapers = [TW K];
    display(['w  (halfwidth) used       : ' num2str(W,'%0.5f ')])
    display(['tw (time-halfwidth) used  : ' num2str(TW)])
end
if ~isempty(fpass)
    param.fpass = fpass;
end
mdkp = [];
[~,f] = mtspectrumc(funTs.vec(:,1), param);

%% Run the decomposition
tic
param.scaleFlag = scaleFlag;
if size(funTs.vec,4)>1
    [sv,sp,fm,tapers,tvec,f] = spsvd_runAugmtd(funTs.vec,param,mdkp,M);
else
    [sv,sp,fm,tapers,tvec,f] = spsvd(funTs.vec,param,mdkp);
end
toc

%% Output
if ~exist('mask','var') || isempty(mask)
    svdStruct.mask = funTs.vol2vec;
else
    svdStruct.mask = mask;
end
if exist('vecNorm','var')
    svdStruct.normFact = vecNorm;
else
    svdStruct.normFact = [];
end
if ~isfield(funTs,'imMean') || isempty(funTs.imMean)
    if isempty(funTs.vol)
        funTs.imMean = nan(size(funTs.vol2vec));
        funTs.imMean(funTs.vol2vec) = mean(funTs.vec,1);
    else
        funTs.imMean = mean(funTs.vol(:,:,:,:),4);
    end
end
svdStruct.imMean = funTs.imMean;
svdStruct.dim = strjoin({'space/taper' 'freq/time' 'modes'},' X ');
svdStruct.sv = permute(sv,[3 1 2]);
svdStruct.sp = sp; %spatial singular vectors
svdStruct.fm = fm; %taper singular vectors
svdStruct.tp = permute(tapers,[2 1]); %tapers
svdStruct.tv = permute(tvec(:,1),[2 1]); %time vector
svdStruct.f = f; %freq vector
svdStruct.tf = 'tp.*exp(-f(i)*tv)'; %taper frequency time vector
svdStruct.c = sv(:,1)'.^2./sum(sv.^2,2)';
svdStruct.w = Wreal;
svdStruct.param = param;

% figure('WindowStyle','docked');
% plot(svdStruct.f,svdStruct.c)