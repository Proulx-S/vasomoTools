function funPsd = runPSD3(funTs,W,K,mask)
%Wrapper of the Chronux's mtspectrumc function for multitaper estimation of
%pds spectra, compatible with MRI data imported by MRIread.m.
%    Parameterization is simplified to the halfbandwidth parameter W only.
%    Includes minimal detrending to avoid distortions from the low frquency
%    edge.
%    See funPsd.psd for other useful parameters.
%    funPsd.tr reflects the frequency resolution in Hz*1000
if ~exist('W','var'); W = []; end
if ~exist('K','var'); K = []; end
if ~exist('mask','var'); mask = []; end

%% Detrend time series (detrend up to order-2 polynomial, since this is the highest order not fitting a sinwave)
tic
funTs2 = cell(size(funTs));
for i = 1:length(funTs)
    [funTs2{i},imMean] = dtrnd4psd(funTs(i));
    if ~isfield(funTs2{i},'imMean') || isempty(funTs2{i}.imMean)
        funTs2{i}.imMean = imMean; clear imMean
    end
end
funTs = [funTs2{:}];
toc

%% Put runs in same structure
if length(funTs)>1
    funTs2 = funTs(1);
    funTs2.vol = cat(6,funTs.vol);
    funTs2.imMean = cat(6,funTs.imMean);
    funTs2.nruns = length(funTs);
    funTs = funTs2;
end

%% Compute PSD
funPsd = doIt(funTs,W,K,mask);

%% Add image mean
funPsd.imMean = funTs.imMean;

%% Split back in different structures
nruns = size(funPsd.vec,4);
if nruns>1
    funPsd2 = funPsd;
    funPsd2.vec = [];
    funPsd2.nruns = 1;
    funPsd2 = repmat(funPsd2,[1 nruns]);
    for i = 1:nruns
        funPsd2(i).vec = funPsd.vec(:,:,:,i);
        % funPsd2(i).vecMean = funPsd.vecMean(:,:,:,i);
        funPsd2(i).tr = funPsd.tr;
    end
    funPsd = funPsd2; clear funPsd2
end

function funPsd = doIt(funTs,W,K,mask)
if exist('mask','var') && ~isempty(mask)
    funTs = vol2vec(funTs,mask,1);
else
    funTs = vol2vec(funTs);
end
tr = funTs.tr/1000;

if any(all(funTs.vec==0,1))
    warning('Some voxels are all 0s. Adjust your mask to avoid later problems')
end

%% Set parameters
Wflag = exist('W','var') && ~isempty(W);
Kflag = exist('K','var') && ~isempty(K);
if Wflag && Kflag
    error('Cannot specify both W and K');
elseif Wflag
    T = tr.*funTs.nframes;
    TW = T*W;
    K = round(TW*2-1);
    TW = (K+1)/2;
    param.tapers = [TW K];
    Wreal = TW/T;
    display(['w  (halfwidth) requested  : ' num2str(W,'%0.5f ')])
    display(['w  (halfwidth) used       : ' num2str(Wreal,'%0.5f ')])
    display(['tw (time-halfwidth) used  : ' num2str(TW)])
    display(['k  (number of tapers) used: ' num2str(K)])
elseif Kflag
    TW = (K+1)/2;
    T = tr.*funTs.nframes;
    W = TW/T; Wreal = W;
    param.tapers = [TW K];
    display(['k  (number of tapers) requested  : ' num2str(K)])
    display(['w  (halfwidth) used       : ' num2str(W,'%0.5f ')])
    display(['tw (time-halfwidth) used  : ' num2str(TW)])
end


%% Perform the multitaper PSD estimation
funPsd = funTs; funPsd.vec = []; funPsd.vecMean = [];
param.Fs = 1/tr;
param.complex = 1;
% [funPsd.vec,f] = mtspectrumc(funTs.vec, param);
if size(funTs.vec,4)>1
    [~,~,f] = mtspectrumc2(funTs.vec(:,1), param);
    sz = size(funTs.vec); sz(1) = length(f); sz(3) = param.tapers(2);
    vec = nan(sz);
    parfor i = 1:size(funTs.vec,4)
        [~,vec(:,:,:,i),~] = mtspectrumc2(funTs.vec(:,:,:,i), param);
    end
    funPsd.vec = vec;
else
    [~,funPsd.vec,f] = mtspectrumc2(funTs.vec, param);
end

% % [funPsd.vec,vecC,f] = mtspectrumc2(funTs.vec(:,:,:,i), param);
% if size(funTs.vec,4)>1
%     vec = repmat(funPsd.vec,[1 1 1 size(funTs.vec,4)]);
%     vecC = repmat(vecC,[1 1 1 size(funTs.vec,4)]);
%     f = repmat(f,[1 1 1 size(funTs.vec,4)]);
%     parfor i = 2:size(funTs.vec,4)
%         [vec(:,:,:,i),vecC(:,:,:,i),f(:,:,:,i)] = mtspectrumc2(funTs.vec(:,:,:,i), param);
%     end
%     funPsd.vec = vec;
% end
funPsd.nframes = size(funPsd.vec,1);
funPsd.ntapers = param.tapers(2);
funPsd.tr = mode(diff(f))*1000;


%% Output some stuff
if isfield(funPsd,'psd')
    funPsd = rmfield(funPsd,'psd');
end
funPsd.psd.dim = strjoin({'space' 'freq' 'taper' 'run'},' X ');
tmp = diff(f,[],4); if any(tmp(:)); dbstack; error('X'); end
funPsd.psd.f = f(:,:,:,1);
funPsd.psd.w = Wreal;
funPsd.psd.tw = TW;
funPsd.psd.mask = funTs.vol2vec;
funPsd.psd.maskLabel = funPsd.vol2vecFlag;
funPsd.psd.param = param;
