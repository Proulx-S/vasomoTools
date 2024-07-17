function [svdStruct,svdTs,svdPsd,funPsd] = MTsvdProulx(funTs,fpass,modeOrder,mask,vecNorm)
% Similar to Mitra 1997. A single svd is run on data tapered for
% sensitivity over user-defined frequency band (fpass).
tsMean = mean(funTs.vol,4);
if ~exist('mask','var') || isempty(mask)
    funTs = vol2vec(funTs);
    mask = funTs.vol2vec;
else
    funTs = vol2vec(funTs,mask,1);
end

%% Apply timeseries normalization
% normFact (vector) based on psd so we need to use its square root here
if exist('normFact','var') && exist('vecNorm','var') && ~isempty(vecNorm)
    funTs.vec = funTs.vec ./ sqrt(vecNorm(logical(mask(funTs.vol2vec))));
end

%% Detrend time series (detrend up to order-2 polynomial, since this is the highest order not fitting a sinwave)
funTs = dtrnd4psd(funTs);

%% Scale variance

%% Set multitaper parameters
tr = funTs.tr/1000;
nFrame = funTs.nframes;
param.tapers = [];
param.Fs = 1/tr;
%%% Set parameters for each user-defined frequency band
fpassOrig = fpass;
paramOrig = param;
for bandInd = 1:size(fpassOrig,1)
    fpass_targ = fpassOrig(bandInd,:);
    W_targ = diff(fpass_targ)/2;
    T = tr.*nFrame;
    TW_targ = T*W_targ;
    K_actual = round(TW_targ*2-1);
    TW_actual = (K_actual+1)/2;
    W_actual = TW_actual/T;

    paramCur = paramOrig;
    paramCur.tapers = [TW_actual K_actual];
    [~,f] = mtspectrumc(funTs.vec(:,1), paramCur);

    f0_targ = mean(fpass_targ); [~,b] = min(abs(f - f0_targ));
    f0_actual = f(b);
    fpass_actual = f0_actual+[-1 1].*W_actual;

    paramCur.fpass = [f0_actual f0_actual];
    mdkp = [];

    param.tapers(bandInd,:) = paramCur.tapers;
    param.fpass(bandInd,:) = paramCur.fpass;
    param.BW(bandInd,1) = W_actual;
end
% display(['frequency band requested: fpass=[' num2str(fpass,'%0.5f ') ']'])
% display(['frequency band used     : fpass=[' num2str(fpassReal,'%0.5f ') ']'])

%% Run the decomposition
tic
[sv,sp,fm,MTS] = spsvd3(funTs.vec,param);
toc

%% Output
svdStruct.mask = mask;
if exist('vecNorm','var')
    svdStruct.normFact = vecNorm;
else
    svdStruct.normFact = [];
end
svdStruct.tsMean = tsMean;
svdStruct.dim = strjoin({'space/taper' 'freq/time' 'modes'},' X ');
svdStruct.sv = permute(sv,[3 1 2]);
svdStruct.sp = sp; %spatial singular vectors
svdStruct.fm = fm; %taper singular vectors
svdStruct.MTS.proj = permute(MTS.proj,[2 1 3]);
svdStruct.MTS.tpInd = permute(MTS.tpInd,[2 1 3]);
svdStruct.MTS.bandInd = permute(MTS.bandInd,[2 1 3]);
svdStruct.MTS.bandFreq = permute(MTS.bandFreq,[1 3 2]);
svdStruct.MTS.bandBw = param.BW';
svdStruct.c = svdStruct.sv(1,:,:,:,:).^2./sum(svdStruct.sv.^2,1);
svdStruct.param = param;
svdStruct.param.info = 'band X :';

%% Reconstruct timeseries of the first component
if nargout>=2
    sp = svdStruct.sp(:,:,modeOrder); % spatial sv (space x 1 x mode)
    fm = svdStruct.fm(:,:,modeOrder); % taper sv (taper x 1 x mode)
    tp = svdStruct.MTS.proj; % tapers (taper x time)
    f = svdStruct.MTS.bandFreq; % frequencies (1 x freq)
    fInd = svdStruct.MTS.bandInd; % frequency indices (taper x 1)
    tsRec = fm'*tp; % WARNING! taking the conjugate here with the ' operator
    svdTs = funTs; svdTs.vol = permute(tsRec,[1 3 4 2]); svdTs.volMean = mean(tsRec);
    svdTs.volsize([1 2]) = 1; svdTs.height = 1; svdTs.width = 1; svdTs.nvoxels = 1;
end
%% Compute psd of the reconstructed timeseries
if nargout>=3
    W = [];
    K = 1;
    cFlag = 1;
    svdPsd = runPSD(svdTs,W,K,cFlag);
end
%% For comparison, compute psd of original timeseries
if nargout>=4
    funPsd = vec2vol(runPSD(funTs,W,K,cFlag));
end

