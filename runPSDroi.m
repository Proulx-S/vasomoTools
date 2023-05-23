function psdRoi = runPSDroi(funTs,W,mask,maskLabel,vecNorm)
% Wrapper of the Chronux's mtspectrumc function for multitaper estimation
% of pds spectra, compatible with MRI data imported by MRIread.m.
%   Parameterization is simplified to the halfbandwidth parameter W only.
%   Includes minimal detrending to avoid distortions from the low frquency
%   edge.
%   funTs is optionally normalized by dividing by the square of vecNorm,
%   which should be obtained from previously computed voxel-wise psd.
%   psd are averaged within mask.
if ~exist('maskLabel','var')
    maskLabel = 'custom';
end

%% Apply normalization and mask
if isfield(funTs,'vol') && ~isempty(funTs.vol)
    if isfield(funTs,'vol2vec') && ~isempty(funTs.vol2vec)
        error('code that')
    else
        if exist('vecNorm','var')
            tmp = vol2vec(funTs);
            volNorm = nan(size(tmp.vol2vec));
            volNorm(tmp.vol2vec) = vecNorm;
        else
            error('code that')
        end
    end
else
    error('code that')
end
vecNorm = volNorm(logical(mask))';
funTs = vol2vec(funTs,mask,1);
if any(all(funTs.vec==0,1))
    warning('Some voxels are all 0s. Adjust your mask to avoid later problems')
end
funTs.vec = funTs.vec ./ sqrt(vecNorm);


%% Set parameters
tr = funTs.tr/1000;
nFrame = funTs.nframes;

T = tr.*nFrame;
TW = T*W;
K = round(TW*2-1);
TW = (K+1)/2;
param.tapers = [TW K];

%% Display actual half-widht used
Wreal = TW/T;
display(['w  (halfwidth) requested  : ' num2str(W,'%0.5f ')])
display(['w  (halfwidth) used       : ' num2str(Wreal,'%0.5f ')])
display(['tw (time-halfwidth) used  : ' num2str(TW)])
display(['k  (number of tapers) used: ' num2str(K)])

%% Detrend time series (detrend up to order-2 polynomial, since this is the highest order not fitting a sinwave)
funTs = dtrnd4psd(funTs);

%% Perform the multitaper PSD estimation
param.Fs = 1/tr;
param.err = [1 0.05];
param.trialave = 1;
[psd,f,psdErr] = mtspectrumc(funTs.vec, param);

%% Output some stuff
psdRoi.dim = strjoin({'freq' 'errBound' 'roi'},' X ');
psdRoi.psd = psd;
psdRoi.psdErr = psdErr';
psdRoi.f = f';
psdRoi.w = Wreal;
psdRoi.tw = TW;
psdRoi.param = param;
psdRoi.mask = mask;
psdRoi.maskLabel = maskLabel;
