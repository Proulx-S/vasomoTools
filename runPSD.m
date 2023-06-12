function funPsd = runPSD(funTs,W)
%Wrapper of the Chronux's mtspectrumc function for multitaper estimation of
%pds spectra, compatible with MRI data imported by MRIread.m.
%    Parameterization is simplified to the halfbandwidth parameter W only.
%    Includes minimal detrending to avoid distortions from the low frquency
%    edge.
%    See funPsd.psd for other useful parameters.
%    funPsd.tr reflects the frequency resolution in Hz*1000
tMean = mean(funTs.vol,4);
funTs = vol2vec(funTs);
tr = funTs.tr/1000;

if any(all(funTs.vec==0,1))
    warning('Some voxels are all 0s. Adjust your mask to avoid later problems')
end

%% Set parameters
T = tr.*funTs.nframes;
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
funTs.vec = funTs.vec - mean(funTs.vec,1);
% funTs = dtrnd4psd(funTs);

%% Perform the multitaper PSD estimation
funPsd = funTs; funPsd.vec = [];
param.Fs = 1/tr;
[funPsd.vec,f] = mtspectrumc(funTs.vec, param);

funPsd.nframes = size(funPsd.vec,1);
funPsd.tr = mode(diff(f))*1000;

%% Output some stuff
if isfield(funPsd,'psd')
    funPsd = rmfield(funPsd,'psd');
end
funPsd.psd.dim = strjoin({'space' 'freq'},' X ');
funPsd.psd.f = f;
funPsd.psd.w = Wreal;
funPsd.psd.tw = TW;
funPsd.psd.mask = funTs.vol2vec;
funPsd.psd.maskLabel = funPsd.vol2vecFlag;
funPsd.psd.param = param;

funPsd.tMean = tMean;

% if nargout>1
%     psdStruct = funPsd.psd;
%     psdStruct.psd = funPsd.vec';
%     psdStruct = setNiceFieldOrder(psdStruct,{'dim' 'psd' 'f'});
% end

