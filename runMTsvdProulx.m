function svdStruct = runMTsvdProulx(funTs,fpass,mask,scaleFlag)
% Similar to Mitra 1997. A single svd is run on data tapered for
% sensitivity over user-defined frequency band (fpass).
if ~exist('scaleFlag','var') || isempty(scaleFlag); scaleFlag = 0; end
if ~isfield(funTs,'imMean') || isempty(funTs.imMean)
    imMean = cat(4,funTs.vol);
    imMean = mean(imMean(:,:,:,:),4);
else
    imMean = funTs.imMean;
end

% %% Detrend time series (detrend up to order-2 polynomial, since this is the highest order not fitting a sinwave)
% tic
% funTs2 = cell(size(funTs));
% parfor i = 1:length(funTs)
%     funTs2{i} = dtrnd4psd(funTs(i));
% end
% funTs = [funTs2{:}];
% toc

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
svdStruct = doIt(funTs,fpass,mask,scaleFlag);
svdStruct.imMean = imMean;
svdStruct.param.type = 'proulx';


function svdStruct = doIt(funTs,fpass,mask,scaleFlag)
funTs = vol2vec(funTs,mask);

%% Remove mean
for i = 1:length(funTs)
    funTs(i).vec = funTs(i).vec - mean(funTs(i).vec,1);
end

%% Set multitaper parameters
tr = funTs.tr/1000;
nFrame = funTs.nframes;
param.tapers = [];
param.Fs = 1/tr;
paramTmp = param; paramTmp.tapers = [0.5 1];
[~,f] = mtspectrumc(funTs.vec(:,1),paramTmp);
fpass(isinf(fpass)) = f(end);

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
    
    f0_targ = mean(fpass_targ); [~,b] = min(abs(f - f0_targ));
    f0_actual = f(b);
    fpass_actual = f0_actual+[-1 1].*W_actual;

    paramCur.fpass = [f0_actual f0_actual];
    mdkp = [];

    param.tapers(bandInd,:) = paramCur.tapers;
    param.fpass(bandInd,:) = paramCur.fpass;
    param.BW(bandInd,1) = W_actual;
end

%% Run the decomposition
tic
param.scale = scaleFlag;
[sp,sv,fm,MTS] = spsvd_freqAugmntd_runAugmntd(funTs.vec,param);
toc
%%% Scale tapers by Fs/2, because for some reason that is what brings the
%%% reconstruction closest to the original data. Residual reconstruction
%%% error concentrates at the very beginning and end of timeseries, as
%%% expected.
MTS = rmfield(MTS,'projScale');
MTS.proj = MTS.proj./(param.Fs/2);

% sp = squeeze(sp);
% sv = diag(squeeze(sv));
% fm = squeeze(fm);
% A = sp*sv*fm';
% data = real( MTS.proj*A' );
% whos sp sv fm A data
% 
% A2 = u*s*v';
% whos u s v A
% data2 = real( MTS.proj*A2' );
% 
% figure('WindowStyle','docked');
% tmp1 = data(20:end-19,5:end-4);
% tmp2 = data2(20:end-19,5:end-4);
% scatter(tmp1(:),tmp2(:))
% drawnow
% 
% figure('WindowStyle','docked');
% tmp1 = funTs.vec(20:end-19,5:end-4);
% tmp2 = data(20:end-19,5:end-4);
% scatter(tmp1(:),tmp2(:))
% drawnow
% 
% %% Output
% svdStruct.mask = funTs.vol2vec;
svdStruct.vol2vec = funTs.vol2vec;
svdStruct.dim = strjoin({'space/taper' 'freq/time' 'modes'},' X ');
svdStruct.sv = sv;
svdStruct.sp = sp; %spatial singular vectors
svdStruct.fm = fm; %taper singular vectors
svdStruct.MTS.proj = permute(MTS.proj,[2 1 3]);
svdStruct.MTS.tpInd = permute(MTS.tpInd,[2 1 3]);
svdStruct.MTS.bandInd = permute(MTS.bandInd,[2 1 3]);
svdStruct.MTS.bandFreq = permute(MTS.bandFreq,[2 1]);
svdStruct.MTS.bandBw = param.BW;
svdStruct.MTS.K = MTS.K;
svdStruct.MTS.NRUN = MTS.NRUN;
svdStruct.MTS.N = MTS.N;
svdStruct.MTS.NCHAN = MTS.NCHAN;
svdStruct.MTS.info = 'taper/freq x time';
svdStruct.c = svdStruct.sv(1,:,:,:,:).^2./sum(svdStruct.sv.^2,1);
svdStruct.param = param;
svdStruct.param.info = 'band X :';

%% Test reconstruction accuracy
% sp = squeeze(svdStruct.sp);
% sv = diag(squeeze(svdStruct.sv));
% fm = squeeze(svdStruct.fm);
% A = sp*sv*fm';
% data = real( permute(svdStruct.MTS.proj,[2 1 3])*A' );
% randInd = randperm(numel(data),1000);
% % figure('WindowStyle','docked');
% scatter(data(randInd),funTs.vec(randInd))
% grid on
% drawnow