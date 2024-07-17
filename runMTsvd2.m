function svdStruct = runMTsvd2(anaType,funTs,fpass,W,K,mask,vecNorm)
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
% funTs.vec = funTs.vec - mean(funTs.vec,2);

%% Set multitaper parameters
tr = funTs.tr/1000;
nFrame = funTs.nframes;
param.tapers = [];
param.Fs = 1/tr;
switch anaType
    case 'svdKlein'
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
        if ~isempty(fpass)
            param.fpass = fpass;
        end
        mdkp = [];
        [~,f] = mtspectrumc(funTs.vec(:,1), param);

    case 'svdMitra'
        %%% Set parameters for the user-defined frequency band
        W = diff(fpass)/2;
        T = tr.*nFrame;
        TW = T*W;
        K = round(TW*2-1);
        TW = (K+1)/2;
        param.tapers = [TW K];
        [~,f] = mtspectrumc(funTs.vec(:,1), param);
        f0 = fpass(1)+W; [~,b] = min(abs(f - f0)); f0 = f(b);
        param.fpass = [f0 f0];
        mdkp = [];
        %%% Display actual frequency band used
        Wreal = W;
        fpassReal = f0+[-1 1].*(TW/T);
        display(['frequency band requested: fpass=[' num2str(fpass,'%0.5f ') ']'])
        display(['frequency band used     : fpass=[' num2str(fpassReal,'%0.5f ') ']'])
    case 'svdProulx'
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
    otherwise
        error('Invalid anaType. Choose one of ''svdMitra'', ''svdKlein'' or ''svdProulx''')
end

%% Run the decomposition
tic
switch anaType
    case 'svdProulx'
        [sv,sp,fm,MTS] = spsvd3(funTs.vec,param);
    otherwise
        [sv,sp,fm,tapers,tvec,f,data] = spsvd(funTs.vec,param,mdkp);
end
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

