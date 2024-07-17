function funPsd = runPSD4(funTs,W,K,mask,runFlag)
%Wrapper of the Chronux's mtspectrumc function for multitaper estimation of
%pds spectra, compatible with MRI data imported by MRIread.m.
%    Parameterization is simplified to the halfbandwidth parameter W only.
%    Includes minimal detrending to avoid distortions from the low frquency
%    edge.
%    See funPsd.psd for other useful parameters.
%    funPsd.tr reflects the frequency resolution in Hz*1000
if ~exist('W','var'); W = []; end
if ~exist('K','var'); K = []; end
if isempty(K) && isempty(W); K = 1; end
if ~exist('mask','var'); mask = []; end
if ~exist('runFlag','var')
    if length(funTs)>1
        runFlag = 'md';
    else
        runFlag = '';
    end
end


if ~strcmp(runFlag,'md')
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
else
% if strcmp(runFlag,'md')
    %% Normalize each run to a noise floor of 1/Hz
    thermalNoiseRange = [0.5 inf];
    modeToRemove = 1:10;
    funTsNorm = cell(size(funTs));
    for runInd = 1:length(funTs)
        % funPsd = runPSD4(funTs(runInd));
        
        svdProulx = runMTsvdProulx(funTs(runInd),thermalNoiseRange);
        % [nnz(svdProulx.vol2vec); numel(svdProulx.vol2vec)]
        recTs = reconMTsvdProulx_spaceXtime(svdProulx,funTs(runInd),modeToRemove,[]);
        recPsd = runPSD4(recTs);
        funTsNorm{runInd} = normPSD2(recPsd,funTs(runInd),thermalNoiseRange);

        % modeTs = reconMTsvdProulx_modeTs(svdProulx,funTs(runInd));
        % funPsd = runPSD4(funTs(runInd));
        % modePsd = reconMTsvdProulx_modePsd(svdProulx,funPsd);
        % modeOrder = 1;
        % viewMTsvdProulx(svdProulx,modePsd,modeTs,modeOrder,funPsd,funTs(runInd));
        % viewPSD3(funPsd)
        % viewPSD3(recPsd)
        
        % viewPSD3(funPsd,1)
        % recTs = reconMTsvdProulx_spaceXtime(svdProulx,funTs(runInd),1:5,[]);
        % recPsd = runPSD4(recTs);
        % viewPSD3(recPsd,1)
        % recPsdNorm = normPSD2(recPsd,[],[0.5 inf]);
        % viewPSD3(recPsdNorm,1)
        % funPsdNorm = normPSD2(recPsd,funPsd,[0.5 inf]);
        % viewPSD3(funPsdNorm,1)
        % funTsNormPsd = runPSD4(funTs(runInd));
        % viewPSD3(funTsNormPsd,1)
    end
    funTsNorm = [funTsNorm{:}];
    % viewPSD3(runPSD4(funTs(4)))
    % viewPSD3(runPSD4(funTsNorm(4)))
    % normPSD2(recPsd,runPSD4(funTs(1)),thermalNoiseRange)
    funTs = funTsNorm;% clear funTsNorm
end

%% Put runs in same structure
if length(funTs)>1
    funTs2 = funTs(1);
    if isfield(funTs,'normFac')
        funTs2.normFac = cat(6,funTs.normFac);
        funTs2.vol2vec = all(~isnan(funTs2.normFac),6);
    end

    % %%% Remove global mean from each run
    % if strcmp(runFlag,'md')
    %     for runInd = 1:length(funTs)
    %         tmp = permute(funTs(runInd).vol,[4 1 2 3]);
    %         funTs(runInd).vol = funTs(runInd).vol - permute(mean(tmp(:,funTs2.vol2vec),2),[2 3 4 1]);
    %         % funTsNorm(runInd).vol = funTsNorm(runInd).vol - permute(mean(tmp(:,funTs2.vol2vec),2),[2 3 4 1]);
    %     end
    % end
    
    funTs2.vol = cat(6,funTs.vol);

    if isfield(funTs,'imMean')
        funTs2.imMean = cat(6,funTs.imMean);
    else
        funTs2.imMean = mean(funTs2.vol,4);
    end

    funTs2.nruns = length(funTs);
    if isfield(funTs,'acqTime')
        funTs2.acqTime = cat(6,funTs.acqTime);
    end

    funTs = funTs2; clear funTs2
end
if strcmp(runFlag,'md')
    if ~isfield(funTs,'nruns')
        funTs.nruns = 1;
    end
    tr = funTs.tr/1000;
    acqTimeSec = datevec(funTs.acqTime); acqTimeSec = etime(acqTimeSec,acqTimeSec(1,:));
    acqTimeSec = round(acqTimeSec/tr)*tr;
    acqTimeSec = permute(acqTimeSec,[2 3 4 5 6 1]);
    funTs.t = repmat(permute((0:funTs.nframes-1)*tr,[1 3 4 2]),[1 1 1 1 1 funTs.nruns]);
    funTs.t = funTs.t + acqTimeSec;
    funTs.vol = funTs.vol(:,:,:,:);
    funTs.t = funTs.t(:,:,:,:);
    % funTs.nframes = size(funTs.vol,4);
    funTs.acqTime = acqTimeSec;

    %%% Detrend the full session
    [funTs,~] = dtrnd4psd(funTs);
    % funTs.vol = funTs.vol - mean(funTs.vol,4);
end
% funTs = vol2vec(funTs);
% plot(squeeze(funTs.t),squeeze(mean(funTs.vec,2)))
%% Compute PSD
funPsd = doIt(funTs,W,K,mask);

%% Add image mean
funPsd.imMean = funTs.imMean;

%% Split back in different structures
if size(funPsd.vec,4)>1
    funPsd2 = funPsd;
    funPsd2.vec = [];
    funPsd2.nruns = 1;
    funPsd2 = repmat(funPsd2,[1 size(funPsd.vec,4)]);
    for i = 1:size(funPsd.vec,4)
        funPsd2(i).vec = funPsd.vec(:,:,:,i);
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
if ~isfield(funTs,'nruns')
    funTs.nruns = 1;
end
if Wflag && Kflag
    error('Cannot specify both W and K');
elseif Wflag
    if isfield(funTs,'t')
        dbstack; error('code that')
    end
    T = tr.*funTs.nframes*funTs.nruns;
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
    T = tr.*funTs.nframes*funTs.nruns;
    W = TW/T; Wreal = W;
    param.tapers = [TW K];
    display(['k  (number of tapers) requested  : ' num2str(K)])
    display(['w  (halfwidth) used       : ' num2str(W,'%0.5f ')])
    display(['tw (time-halfwidth) used  : ' num2str(TW)])
end


%% Perform the multitaper PSD estimation
funPsd = funTs; funPsd.vec = [];
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
    funPsd.nframes = size(funPsd.vec,1);
else
    if ~isfield(funTs,'t')
        funTs.t = [];
    end
    % Compute spectra, using missing data Slepian tapers only if necessary
    [~,funPsd.vec,f] = mtspectrumc3(funTs.vec, param,funTs.t);
    
    if ~isfield(funPsd,'nruns') || isempty(funPsd.nruns); funPsd.nruns = 1; end
    funPsd.nframes = size(funPsd.vec,1)/funPsd.nruns;
    % figure('WindowStyle','docked');
    % plot(f,mean(mean(conj(funPsd.vec).*funPsd.vec,3),2))
    % [~,funPsd.vec,f] = mtspectrumc2(funTs.vec, param);
end
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
% funPsd.psd.maskLabel = funPsd.vol2vecFlag;
funPsd.psd.param = param;
