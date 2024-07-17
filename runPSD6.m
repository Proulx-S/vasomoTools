function funPsd = runPSD6(funTs,W,K,mask,memFlag,skipSVD,skipPSD,verbose)
% Wrapper for the Chronux's mtspectrumc function for multitaper estimation of
% pds spectra, compatible with MRI data imported by MRIread.m.
%    Other features are implemented on top or adapted from of Chronux:
%    -missing data slepian tapers catenate consecutive runs
%    -...
%
%    Parameterization is simplified to use either the halfbandwidth parameter W or the number of tapers.
%    Alternatively, precomputed tapers can be input as K to skip to save
%    time (useful when missing data slepian are used).
%    See funPsd.psd for other useful parameters.
%    funPsd.tr reflects the frequency resolution in Hz*1000
if ~exist('verbose','var') || isempty(verbose); verbose = true; end
if ~exist('memFlag','var') || isempty(memFlag); memFlag = false; end
if ~exist('W','var'); W = []; end
if ~exist('K','var'); K = []; end
if isempty(K) && isempty(W); K = 1; end
if length(K)>1
    tp = K; K = size(tp,2);
else
    tp = [];
end
if ~exist('skipSVD','var') || isempty(skipSVD); skipSVD = false; end
if funTs(1).nvoxels==1
    if verbose; disp('only one timeseries, skipping SVD'); end
    skipSVD = true;
end
if ~exist('skipPSD','var') || isempty(skipPSD); skipPSD = false; end
if ~exist('mask','var'); mask = []; end

if exist('mask','var') && ~isempty(mask)
    funTs = applyMask(funTs,mask);
end
tr = funTs.tr/1000;

if any(all(funTs.vec==0,1))
    warning('Some voxels are all 0s. Adjust your mask to avoid later problems')
end

%% Set parameters
Wflag = exist('W','var') && ~isempty(W);
Kflag = exist('K','var') && ~isempty(K);
tpFlag = ~isempty(tp); if tpFlag; Wflag = false; Kflag = false; end
if ~isfield(funTs,'nruns')
    funTs.nruns = 1;
end
T = tr.*funTs.nframes*funTs.nruns;
if Wflag && Kflag
    error('Cannot specify both W and K');
elseif Wflag
    Wreq = W;
    TW = T*W;
    K = round(2*TW-1); if K==0; K=1; end
    TW = (K+1)/2;
    W = TW/T;
    param.tapers = [TW K];

    if verbose
        display(['W  (halfwidth)       : ' num2str(Wreq,'%0.5f ')])
        display(['W  (halfwidth)       : ' num2str(W,'%0.5f ')])
        display(['TW (time-halfwidth)  : ' num2str(TW)])
        display(['K  (number of tapers): ' num2str(K)])
    end
elseif Kflag
    TW = (K+1)/2;
    T = tr.*funTs.nframes*funTs.nruns;
    W = TW/T;
    param.tapers = [TW K];
    if verbose
        display(['K  (number of tapers): ' num2str(K)])
        display(['W  (halfwidth)       : ' num2str(W,'%0.5f ')])
        display(['TW (time-halfwidth)  : ' num2str(TW)])
    end
elseif tpFlag
    TW = (K+1)/2;
    T = tr.*funTs.nframes*funTs.nruns;
    W = TW/T;
    param.tapers = [TW K];
    if verbose
        display(['using precomputed tapers']);
        display(['K  (number of tapers): ' num2str(K)])
        display(['W  (halfwidth)       : ' num2str(W,'%0.5f ')])
        display(['TW (time-halfwidth)  : ' num2str(TW)])
    end
end

funPsd = funTs;
[funPsd.vol] = deal([]);
[funPsd.vec] = deal([]);
tmp = strsplit(funTs(1).vecInfo,' x '); tmp{1} = 'freq/time'; tmp = strjoin(tmp,' x ');
[funPsd.vecInfo] = deal(tmp);
tmp = strsplit(funTs(1).volInfo,' x '); tmp{4} = 'freq/time'; tmp = strjoin(tmp,' x ');
[funPsd.volInfo] = deal(tmp);
param.Fs = 1/tr;
param.complex = 1;
for sInd = 1:length(funTs)
    if ~isfield(funTs(sInd),'t')
        funTs(sInd).t = [];
    end
    
    if ~skipPSD
        if verbose; disp('for psd'); end
        [~,f,funPsd(sInd).tp] = mtspectrumc4(funTs(sInd).vec(:,1), param,funTs(sInd).t,tp);
        funPsd(sInd).tp = permute(funPsd(sInd).tp,[3 4 5 1 2]);
    else
        dbstack; error('code that')
    end
    
    if ~skipSVD
        if verbose; disp('for svd'); end
        paramSVD = param;
        paramSVD.tapers(2) = paramSVD.tapers(2)*2;
        paramSVD.tapers(1) = (paramSVD.tapers(2)+1)/2;
        [~,~,funPsd(sInd).svd.tp] = mtspectrumc4(funTs(sInd).vec(:,1), paramSVD,funTs(sInd).t,tp);
        funPsd(sInd).svd.tp = permute(funPsd(sInd).svd.tp,[2 3 1]);
    end

    if ~isempty(tp); tp = []; end
    sz = size(funTs(sInd).vec); sz(1) = length(f);
    if memFlag % using this saves memory by averaging power across tapers but looses phase information
        avTapers = 1;
        % anticipate memory needs
        [NT,C] = size(funTs(sInd).vec(:,:,:,1));
        pad = 0;
        NFFT=max(2^(nextpow2(NT)+pad),NT);
        % accordingly choose number of computation blocs
        nBloc = ceil(NFFT*C*K/3e10);
        sz(3) = 1;
    else
        avTapers = 0;
        nBloc = 1;
        sz(3) = size(funPsd(sInd).tp,2);
    end
    funPsd(sInd).vec = nan(sz,class(funTs(sInd).vec));
    
    for runInd = 1:size(funTs(sInd).vec,4)
        if ~isfield(funTs(sInd),'t'); funTs(sInd).t = []; end
        
        %% Compute power spectra
        if ~skipPSD
            if verbose; disp('power spectra: computing'); end
            cohF = [];
            nMode = [];
            [funPsd(sInd).vec(:,:,:,runInd),f] = fastMtPSD(funPsd(sInd).tp,funTs(sInd).vec(:,:,:,runInd),param.Fs,[],nBloc,avTapers);
            funPsd(sInd).f = permute(f,[1 3 4 2]);
            if verbose; disp('power spectra: done'); end
        end

        %% Compute coherence spectra
        if ~skipSVD
            if verbose; disp('coherence spectra: computing'); end
            [funPsd(sInd).svd.u,...
                funPsd(sInd).svd.s,...
                funPsd(sInd).svd.v,...
                funPsd(sInd).svd.coh,...
                funPsd(sInd).svd.f]...
                = fastKleinMtSVD(funPsd(sInd).svd.tp,funTs(sInd).vec(:,:,:,runInd),param.Fs,funTs(sInd).t,cohF,nMode);
            toc
            funPsd(sInd).svd.s = permute(funPsd(sInd).svd.s,[2 1 3]);
            funPsd(sInd).svd.coh = permute(funPsd(sInd).svd.coh,[2 1 3]);
            funPsd(sInd).svd.tp = 'same as for psd';
            funPsd(sInd).svd.info = strjoin({'tapers/vox' 'mode' 'freq/time'},' x ');
            if verbose; disp('coherence spectra: done'); end
        % else
        %     dbstack; error('code that')
        end
    end

    %% Output some stuff
    funPsd(sInd).nframes = size(funPsd(sInd).vec,1)/funPsd(sInd).nruns;
    funPsd(sInd).ntapers = param.tapers(2);
    funPsd(sInd).tr = mode(diff(f))*1000;
    funPsd(sInd).tp = permute(funPsd(sInd).tp,[3 4 5 1 2]);
    if isfield(funPsd(sInd),'psd')
        funPsd(sInd) = rmfield(funPsd(sInd),'psd');
    end
    funPsd(sInd).psd.dim = strjoin({'space' 'freq' 'taper' 'run'},' X ');
    tmp = diff(f,[],4); if any(tmp(:)); dbstack; error('X'); end
    funPsd(sInd).psd.f = f(:,:,:,1);
    funPsd(sInd).psd.w = W;
    funPsd(sInd).psd.tw = TW;
    funPsd(sInd).psd.mask = funTs(sInd).vol2vec;
    if isfield(param,'tp')
        funPsd(sInd).psd.param = rmfield(param,'tp');
    else
        funPsd(sInd).psd.param = param;
    end
end