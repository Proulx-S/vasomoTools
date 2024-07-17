function [funPsd,tp] = runPSD5(funTs,W,K,mask,memFlag)
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
if ~exist('memFlag','var'); memFlag = false; end
if ~exist('W','var'); W = []; end
if ~exist('K','var'); K = []; end
if isempty(K) && isempty(W); K = 1; end
if length(K)>1
    tp = K; K = size(tp,2);
else
    tp = [];
end
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

    display(['W  (halfwidth)       : ' num2str(Wreq,'%0.5f ')])
    display(['W  (halfwidth)       : ' num2str(W,'%0.5f ')])
    display(['TW (time-halfwidth)  : ' num2str(TW)])
    display(['K  (number of tapers): ' num2str(K)])
elseif Kflag
    TW = (K+1)/2;
    T = tr.*funTs.nframes*funTs.nruns;
    W = TW/T;
    param.tapers = [TW K];
    display(['K  (number of tapers): ' num2str(K)])
    display(['W  (halfwidth)       : ' num2str(W,'%0.5f ')])
    display(['TW (time-halfwidth)  : ' num2str(TW)])
elseif tpFlag
    TW = (K+1)/2;
    T = tr.*funTs.nframes*funTs.nruns;
    W = TW/T;
    param.tapers = [TW K];
    display(['using precomputed tapers']);
    display(['K  (number of tapers): ' num2str(K)])
    display(['W  (halfwidth)       : ' num2str(W,'%0.5f ')])
    display(['TW (time-halfwidth)  : ' num2str(TW)])
end

funPsd = funTs;
[funPsd.vol] = deal([]);
[funPsd.vec] = deal([]);
tmp = strsplit(funTs(1).vecInfo,' x '); tmp{1} = 'freq'; tmp = strjoin(tmp,' x ');
[funPsd.vecInfo] = deal(tmp);
tmp = strsplit(funTs(1).volInfo,' x '); tmp{1} = 'freq'; tmp = strjoin(tmp,' x ');
[funPsd.volInfo] = deal(tmp);
param.Fs = 1/tr;
param.complex = 1;
for sInd = 1:length(funTs)
    if ~isfield(funTs(sInd),'t')
        funTs(sInd).t = [];
    end
    [~,f,tp] = mtspectrumc4(funTs(sInd).vec(:,1), param,funTs(sInd).t,tp);
    sz = size(funTs(sInd).vec); sz(1) = length(f);

    if memFlag % using this saves memory by averaging power across tapers but looses phase information
        % anticipate memory needs
        [NT,C] = size(funTs(sInd).vec(:,:,:,1));
        pad = 0;
        NFFT=max(2^(nextpow2(NT)+pad),NT);
        % accordingly choose number of computation blocs
        nBloc = ceil(NFFT*C*K/3e10);
    else
        sz(3) = size(tp,2);
    end
    funPsd(sInd).vec = nan(sz,class(funTs(sInd).vec));
    for runInd = 1:size(funTs(sInd).vec,4)
        if ~isfield(funTs(sInd),'t'); funTs(sInd).t = []; end
        [funPsd(sInd).vec(:,:,:,runInd),~] = mtspectrumc4(funTs(sInd).vec(:,:,:,runInd),param,[],tp,memFlag);
    end
    funPsd(sInd).nframes = size(funPsd(sInd).vec,1)/funPsd(sInd).nruns;
    % if size(funTs(sInd).vec,4)>1
    %     %% ComputeRun by run
    %     [~,~,f] = mtspectrumc2(funTs(sInd).vec(:,1), param);
    %     sz = size(funTs(sInd).vec); sz(1) = length(f); sz(3) = param.tapers(2);
    %     vec = nan(sz);
    %     parfor i = 1:size(funTs(sInd).vec,4)
    %         [~,vec(:,:,:,i),f] = mtspectrumc3(funTs(sInd).vec(:,:,:,i), param);
    %     end
    %     funPsd(sInd).vec = vec;
    %     funPsd(sInd).nframes = size(funPsd(sInd).vec,1);
    % else
    %     if ~isfield(funTs(sInd),'t')
    %         funTs(sInd).t = [];
    %     end
    %     % Compute spectra, using missing data Slepian tapers only if necessary
    %     [~,funPsd(sInd).vec,f] = mtspectrumc3(funTs(sInd).vec, param,funTs(sInd).t);
    % 
    %     if ~isfield(funPsd(sInd),'nruns') || isempty(funPsd(sInd).nruns); funPsd(sInd).nruns = 1; end
    %     funPsd(sInd).nframes = size(funPsd(sInd).vec,1)/funPsd(sInd).nruns;
    %     % figure('WindowStyle','docked');
    %     % plot(f,mean(mean(conj(funPsd.vec).*funPsd.vec,3),2))
    %     % [~,funPsd.vec,f] = mtspectrumc2(funTs.vec, param);
    % end
    funPsd(sInd).ntapers = param.tapers(2);
    funPsd(sInd).tr = mode(diff(f))*1000;

    %% Output some stuff
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