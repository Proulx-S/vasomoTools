function funPsd = runPSD7(funTs,W,K,win,mask,memFlag,skipSVD,skipPSD,verbose)
% Wrapper for the Chronux's mtspectrumc function for multitaper estimation of
% pds spectra, compatible with MRI data imported by MRIread.m.
%
%    win [int int]: time window width and step size, in seconds
%
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
if ~exist('win','var') || isempty(win); param.win = inf; else, param.win = win; end; clear win
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

if isfield(funTs,'vec') && ~isempty(funTs.vec)
    tmp = all(funTs.vec==0,1);
else
    tmp = all(funTs.vol==0,4);
end
if any(tmp(:))
    warning('Some voxels are all 0s. Adjust your mask to avoid later problems')
end

%% Set parameters
Wflag = exist('W','var') && ~isempty(W);
Kflag = exist('K','var') && ~isempty(K);
tpFlag = ~isempty(tp); if tpFlag; Wflag = false; Kflag = false; end
if ~isfield(funTs,'nruns')
    funTs.nruns = 1;
end
if param.win(1)==inf
    T = tr.*funTs.nframes;
    param.win(1) = T;
    param.win(2) = inf;
else
    param.win = ceil(param.win./tr).*tr;
    if length(param.win)==1; param.win(2) = ceil(param.win(1)/2/tr)*tr; end
    T = param.win(1);
end
if param.win(2)==inf
    if verbose; disp('Param for psd. Svd uses twice K of psd'); end
else
    if verbose; disp('Param for time-reolved psd, full psd uses ceil(K/2). Svd uses twice K of psd'); end
end
if Wflag && Kflag
    error('Cannot specify both W and K');
elseif Wflag
    [TW,W,K] = W2K(T,W);
elseif Kflag || tpFlag
    if tpFlag && verbose
        disp('using precomputed tapers');
    end
    [TW,W,K] = K2W(T,K,verbose);
end
param.tapers = [TW K];


%% xgram parameters
param.win = param.win/tr;
param.win(3) = ceil(funTs.nframes / (param.win(2)));
allWin = repmat(1:param.win(1),[param.win(3) 1]);
allWin = allWin + (((1:param.win(3))-1)*param.win(2))'; % win x t
allWin(any(allWin>funTs.nframes,2),:) = [];
allWin(end+1,:) = (funTs.nframes-param.win(1)+1:funTs.nframes)';
param.win(3) = [];
param.win = param.win*tr;
if verbose && param.win(2)~=inf
    disp(['win(1) (window width): ' num2str(param.win(1)) 'sec'])
    disp(['win(2) (step size)   : ' num2str(param.win(2)) 'sec'])
end




funPsd = funTs;
[funPsd.vol] = deal([]);
[funPsd.vec] = deal([]);
if ~isfield(funTs,'volInfo'); [funTs.volInfo] = deal(strjoin({'X' 'Y' 'Z' 'freq/time' 'taper' 'run'},' x ')); end
if ~isfield(funTs,'vecInfo'); [funTs.vecInfo] = deal(strjoin({'freq/time' 'vox' 'taper' 'run'},' x ')); end
tmp = strsplit(funTs(1).vecInfo,' x '); tmp{1} = 'freq/time'; tmp = strjoin(tmp,' x ');
[funPsd.vecInfo] = deal(tmp);
tmp = strsplit(funTs(1).volInfo,' x '); tmp{4} = 'freq/time'; tmp = strjoin(tmp,' x ');
[funPsd.volInfo] = deal(tmp);
param.Fs = 1/tr;
param.complex = 1;
funTs = vol2vec(funTs);

for sInd = 1:length(funTs)
    %% Get tapers
    if ~isfield(funTs(sInd),'t') || isempty(funTs(sInd).t)
        funTs(sInd).t = (0:funTs(sInd).tr/1000:(funTs(sInd).nframes-1)*funTs(sInd).tr/1000)';
    end
    if ~skipPSD
        paramPsdGram = param;
        if verbose && ~isempty(tp); disp('getting taper for psd'); end
        if param.win(2)==inf
            %%% full spectrum
            paramPsd = param;
            [~,funPsd(sInd).psd.f,funPsd(sInd).psd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramPsd,funTs(sInd).t(:,:,:,:,:,1),tp);
            funPsd(sInd).psd.tp = permute(funPsd(sInd).psd.tp,[3 4 5 1 2]);
        else
            %%% full spectrum
            paramPsd = param;
            % [paramPsd.tapers(1),~,paramPsd.tapers(2)] = K2W(funTs(sInd).nframes.*tr,ceil(param.tapers(2)/2),verbose);
            [~,funPsd(sInd).psd.f,funPsd(sInd).psd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramPsd,funTs(sInd).t(:,:,:,:,:,1),tp);
            funPsd(sInd).psd.tp = permute(funPsd(sInd).psd.tp,[3 4 5 1 2]);
            %%% time-resolved spectrum
            paramPsdGram = param;
            [~,funPsd(sInd).psdGram.f,funPsd(sInd).psdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramPsdGram.win(1)/tr,1,:,1), paramPsdGram,funTs(sInd).t(:,:,:,1:paramPsdGram.win(1)/tr,:,1),tp);
            funPsd(sInd).psdGram.tp = permute(funPsd(sInd).psdGram.tp,[3 4 5 1 2]);
        end
        % % use time windows
        % [~,funPsd(sInd).psdGram.f,funPsd(sInd).psdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramPsdGram.win(1)/tr,1,:,1), paramPsdGram,funTs(sInd).t(:,:,:,1:paramPsdGram.win(1)/tr,:,1),tp);
        % funPsd(sInd).psdGram.tp = permute(funPsd(sInd).psdGram.tp,[3 4 5 1 2]);
        % if paramPsdGram.win(2)~=inf
        %     % prepare for also performng psd over the full timeseries,
        %     % using half as many tapers
        %     paramPsd = paramPsdGram;
        %     [paramPsd.tapers(1),~,paramPsd.tapers(2)] = K2W(funTs(sInd).nframes.*tr,ceil(K/2),verbose);
        %     [~,funPsd(sInd).psd.f,funPsd(sInd).psd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramPsd,funTs(sInd).t(:,:,:,:,:,1),tp);
        %     funPsd(sInd).psd.tp = permute(funPsd(sInd).psd.tp,[3 4 5 1 2]);
        % else
        %     paramPsd = param;
        %     funPsd(sInd).psd.tp = funPsd(sInd).psdGram.tp;
        % end
    else
        dbstack; error('code that')
    end

    if ~skipSVD
        if verbose && ~isempty(tp); disp('getting taper for svd'); end
        if param.win(2)==inf
            %%% full coherence spectrum
            paramSvd = param;
            % [paramSvd.tapers(1),~,paramSvd.tapers(2)] = K2W(funTs(sInd).nframes.*tr,ceil(param.tapers(2)*2),verbose);
            [~,funPsd(sInd).svd.f,funPsd(sInd).svd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramSvd,funTs(sInd).t(1,1,1,:,1,1),tp);
            funPsd(sInd).svd.tp = permute(funPsd(sInd).svd.tp,[2 3 1]);
        else
            %%% full coherence spectrum
            paramSvd = param;
            [~,funPsd(sInd).svd.f,funPsd(sInd).svd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramSvd,funTs(sInd).t(1,1,1,:,1,1),tp);
            funPsd(sInd).svd.tp = permute(funPsd(sInd).svd.tp,[2 3 1]);
            %%% time-resolved coherence spectrum
            paramSvdGram = param;
            % [paramSvdGram.tapers(1),~,paramSvdGram.tapers(2)] = K2W(funTs(sInd).nframes.*tr,ceil(param.tapers(2)*2),verbose);
            [~,funPsd(sInd).svdGram.f,funPsd(sInd).svdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramSvdGram.win(1)/tr,1,:,1),paramSvdGram,funTs(sInd).t(1,1,1,1:paramSvdGram.win(1)/tr,1,1),tp);
            funPsd(sInd).svdGram.tp = permute(funPsd(sInd).svdGram.tp,[2 3 1]);
        end

        % paramSvdGram = paramPsdGram;
        % paramSvdGram.tapers(2) = paramSvdGram.tapers(2)*2;
        % paramSvdGram.tapers(1) = (paramSvdGram.tapers(2)+1)/2;
        % % use time windows
        % [~,funPsd(sInd).svdGram.f,funPsd(sInd).svdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramSvdGram.win(1)/tr,1,:,1),paramSvdGram,funTs(sInd).t(1,1,1,1:paramSvdGram.win(1)/tr,1,1),tp);
        % funPsd(sInd).svdGram.tp = permute(funPsd(sInd).svdGram.tp,[2 3 1]);
        % if paramSvdGram.win(2)~=inf
        %     % prepare for also performng svd over the full timeseries,
        %     % using half as many tapers
        %     paramSvd = paramSvdGram;
        %     [paramSvd.tapers(1),~,paramSvd.tapers(2)] = K2W(funTs(sInd).nframes.*tr,ceil(paramSvdGram.tapers(2)/2),verbose);
        %     [~,funPsd(sInd).svd.f,funPsd(sInd).svd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramSvd,funTs(sInd).t(1,1,1,:,1,1),tp);
        %     funPsd(sInd).svd.tp = permute(funPsd(sInd).svd.tp,[2 3 1]);
        % end
    end

    if ~isempty(tp); tp = []; end
    
    szPsd = size(funTs(sInd).vec);
    szSvd = size(funTs(sInd).vec);
    if param.win(2)~=inf
        szPsdGram = size(funTs(sInd).vec);
        szSvdGram = size(funTs(sInd).vec);
    end
    %%% Number of frequency points
    szPsd(1) = length(funPsd(sInd).psd.f);
    szSvd(1) = length(funPsd(sInd).svd.f);
    if param.win(2)~=inf
        szPsdGram(1) = length(funPsd(sInd).psdGram.f);
        szSvdGram(1) = length(funPsd(sInd).svdGram.f);
    end
    %%% Number of tapers to save
    if memFlag % using this saves memory by averaging power across tapers but looses phase information
        avTapers = 1;
        % anticipate memory needs
        [NT,C] = size(funTs(sInd).vec(:,:,:,1));
        pad = 0;
        NFFT=max(2^(nextpow2(NT)+pad),NT);
        % accordingly choose number of computation blocs
        nBloc = ceil(NFFT*C*K/3e10);
        szPsd(3) = 1;
    else
        avTapers = 0;
        nBloc = 1;
        szPsd(3) = size(funPsd(sInd).psd.tp,2);
    end
    szSvd(3) = 1;
    if param.win(2)~=inf
        szSvdGram(3) = 1;
        szPsdGram(3) = 1;
    end

    %%% Number of runs in the same structure
    szPsd(4) = size(funTs(sInd).vec,4);
    funPsd(sInd).psd.vec = nan(szPsd,class(funTs(sInd).vec));
    szSvd = szPsd; szSvd(1) = length(funPsd(sInd).psd.f);
    
    %%% Number of time windows
    if param.win(2)~=inf
        szPsdGram(5) = size(allWin,1); %number of windows
        funPsd(sInd).psdGram.vec = nan(szPsdGram([]),class(funTs(sInd).vec)); % freq taper vox run
        szSvdGram(5) = size(allWin,1); %number of windows
    end


    sz = [length(funPsd(sInd).psd.f) paramPsd.tapers(2) size(funTs(sInd).vec,2) size(funTs(sInd).vec,4)];
    if avTapers; sz(2) = 1; end
    funPsd(sInd).psd.vec = nan(sz,class(funTs(sInd).vec)); % freq taper vox run

    if param.win(2)~=inf
        sz = [length(funPsd(sInd).psdGram.f) paramPsdGram.tapers(2) size(funTs(sInd).vec,2) size(funTs(sInd).vec,4)];
        sz(2) = 1; sz(5) = size(allWin,1);
        funPsd(sInd).psdGram.vec = nan(sz,class(funTs(sInd).vec)); % freq taper vox run
    end

    for runInd = 1:szPsd(4)
        if verbose && szPsd(4)>1; disp(['---Run ' num2str(runInd) '/' num2str(szPsd(4)) '---']); end
        %% Compute power spectra
        if ~skipPSD
            %%% over the full time series
            if verbose; disp('full power spectra: computing'); end
            [funPsd(sInd).psd.vec(:,:,:,runInd),funPsd(sInd).psd.f] = fastMtPSD(funPsd(sInd).psd.tp,funTs(sInd).vec(:,:,:,runInd),paramPsd.Fs,[],nBloc,avTapers);
            funPsd(sInd).psd.T(:,:,:,runInd) = size(funTs(sInd).vec(:,:,:,runInd),1).*tr;
            if verbose; disp('full power spectra: done'); end
            %%% over each time window
            if param.win(2)~=inf
                for winInd = 1:size(allWin,1)
                    if verbose; disp(['time-resolved power spectra: computing window ' num2str(winInd) '/' num2str(size(allWin,1))]); end
                    if runInd==1 && winInd==1; curVerbose = 2; else curVerbose = 1; end
                    [funPsd(sInd).psdGram.vec(:,:,:,runInd,winInd),funPsd(sInd).psdGram.f] = fastMtPSD(funPsd(sInd).psdGram.tp,funTs(sInd).vec(allWin(winInd,:),:,:,runInd)-mean(funTs(sInd).vec(allWin(winInd,:),:,:,runInd),1),paramPsdGram.Fs,[],nBloc,1,curVerbose);
                    funPsd(sInd).psdGram.T(:,:,:,runInd,winInd) = size(funTs(sInd).vec(allWin(winInd,:),:,:,runInd),1).*tr;
                    funPsd(sInd).psdGram.tWin(1,1,1,runInd,winInd,1) = mean(funTs(sInd).t(1,1,1,allWin(winInd,[1 end]),1,runInd));
                end
                if verbose; disp('time-resolved power spectra: done'); end
            end
        end

        %% Compute coherence spectra
        if ~skipSVD
            cohF = [0 inf];
            nMode = 3;

            %%% Over the full time series
            if verbose; disp('full coherence spectra: computing'); end
            [funPsd(sInd).svd.u(:,:,:,runInd),...
                funPsd(sInd).svd.s(:,:,:,runInd),...
                funPsd(sInd).svd.v(:,:,:,runInd),...
                funPsd(sInd).svd.coh(:,:,:,runInd),...
                funPsd(sInd).svd.f]...
                = fastKleinMtSVD(funPsd(sInd).svd.tp,funTs(sInd).vec(:,:,:,runInd),paramSvd.Fs,funTs(sInd).t(:,:,:,:,:,runInd),cohF,nMode);
            funPsd(sInd).svd.T(:,:,:,runInd) = size(funTs(sInd).vec(:,:,:,runInd),1).*tr;
            if verbose; disp('full coherence spectra: done'); end

            %%% Over each time window
            if param.win(2)~=inf
                for winInd = 1:size(allWin,1)
                    if verbose; disp(['time-resolved coherence spectra: computing window ' num2str(winInd) '/' num2str(size(allWin,1))]); end
                    [funPsd(sInd).svdGram.u(:,:,:,runInd,winInd),...
                        funPsd(sInd).svdGram.s(:,:,:,runInd,winInd),...
                        funPsd(sInd).svdGram.v(:,:,:,runInd,winInd),...
                        funPsd(sInd).svdGram.coh(:,:,:,runInd,winInd),...
                        funPsd(sInd).svdGram.f]...
                        = fastKleinMtSVD(funPsd(sInd).svdGram.tp,funTs(sInd).vec(allWin(winInd,:),:,:,runInd)-mean(funTs(sInd).vec(allWin(winInd,:),:,:,runInd),1),paramSvdGram.Fs,funTs(sInd).t(:,:,:,allWin(winInd,:),:,runInd),cohF,nMode);
                    funPsd(sInd).svdGram.T(:,:,:,runInd,winInd) = size(funTs(sInd).vec(allWin(winInd,:),:,:,runInd),1).*tr;
                    funPsd(sInd).svdGram.tWin(1,1,1,runInd,winInd,1) = mean(funTs(sInd).t(1,1,1,allWin(winInd,[1 end]),1,runInd));
                end
                if verbose; disp('time-resolved coherence spectra: done'); end
            end
        end
    end

    %% Sort outputs
    funPsd(sInd).psd.f = permute(funPsd(sInd).psd.f,[2 1 3 4 5 6]);
    funPsd(sInd).psd.tp = permute(funPsd(sInd).psd.tp,[4 5 1 2 3 6]);
    funPsd(sInd).psd.vec;
    funPsd(sInd).psd.K = paramPsd.tapers(2);
    [~,funPsd(sInd).psd.W,~] = K2W(funPsd(sInd).psd.T,funPsd(sInd).psd.K,0);
    funPsd(sInd).psd.param = paramPsd;
    funPsd(sInd).psd.info = strjoin({'freq/time' 'tapers' 'vox' 'run' 'timeWindow' 'mode'},' x ');
    funPsd(sInd).psd.mask = funTs(sInd).vol2vec;

    if param.win(2)~=inf
        funPsd(sInd).psdGram.f = permute(funPsd(sInd).psdGram.f,[2 1 3 4 5 6]);
        funPsd(sInd).psdGram.tp = permute(funPsd(sInd).psdGram.tp,[4 5 1 2 3 6]);
        funPsd(sInd).psdGram.vec;
        % funPsd(sInd).psdGram.t;
        % funPsd(sInd).psdGram.tWin = mean(funPsd(sInd).psdGram.t,1);
        funPsd(sInd).psdGram.lWin = paramPsdGram.win(1);
        funPsd(sInd).psdGram.K = paramPsdGram.tapers(2);
        [~,funPsd(sInd).psdGram.W,~] = K2W(funPsd(sInd).psdGram.T,funPsd(sInd).psdGram.K,0);
        funPsd(sInd).psdGram.param = paramPsdGram;
        funPsd(sInd).psdGram.info = strjoin({'freq/time' 'tapers' 'vox' 'run' 'timeWindow' 'mode'},' x ');
        funPsd(sInd).psdGram.mask = funTs(sInd).vol2vec;
    end
    
    funPsd(sInd).svd.f;
    funPsd(sInd).svd.tp = permute(funPsd(sInd).svd.tp,[2 4 3 5 6 1]);
    funPsd(sInd).svd.u;
    funPsd(sInd).svd.s = permute(funPsd(sInd).svd.s,[2 1 3 4 5 6]);
    funPsd(sInd).svd.v = permute(funPsd(sInd).svd.v,[6 2 3 4 5 1]);
    funPsd(sInd).svd.coh = permute(funPsd(sInd).svd.coh,[2 1 3 4 5 6]);
    funPsd(sInd).svd.K = paramSvd.tapers(2);
    [~,funPsd(sInd).svd.W,~] = K2W(funPsd(sInd).svd.T,funPsd(sInd).svd.K,0);
    funPsd(sInd).svd.param = paramSvd;
    funPsd(sInd).svd.info = strjoin({'vox' 'mode' 'freq/time' 'run' 'timeWindow' 'tapers'},' x ');
    funPsd(sInd).svd.mask = funTs(sInd).vol2vec;

    if param.win(2)~=inf
        funPsd(sInd).svdGram.f;
        funPsd(sInd).svdGram.tp = permute(funPsd(sInd).svdGram.tp,[2 4 3 5 6 1]);
        funPsd(sInd).svdGram.u;
        funPsd(sInd).svdGram.s = permute(funPsd(sInd).svdGram.s,[2 1 3 4 5 6]);
        funPsd(sInd).svdGram.v = permute(funPsd(sInd).svdGram.v,[6 2 3 4 5 1]);
        funPsd(sInd).svdGram.coh = permute(funPsd(sInd).svdGram.coh,[2 1 3 4 5 6]);
        % funPsd(sInd).svdGram.t;
        % funPsd(sInd).svdGram.tWin = mean(funPsd(sInd).svdGram.t,1);
        funPsd(sInd).svdGram.lWin = paramSvdGram.win(1);
        funPsd(sInd).svdGram.K = paramSvdGram.tapers(2);
        [~,funPsd(sInd).svdGram.W,~] = K2W(funPsd(sInd).svdGram.T,funPsd(sInd).svdGram.K,0);
        funPsd(sInd).svdGram.param = paramSvdGram;
        funPsd(sInd).svdGram.info = strjoin({'vox' 'mode' 'freq/time' 'run' 'timeWindow' 'tapers'},' x ');
        funPsd(sInd).svdGram.mask = funTs(sInd).vol2vec;
    end
    
    %% For backward compatibility
    funPsd(sInd).vec = permute(funPsd(sInd).psd.vec,[1 3 2 4]); funPsd(sInd).psd.vec = [];
    funPsd(sInd).f = permute(funPsd(sInd).psd.f,[2 3 4 5 6 1]);
    funPsd(sInd).K = funPsd(sInd).psd.K;
    funPsd(sInd).T = funPsd(sInd).psd.T;
    funPsd(sInd).W = funPsd(sInd).psd.W;
    funPsd(sInd).tr = mode(diff(funPsd(sInd).f))*1000;
    funPsd(sInd).nfreq = size(funPsd(sInd).vec,1);
    if funPsd(sInd).nframes~=size(funTs(sInd).vec,1); dbstack; warning('something wrong with nframes'); end
    if funPsd(sInd).nruns~=size(funPsd(sInd).vec,4); dbstack; warning('something wrong with nframes'); end
    funPsd(sInd).ntapers = size(funPsd(sInd).vec,3);
end
    