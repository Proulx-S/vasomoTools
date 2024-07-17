function funPsd = runFullMT2(funTs,W,K,win,onsets,ondurs,mask,memFlag,skipSVD,skipPSD,verbose,taperPerm,phaseRand)
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
%    Alternatively, precomputed tapers can be input as K to save
%    time (useful when missing data slepian are used).
%    See funPsd.psd for other useful parameters.
%    funPsd.tr reflects the frequency resolution in Hz*1000
if ~exist('verbose','var') || isempty(verbose); verbose = true; end
if ~exist('memFlag','var') || isempty(memFlag); memFlag = false; end
if ~exist('W','var'); W = []; end
if ~exist('K','var'); K = []; end
if isempty(K) && isempty(W); K = 1; end
if ~exist('skipSVD','var') || isempty(skipSVD); skipSVD = false; end
if ~exist('skipPSD','var') || isempty(skipPSD); skipPSD = false; end
if ~exist('mask','var'); mask = []; end
if ~exist('onsets','var'); onsets = []; end
if ~exist('ondurs','var'); ondurs = []; end

if isempty(onsets); onsets = funTs.dsgn.onsets; end
if isempty(ondurs); ondurs = funTs.dsgn.ondurs; end

% if K==1; skipSVD = true; end

if iscell(funTs)
    for I = 1:numel(funTs)
        funPsd{I} = runFullMT2(funTs{I},W,K,win,onsets,ondurs,mask,memFlag,skipSVD,skipPSD,verbose,taperPerm,phaseRand);
    end
elseif isstruct(funTs)
    for I = 1:numel(funTs)
        cohFperm = [0 1.1];
        funPsd(I) = doIt(funTs(I),W,K,win,onsets,ondurs,mask,memFlag,skipSVD,skipPSD,verbose,taperPerm,phaseRand,[],cohFperm);

        % % taperPerm = 2^7;
        % % phaseRand = 0;
        % for K = [4 8]
        %     cohFperm = [0 1.1];
        %     funPsd(I) = doIt(funTs(I),W,K,win,mask,memFlag,skipSVD,skipPSD,verbose,taperPerm,phaseRand,[],cohFperm);
        %     funPsd(I).svd.nPerm = taperPerm;
        %     % funTs(2:end) = [];
        %     % save tmp2 -v7.3
        %     close all
        %     % f0List = [0 0.00301408 0.0572676 0.072338 0.180845 0.186873 0.277296 0.72338 0.747492 0.940394];
        %     % f0List = [1.33825 1.46183 1.46484];
        %     % f0List = [1.33222 1.36538 1.46183];
        %     % f0List = [0.093465 0.271267 0.60583];
        %     f0List = [0.0542535 0.60583];
        %     plotPerm4(funPsd(I).svd,f0List)
        %     % plotPerm4(funPsd(I).svd)
        % end
        % keyboard
        % % funTs(2:end) = [];
        % % save tmp -v7.3
    end
else
    dbstack; error('this should not happen')
end





function funPsd = doIt(funTs,W,K,win,onsetList,durList,mask,memFlag,skipSVD,skipPSD,verbose,taperPerm,phaseRand,cohF,cohFperm)
if length(K)>1; tp = K; K = size(tp,2); else tp = []; end
if ~exist('win','var');             win = []; end
if ~exist('onsetList','var'); onsetList = []; end
if ~exist('durList','var');     durList = []; end

if isempty(win); win = [inf 0]; end



if funTs.nvoxels==1; if verbose; disp('only one timeseries, skipping SVD'); end; skipSVD = true; end

param.win = win;
param.onsetList = onsetList;
param.durList = durList;

if param.win(1)==inf; skipGram = true; skipTrialGram = true; else skipGram = false; skipTrialGram = false; end
if isempty(param.onsetList); skipTrialGram = true; end


% if length(param.win)>2
%     param.onsetList = param.win(3:end)';
%     param.win(3:end) = [];
% else
%     param.onsetList = [];
% end


%% Mask
if ~isempty(mask)
    funTs = applyMask(funTs,mask);
end

%% Assert
if isfield(funTs,'vec') && ~isempty(funTs.vec); tmp = all(funTs.vec==0,1); else tmp = all(funTs.vol==0,4); end
if any(tmp(:)); warning('Some voxels are all 0s. Adjust your mask to avoid later problems'); end

%% Set parameters
tr = funTs.tr/1000;
Wflag = ~isempty(W);
Kflag = ~isempty(K);
tpFlag = ~isempty(tp); if tpFlag; Wflag = false; Kflag = false; end
if ~isfield(funTs,'nruns'); funTs.nruns = 1; end
%%% Window size
if param.win(1)==inf
    %%%% single-window over the full time series
    T = tr.*funTs.nframes;
    param.win(1) = funTs.nframes;
    param.win(2) = 0;
else
    %%%% multiple time windows
    % param.win = ceil(param.win./tr).*tr;
    if length(param.win)==1
        param.win(2) = 1;
    end
    T = param.win(1)*tr;
end
% if param.win(2)==inf
%     if verbose; disp('Param for psd. Svd uses twice K of psd'); end
% else
%     if verbose; disp('Param for time-reolved psd, full psd uses ceil(K/2). Svd uses twice K of psd'); end
% end
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
if ~skipGram
    param.win(3) = ceil(funTs.nframes / (param.win(2)));
    allWin = repmat(1:param.win(1),[param.win(3) 1]);
    allWin = allWin + (((1:param.win(3))-1)*param.win(2))'; % win x t
    allWin(any(allWin>funTs.nframes,2),:) = [];
    allWin(end+1,:) = (funTs.nframes-param.win(1)+1:funTs.nframes)';
    param.win(3) = [];
    param.win = param.win*tr;
    if verbose && param.win(2)~=inf
        disp(['win(1) (window width): ' num2str(param.win(1),'%0.3f') 'sec or ' num2str(param.win(1)/tr) 'vol'])
        disp(['win(2) (step size)   : ' num2str(param.win(2),'%0.3f') 'sec or ' num2str(param.win(2)/tr) 'vol'])
    end
else
    allWin = [];
end

allWin = unique(allWin,'rows');


%% trial-locked xgram parameters
if ~skipTrialGram
    allWin; % [win X timeIndex]
    onsetList = param.onsetList;
    winSz = param.win(1)./tr;
    n = max(allWin(:));
    nWin = size(allWin,1);
    nTrial = size(onsetList,1);
    allWin2 = repmat({zeros(n,nWin)},[nTrial 1]);
    for winInd = 1:nWin
        for trialInd = 1:nTrial
            if trialInd == 1
                allWin2{trialInd}(allWin(winInd,:),winInd) = 1;
            else
                offsetInd = floor((onsetList(trialInd) - onsetList(1)) ./ tr);
                tInd = allWin(winInd,:) + offsetInd;
                tInd(tInd>n) = [];
                allWin2{trialInd}(tInd,winInd) = 1;
            end
        end
    end
    allWin3 = any(cat(3,allWin2{:}),3); % [win X time]

    %%% remove windows exceeding timeseries
    endInd = find(allWin3(end,:)==1,1)+1;
    for trialInd = 1:nTrial
        allWin2{trialInd}(:,endInd:end) = [];
    end
    % allWin3(:,endInd:end) = [];

    %%% remove completely overlapping windows
    endInd = find(sum((allWin2{1}+allWin2{2}(:,1))==2,1)==winSz);
    for trialInd = 1:nTrial
        allWin2{trialInd}(:,endInd:end) = [];
    end
    % allWin3(:,endInd:end) = [];

    % %%% visualize
    % figure('WindowStyle','docked');
    % ht = tiledlayout(1,nTrial+1); ax = {};
    % for trialInd = 1:nTrial
    %     ax{trialInd} = nexttile([1 1]);
    %     imagesc(allWin2{trialInd})
    % end
    % ax{trialInd} = nexttile([1 1]);
    % imagesc(allWin3)

    %%% put back to standard index format
    % winSz2 = unique(sum(allWin3,1));
    % if length(winSz2)~=1; error('badly defined trial-locked windows'); end
    % allWin4 = zeros(size(allWin3,2),winSz2);
    % for winInd = 1:size(allWin3,2)
    %     allWin4(winInd,:) = find(allWin3(:,winInd));
    % end
    % allWinTrialLock = allWin4; % [win X timeIndex]
    allWin3 = cat(3,allWin2{:}); % [win X time X trial]
    winSz2 = unique(sum(allWin3,1));
    if length(winSz2)~=1; dbstack; error('badly defined trial-locked windows'); end
    allWin4 = zeros(winSz2,size(allWin3,2),size(allWin3,3));
    for winInd = 1:prod(size(allWin3,[2 3]))
        allWin4(:,winInd) = find(allWin3(:,winInd));
    end
    allWinTrialLock = permute(allWin4,[2 1 3]); % [win X timeIndex]
    % reshape(allWinTrialLock,size(allWinTrialLock,1),prod(size(allWinTrialLock,[2 3])))
    clear allWin2 allWin3 allWin4 winSz2
else
    allWinTrialLock = [];
end


%% initiate stuff
funPsd = funTs;
[funPsd.vol] = deal([]);
[funPsd.vec] = deal([]);
if ~isfield(funTs,'volInfo'); [funTs.volInfo] = deal(strjoin({'X' 'Y' 'Z' 'time/freq' 'taper/mode' 'run'},' x ')); end
if ~isfield(funTs,'vecInfo'); [funTs.vecInfo] = deal(strjoin({'time/freq' 'vox' 'taper/mode' 'run'},' x ')); end
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
    if skipPSD
        dbstack; error('code that')
    else

        %%% full timeseries
        if verbose; disp('getting taper for full timeseries'); end
        K   = param.tapers(2);
        tr  = funTs(sInd).tr/1000;
        N   = funTs(sInd).nframes;
        [TP.full.tp,TP.full.eigs] = getTapers(K,tr,N);
        TP.full.t = funTs(sInd).t;
        % [~,funPsd(sInd).psd.f,funPsd(sInd).psd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramPsd,funTs(sInd).t(:,:,:,:,:,1),tp);

        %%% time-resolved
        if ~skipGram
            if verbose; disp('getting taper for time-resolved analysis'); end
            K   = param.tapers(2);
            tr  = funTs(sInd).tr/1000;
            N   = size(allWin,2);
            [TP.gram.tp,TP.gram.eigs] = getTapers(K,tr,N);
            TP.gram.t = funTs.t(allWin(1,:));
        else
            TP.gram.tp = [];
            TP.gram.eigs = [];
            TP.gram.t = [];
        end

        %%% trial-locked
        if ~skipTrialGram
            %%%% Regular tapers repeated at each trial
            E = length(param.onsetList);
            TP.trialGram.tp = repmat(TP.gram.tp,[1 1 E]);
            TP.trialGram.eigs = TP.gram.eigs;
            TP.trialGram.t = TP.gram.t;
            
            %%%% Missing data tapers (experimental)
            K   = param.tapers(2);
            tr  = funTs(sInd).tr/1000;
            t   = funTs(sInd).t;
            t   = t(allWinTrialLock(1,:));
            N   = length(t);
            [TP.trialGramMD.tp,TP.trialGramMD.eigs] = getTapers(K,tr,N,t);
            TP.trialGramMD.t = t;
            
            
            % if verbose; disp('getting taper for trial-locked analysis'); end
            % %%% Defining KE cross-trial tapers as K regular taper per E
            % %%% events (trials)
            % K   = param.tapers(2);
            % tr  = funTs(sInd).tr/1000;
            % N   = size(allWin,2);
            % TP.trialGram = getTapers(K,tr,N);
            % E = length(param.onsetList);
            % % eval(['TP.trialGram = blkdiag(' strjoin(repmat({'TP.trialGram'},1,E),',') ');']);
            % eval(['TP.trialGram = reshape(blkdiag(' strjoin(repmat({'TP.trialGram'},1,E),',') '),[N*E K E]);']);
            % % TPx = nan(size(funTs(sInd).t,1),size(TP.trialGram,2));
            % % TPx(allWinTrialLock(1,:),:) = TP.trialGram;
            % % plot(funTs(sInd).t,TPx)
            %
            %
            % %%% Defining missing data tapers
            % % if verbose; disp('getting taper for trial-locked analysis'); end
            % % K   = param.tapers(2);
            % % tr  = funTs(sInd).tr/1000;
            % % N   = false(1,funTs(sInd).nframes); N(allWinTrialLock(1,:)) = true; % here N is a logical vector specifying included data, producing tapers appropriate for timeseries with missing data
            % % TP.trialGram = getTapers(K,tr,N);
        else
            TP.trialGram.tp     = [];
            TP.trialGram.eigs   = [];
            TP.trialGram.t      = [];
            TP.trialGramMD.tp   = [];
            TP.trialGramMD.eigs = [];
            TP.trialGramMD.t    = [];
        end
        TP.full.info = 'time x taper x trial';
        TP.gram.info = 'time x taper x trial';
        TP.trialGram.info = 'time x taper x trial';
        TP.trialGramMD.info = 'time x taper x trial';
    end

    %% Compute (also with taper-level permutation
    skip.psd       = skipPSD;
    skip.svd       = skipSVD;
    skip.gram      = skipGram;
    skip.trialGram = skipTrialGram;
    paramInit = param;
    param.psd          = rmfield(paramInit,{'onsetList' 'durList'});
    param.svd          = rmfield(paramInit,{'onsetList' 'durList'});
    param.psdGram      = rmfield(paramInit,{'onsetList' 'durList'});
    param.svdGram      = rmfield(paramInit,{'onsetList' 'durList'});
    param.psdTrialGram = paramInit;
    param.svdTrialGram = paramInit;
    clear paramInit
    param.psdGram.winInd      = allWin;
    param.svdGram.winInd      = allWin;
    param.psdTrialGram.winInd = allWinTrialLock;
    param.svdTrialGram.winInd = allWinTrialLock;

    param.perm = taperPerm;


    % tic
    [funPsd.psd,funPsd.psdGram,funPsd.psdTrialGram,funPsd.psdTrialGramMD,funPsd.svd,funPsd.svdGram,funPsd.svdTrialGram,funPsd.svdTrialGramMD]...
        = computeAll(funTs,TP,param,skip,verbose);
    % disp('+++++')
    % toc
    % disp('+++++')
    % %% Compute with fourrier domain phase randomization
    % tic
    % if phaseRand
    %     dbstack; error('double check that')
    %     % cohF = [0 0.2];
    %     nf = nnz(funPsd.svd.f>=cohFperm(1) & funPsd.svd.f<=cohFperm(2));
    %     funTsShuf = funTs;
    %     % sz = size(funPsd.svd.u,1:5); sz(3) = nf;
    %     % funPsd.svd.uPhaseRand = zeros([size(funPsd.svd.u,1:5) phaseRand],'single');
    %     sz = size(funPsd.svd.s,1:5); sz(3) = nf;
    %     funPsd.svd.sPhaseRand = zeros([sz phaseRand],'single');
    %     % sz = size(funPsd.svd.v,1:5); sz(3) = nf;
    %     % funPsd.svd.vPhaseRand = zeros([sz phaseRand],'single');
    %     sz = size(funPsd.svd.coh,1:5); sz(3) = nf;
    %     funPsd.svd.cohPhaseRand = zeros([sz phaseRand],'single');
    %     % sz = size(funPsd.svdGram.u,1:5); sz(3) = nf;
    %     % funPsd.svdGram.uPhaseRand = zeros([sz phaseRand],'single');
    %     % sz = size(funPsd.svdGram.s,1:5); sz(3) = nf;
    %     % funPsd.svdGram.sPhaseRand = zeros([sz phaseRand],'single');
    %     % sz = size(funPsd.svdGram.v,1:5); sz(3) = nf;
    %     % funPsd.svdGram.vPhaseRand = zeros([sz phaseRand],'single');
    %     % sz = size(funPsd.svdGram.coh,1:5); sz(3) = nf;
    %     % funPsd.svdGram.cohPhaseRand = zeros([sz phaseRand],'single');
    %     for phaseRandInd = 1:phaseRand
    %         disp(['phase randomization ' num2str(phaseRandInd) '/' num2str(phaseRand)])
    %         [X,Y] = pol2cart(rand(size(funTs(sInd).vec)).*(2*pi),abs(fft(funTs(sInd).vec,[],1)));
    %         funTsShuf(sInd).vec = ifft(complex(X,Y));
    %         paramTmp = param; paramTmp.win = [inf inf];
    %         [~,~,svd,svdGram] = computeAll(funTsShuf,TP,szPsd,1,skipSVD,nBloc,avTapers,paramPsd,paramSvd,0,tr,paramTmp,allWin,paramPsdGram,paramSvdGram,[],cohFperm,[]);
    %         % funPsdShuf = computeAll(funTsShuf,funPsd,szPsd,1,skipSVD,nBloc,avTapers,paramPsd,paramSvd,0,sInd,tr,param,allWin,paramPsdGram,paramSvdGram);
    %         % funPsd.svd.uPhaseRand(:,:,:,:,:,phaseRandInd) = svd.u;
    %         funPsd.svd.sPhaseRand(:,:,:,:,:,phaseRandInd) = svd.s;
    %         % funPsd.svd.vPhaseRand(:,:,:,:,:,phaseRandInd) = svd.v;
    %         funPsd.svd.cohPhaseRand(:,:,:,:,:,phaseRandInd) = svd.coh;
    %         % funPsd.svdGram.uPhaseRand(:,:,:,:,:,phaseRandInd) = svdGram.u;
    %         % funPsd.svdGram.sPhaseRand(:,:,:,:,:,phaseRandInd) = svdGram.s;
    %         % funPsd.svdGram.vPhaseRand(:,:,:,:,:,phaseRandInd) = svdGram.v;
    %         % funPsd.svdGram.cohPhaseRand(:,:,:,:,:,phaseRandInd) = svdGram.coh;
    %     end
    % end
    % disp('+++++')
    % toc
    % disp('+++++')



    %% Sort outputs
    funPsd(sInd).param = param;


    % funPsd(sInd).psd.f = permute(funPsd(sInd).psd.f,[2 1 3 4 5 6]);
    % funPsd(sInd).psd.tp = TP.full;
    % funPsd(sInd).psd.vec;
    % funPsd(sInd).psd.K = param.psd.tapers(2);
    % [~,funPsd(sInd).psd.W,~] = K2W(funPsd(sInd).psd.T,funPsd(sInd).psd.K,0);
    % funPsd(sInd).psd.param = param.psd;
    % funPsd(sInd).psd.info = strjoin({'freq/time' 'tapers' 'vox' 'run' 'timeWindow' 'mode'},' x ');
    % funPsd(sInd).psd.mask = funTs(sInd).vol2vec;
    %
    % if ~skipGram
    %     funPsd(sInd).psdGram.f = permute(funPsd(sInd).psdGram.f,[2 1 3 4 5 6]);
    %     funPsd(sInd).psdGram.tp = TP.gram;
    %     funPsd(sInd).psdGram.vec;
    %     % funPsd(sInd).psdGram.t;
    %     % funPsd(sInd).psdGram.tWin = mean(funPsd(sInd).psdGram.t,1);
    %     funPsd(sInd).psdGram.lWin = param.psdGram.win(1);
    %     funPsd(sInd).psdGram.K = param.psdGram.tapers(2);
    %     [~,funPsd(sInd).psdGram.W,~] = K2W(funPsd(sInd).psdGram.T,funPsd(sInd).psdGram.K,0);
    %     funPsd(sInd).psdGram.param = param.psdGram;
    %     funPsd(sInd).psdGram.info = strjoin({'freq/time' 'tapers' 'vox' 'run' 'timeWindow' 'mode'},' x ');
    %     funPsd(sInd).psdGram.mask = funTs(sInd).vol2vec;
    % end
    %
    % if ~skipSVD
    %     funPsd(sInd).svd.f;
    %     funPsd(sInd).svd.tp = permute(TP.full,[3 4 1 5 6 2]);
    %     funPsd(sInd).svd.u;
    %     funPsd(sInd).svd.s = permute(funPsd(sInd).svd.s,[2 1 3 4 5 6]);
    %     funPsd(sInd).svd.v = permute(funPsd(sInd).svd.v,[6 2 3 4 5 1]);
    %     funPsd(sInd).svd.coh = permute(funPsd(sInd).svd.coh,[2 1 3 4 5 6]);
    %     funPsd(sInd).svd.K = param.svd.tapers(2);
    %     [~,funPsd(sInd).svd.W,~] = K2W(funPsd(sInd).svd.T,funPsd(sInd).svd.K,0);
    %     funPsd(sInd).svd.param = param.svd;
    %     funPsd(sInd).svd.info = strjoin({'vox' 'mode' 'freq/time' 'run' 'timeWindow' 'tapers'},' x ');
    %     funPsd(sInd).svd.mask = funTs(sInd).vol2vec;
    % else
    %     funPsd(sInd).svd = [];
    % end
    %
    % if ~skipGram && ~skipSVD
    %     funPsd(sInd).svdGram.f;
    %     funPsd(sInd).svdGram.tp = permute(TP.gram,[3 4 1 5 6 2]);
    %     funPsd(sInd).svdGram.u;
    %     funPsd(sInd).svdGram.s = permute(funPsd(sInd).svdGram.s,[2 1 3 4 5 6]);
    %     funPsd(sInd).svdGram.v = permute(funPsd(sInd).svdGram.v,[6 2 3 4 5 1]);
    %     funPsd(sInd).svdGram.coh = permute(funPsd(sInd).svdGram.coh,[2 1 3 4 5 6]);
    %     % funPsd(sInd).svdGram.t;
    %     % funPsd(sInd).svdGram.tWin = mean(funPsd(sInd).svdGram.t,1);
    %     funPsd(sInd).svdGram.lWin = param.svdGram.win(1);
    %     funPsd(sInd).svdGram.K = param.svdGram.tapers(2);
    %     [~,funPsd(sInd).svdGram.W,~] = K2W(funPsd(sInd).svdGram.T,funPsd(sInd).svdGram.K,0);
    %     funPsd(sInd).svdGram.param = param.svdGram;
    %     funPsd(sInd).svdGram.info = strjoin({'vox' 'mode' 'freq/time' 'run' 'timeWindow' 'tapers'},' x ');
    %     funPsd(sInd).svdGram.mask = funTs(sInd).vol2vec;
    % else
    %     funPsd(sInd).svdGram = [];
    % end
    %
    %
    % if ~skip.trialGram
    %     if ~skip.psd
    %         funPsd(sInd).psdTrialGram.f = permute(funPsd(sInd).psdTrialGram.f,[2 1 3 4 5 6]);
    %         funPsd(sInd).psdTrialGram.tp = TP.trialGram;
    %         funPsd(sInd).psdTrialGram.vec;
    %         funPsd(sInd).psdTrialGram.lWin = param.psdTrialGram.win(1);
    %         funPsd(sInd).psdTrialGram.K = param.psdTrialGram.tapers(2);
    %         [~,funPsd(sInd).psdTrialGram.W,~] = K2W(funPsd(sInd).psdTrialGram.T,funPsd(sInd).psdTrialGram.K,0);
    %         funPsd(sInd).psdTrialGram.param = param.psdTrialGram;
    %         funPsd(sInd).psdTrialGram.info = strjoin({'freq/time' 'tapers' 'vox' 'run' 'timeWindow' 'mode'},' x ');
    %         funPsd(sInd).psdTrialGram.mask = funTs(sInd).vol2vec;
    %     else
    %         funPsd(sInd).psdTrialGram = [];
    %     end
    %     if ~skip.svd
    %         funPsd(sInd).svdTrialGram.f;
    %         funPsd(sInd).svdTrialGram.tp = permute(TP.trialGram,[3 4 1 5 6 2]);
    %         % funPsd(sInd).svdTrialGram.u;
    %         % funPsd(sInd).svdTrialGram.s = permute(funPsd(sInd).svdTrialGram.s,[2 1 3 4 5 6]);
    %         % funPsd(sInd).svdTrialGram.v = permute(funPsd(sInd).svdTrialGram.v,[6 2 3 4 5 1]);
    %         funPsd(sInd).svdTrialGram.coh = permute(funPsd(sInd).svdTrialGram.coh,[2 1 3 4 5 6]);
    %         funPsd(sInd).svdTrialGram.lWin = param.svdTrialGram.win(1);
    %         funPsd(sInd).svdTrialGram.K = param.svdTrialGram.tapers(2);
    %         [~,funPsd(sInd).svdTrialGram.W,~] = K2W(funPsd(sInd).svdTrialGram.T,funPsd(sInd).svdTrialGram.K,0);
    %         funPsd(sInd).svdTrialGram.param = param.svdTrialGram;
    %         funPsd(sInd).svdTrialGram.info = strjoin({'vox' 'mode' 'freq/time' 'run' 'timeWindow' 'tapers'},' x ');
    %         funPsd(sInd).svdTrialGram.mask = funTs(sInd).vol2vec;
    %     else
    %         funPsd(sInd).svdTrialGram = [];
    %     end
    % end
    %
    % %% For backward compatibility
    % funPsd(sInd).vec = permute(funPsd(sInd).psd.vec,[1 3 2 4]); funPsd(sInd).psd.vec = [];
    % funPsd(sInd).f = permute(funPsd(sInd).psd.f,[2 3 4 5 6 1]);
    % funPsd(sInd).K = funPsd(sInd).psd.K;
    % funPsd(sInd).T = funPsd(sInd).psd.T;
    % funPsd(sInd).W = funPsd(sInd).psd.W;
    % funPsd(sInd).tr = mode(diff(funPsd(sInd).f))*1000;
    % funPsd(sInd).nfreq = size(funPsd(sInd).vec,1);
    % funPsd(sInd).nframes=size(funTs(sInd).vec,1);
    % funPsd(sInd).nruns=size(funPsd(sInd).vec,4);
    % % if funPsd(sInd).nframes~=size(funTs(sInd).vec,1); dbstack; warning('something wrong with nframes'); end
    % % if funPsd(sInd).nruns~=size(funPsd(sInd).vec,4); dbstack; warning('something wrong with nframes'); end
    % funPsd(sInd).ntapers = size(funPsd(sInd).vec,3);
end


function [psd,psdGram,psdTrialGram,psdTrialGramMD,svd,svdGram,svdTrialGram,svdTrialGramMD] = computeAll(funTs,TP,param,skip,verbose)
if ~exist('nShuf','var');         nShuf = []; end
if ~exist('nRun','var');           nRun = []; end
if ~exist('cohFrange','var'); cohFrange = []; end
if isempty(nShuf);         nShuf = 0; end
if isempty(nRun);           nRun = 1; end
if isempty(cohFrange); cohFrange = [0 inf]; end
orderWin = -1;
% skip.gram = any(ismember(param.win(2),[0 inf nan]));
for runInd = 1:nRun
    if verbose && nRun>1; disp(['---Run ' num2str(runInd) '/' num2str(nRun) '---']); end

    Fs = param.Fs;
    

    %% %%%%%%%%%%%%%%%%%%%%%
    % Over full timeseires %
    %%%%%%%%%%%%%%%%%%%%% %%
    %[time x trial x run x taper x freq x vox x window x mode]

    %%% tapers
    tp = permute(TP.full.tp,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
    [Nk,Ek,Rk,Kk,Fk,Vk,Wk,M] = size(tp);
    K = Kk;
    N = Nk;

    %%% windows and trials
    W = 1;
    E = 1;

    %%% time
    if isfield(funTs,'t') && ~isempty(funTs.t)
        t = funTs.t; % tapers[time x trial x run x taper x freq x vox x window x mode]
    else
        dbstack; error('X');
        t = permute(linspace(0,(N-1)/Fs,N),[2 1 3 4 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
    end
    [Nt,Et,Rt,Kt,Ft,Vt,Wt,Mt] = size(t);

    %%% freq
    pad = 0;
    NFFT=max(2^(nextpow2(N)+pad),N);
    [f,fInd]=getfgrid(Fs,NFFT,[0 Fs/2]);
    f = permute(f,[1 3 4 5 2 6 7 8]); % frequencies[time x trial x run x taper x freq x vox x window x mode]
    [Nf,Ef,Rf,Kf,Ff,Vf,Wf] = size(f);
    F = Ff;

    %%% channels and runs
    [~,V,~,R] = size(funTs.vec);

    %%% modes
    M = min([V K]);

    %%% allocate
    res.full.PSD    = zeros(1,E,1,1,F,V,W,1  ); % psd       [time x trial x run x taper x freq x vox x window x mode] at each trial
    res.full.COH    = zeros(1,E,1,1,F,1,W,M  ); % coherence [time x trial x run x taper x freq x vox x window x mode] at each trial

    %%% Compute J
    d = funTs.vec;              % [time  x vox x taper x run               ]
    d = permute(d,[3 2 4 1]);   % [taper x vox x run   x time*trial        ]
    d = reshape(d,[1 V R N E]); % [taper x vox x run   x time       x trial]
    d = permute(d,[4 5 3 1 6 2 7 8]);
    tp = tp;
    t = t;
    J = getJ2(d,tp,t,f)/Fs; % [time x trial x run x taper x freq x vox x window]
    %                         [N      E       R     K       F      V     W]

    %%% Compute psd
    res.full.PSD = mean(  conj(J).*J  ,4);


    %%% Compute coherence
    if K > 1
        dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
        prm = [ 6   4   2   1   3   5   7   8  ];
        dim = strjoin(dim(prm),' ');
        j = permute(J,prm); %[V K E N R F W Mk]
        j = reshape(j,[V K E*1*R*F*1*1]); %[V K E*N*R*F*W*Mk]
        [u,s,v] = pagesvd(j,'econ','vector'); % s[Mk V E*N*R*F*W]
        coh     = s.^2./sum(s.^2,1); % coherence[Mk V E*N*R*F*W]
        spSVmag = abs(u);
        % coh = reshape(coh,[M 1 E 1 R F 1]); % [Mk V E N R F W]
        %       1    2   3   4   5   6   7   8
        dim = {'Mk' 'V' 'E' 'N' 'R' 'F' 'W' 'K'};
        prm = [4 3 5 8 6 2 7 1];
        % prm = [ 2    4   3   6   5   7   8   1 ];
        dim = strjoin(dim(prm),' ');
        res.full.COH  = permute(reshape(coh,[M 1 E 1 R F 1]),prm); %[time x trial x run x taper x freq x vox x window x mode]
        res.full.spSV = permute(u,[4 5 6 7 3 1 8 2]); %[time x trial x run x taper x freq x vox x window x mode]

        % %%% Estimate null distribution2
        % if param.perm
        %     allPerm = perms(1:K); allPermN = size(allPerm,1);
        %     j2 = reshape(permute(j,[3 2 1]),[F K*V]);
        % 
        %     coh_pVal = zeros(M,1,F);
        %     spSVmag_pVal = zeros(V,M,F);
        %     parfor perm = 1:param.perm
        % 
        %         [u,s,vPerm] = pagesvd(...
        %             permute(  reshape(  j2(:,curPerm)  ,[F K V])  ,[3 2 1])...
        %             ,'econ','vector'); % s[M V E*N*R*F*W]
        %         cohPerm     = s.^2./sum(s.^2,1); % coherence[M V E*N*R*F*W]
        %         spSVmagPerm = abs(u);
        %         coh_pVal     = coh_pVal + cohPerm>coh;
        %         spSVmag_pVal = spSVmag_pVal + spSVmagPerm>spSVmag;
        %     end
        % 
        % 
        % end
        
        
        %%% Estimate null distribution
        if param.perm
            disp(['doing ' num2str(param.perm) ' random taper permutations'])
            % V = 4; F = 3; K = 10;
            % j = repmat(1:K,[V 1 F]);
            %%%get all possible permutions beforehand for speed
            allPerm = perms(1:K); allPermN = size(allPerm,1);
            %%%initialize inputs and outputs
            j2 = reshape(permute(j,[3 2 1]),[F K*V]); % F K*V
            COH = permute(res.full.COH,[8 6 5 1 2 3 4 7]); % M V F
            spSVmag = permute(abs(res.full.spSV),[6 8 5 1 2 3 4 7]); % V M F
            COH_permMean     = zeros(size(COH));
            spSVmag_permMean = zeros(size(spSVmag));
            COH_pVal     = zeros(size(COH));
            spSVmag_pVal = zeros(size(spSVmag));
            %%%loop over permutations
            parfor perm = 1:param.perm
                % for each voxel, randomly permute tapers, using the same
                % permutation across frequencies
                [COH_pVal_cur,spSVmag_pVal_cur,COH_permMean_cur,spSVmag_permMean_cur] = tmpSvd(j2,allPerm,allPermN,V,K,F,COH,spSVmag);
                COH_pVal     = COH_pVal     + COH_pVal_cur;
                spSVmag_pVal = spSVmag_pVal + spSVmag_pVal_cur;
                COH_permMean     = COH_permMean     + COH_permMean_cur;
                spSVmag_permMean = spSVmag_permMean + spSVmag_permMean_cur;
            end
            %%%summarize across permutations
            COH_permMean     = COH_permMean     ./ param.perm; % M V F
            spSVmag_permMean = spSVmag_permMean ./ param.perm; % V M F
            COH_pVal     = COH_pVal     ./ param.perm; % M V F
            spSVmag_pVal = spSVmag_pVal ./ param.perm; % V M F
            %%%FDR
            COH_pVal = permute(COH_pVal,[3 2 1]); % F V M
            COH_fdr = zeros(size(COH_pVal));
            for i = 1:length(COH_pVal(1,:))
                COH_fdr(:,i) = mafdr(COH_pVal(:,i),'BHFDR',true);
            end
            COH_pVal = permute(COH_pVal,[3 2 1]); % M V F
            COH_fdr  = permute(COH_fdr ,[3 2 1]); % M V F
            
            spSVmag_fdr = zeros(size(spSVmag_pVal)); % V M F
            for i = 1:length(spSVmag_pVal(1,:))
                spSVmag_fdr(:,i) = mafdr(spSVmag_pVal(:,i),'BHFDR',true);
            end
            
            %%%compile results
            res.full.COH_permMean = permute(COH_permMean,[4 5 6 7 3 2 8 1]); %[time x trial x run x taper x freq x vox x window x mode]
            res.full.COH_pVal     = permute(COH_pVal    ,[4 5 6 7 3 2 8 1]); %[time x trial x run x taper x freq x vox x window x mode]
            res.full.COH_fdr      = permute(COH_fdr     ,[4 5 6 7 3 2 8 1]); %[time x trial x run x taper x freq x vox x window x mode]
            res.full.spSV_permMean = permute(spSVmag_permMean,[4 5 6 7 3 1 8 2]); %[time x trial x run x taper x freq x vox x window x mode]
            res.full.spSV_pVal     = permute(spSVmag_pVal    ,[4 5 6 7 3 1 8 2]); %[time x trial x run x taper x freq x vox x window x mode]
            res.full.spSV_fdr      = permute(spSVmag_fdr     ,[4 5 6 7 3 1 8 2]); %[time x trial x run x taper x freq x vox x window x mode]
        end
    else
        res.full.COH  = [];
        res.full.spSV = [];
    end

    %%%
    res.full.f      = f;
    res.full.K      = K;
    res.full.T      = N/Fs;
    res.full.E      = [];
    res.full.win    = [];
    res.full.param  = rmfield(param,{'complex' 'psd' 'svd' 'psdGram' 'svdGram' 'psdTrialGram' 'svdTrialGram'});
    res.full.dim     = [N E R K F V W M];
    res.full.dimInfo = '[N E R K F V W M]';



    %% %%%%%%%%%%%%%%%%%%%%%
    % Over each timewindow %
    %%%%%%%%%%%%%%%%%%%%% %%
    if ~skip.gram
        %[time x trial x run x taper x freq x vox x window x mode]
        %[   7       2     1       3      5     8       20]
        tp = permute(TP.gram.tp,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
        [Nk,Ek,Rk,Kk,Fk,Vk,Wk,M] = size(tp);
        K = Kk;

        %%% windows
        w = permute(param.psdGram.winInd,[2 3 4 5 6 7 1 8]); % windows[time x trial x run x taper x freq x vox x window x mode]
        [Nw,Ew,Rw,Kw,Fw,Vw,Ww,Mw] = size(w);
        W = Ww;
        E = Ew;
        if Nw~=Nk; dbstack; error('X'); end
        N = Nk;

        %%% time
        t = permute(linspace(0,(N-1)/Fs,N),[2 1 3 4 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
        [Nt,Et,Rt,Kt,Ft,Vt,Wt,Mt] = size(t);

        %%% freq
        pad = 0;
        NFFT=max(2^(nextpow2(N)+pad),N);
        [f,fInd]=getfgrid(Fs,NFFT,[0 Fs/2]);
        f = permute(f,[1 3 4 5 2 6 7 8]); % frequencies[time x trial x run x taper x freq x vox x window x mode]
        [Nf,Ef,Rf,Kf,Ff,Vf,Wf] = size(f);
        F = Ff;

        %%% channels and runs
        [~,V,~,R] = size(funTs.vec);

        %%% modes
        M = min([V K]);

        %%% allocate
        res.gram.PSD    = zeros(1,E,1,1,F,V,W,1  ); % psd       [time x trial x run x taper x freq x vox x window x mode] at each trial
        res.gram.COH    = zeros(1,E,1,1,F,1,W,M  ); % coherence [time x trial x run x taper x freq x vox x window x mode] at each trial

        %%% loop over windows
        for wInd = 1:W
            %%% Compute J
            ind  = w(:,:,:,:,:,:,wInd);
            tWin = reshape(funTs.t(ind,:,:,:),size(ind)); tWin = tWin([1 end],:);
            d = funTs.vec(ind,:,:,:);   % [time  x vox x taper x run               ]
            d = permute(d,[3 2 4 1]);   % [taper x vox x run   x time*trial        ]
            d = reshape(d,[1 V R N E]); % [taper x vox x run   x time       x trial]
            d = permute(d,[4 5 3 1 6 2 7 8]); %[time x trial x run x taper x freq x vox x window x mode]
            if orderWin == 0
                d = d - mean(d,1);
            end
            tp = tp;
            t = t; % adjust t here for sub-tr stimulus onsets
            J = getJ2(d,tp,t,f)/Fs; % [time x trial x run x taper x freq x vox x window]
            %                         [N      E       R     K       F      V     W]

            %%% Compute psd at each trial
            res.gram.PSD(:,:,:,:,:,:,wInd)    = mean(  conj(J).*J  ,4);

            %%% Compute coherence at each trial
            if K > 1
                dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
                prm = [ 6   4   2   1   3   5   7   8  ];
                dim = strjoin(dim(prm),' ');
                j = permute(J,prm); %[V K E N R F W Mk]
                j = reshape(j,[V K E*1*R*F*1*1]); %[V K E*N*R*F*W*Mk]
                [u,s,v] = pagesvd(j,'econ','vector'); % s[Mk V E*N*R*F*W]
                coh = s.^2./sum(s.^2,1); % coherence[Mk V E*N*R*F*W]
                coh = reshape(coh,[M 1 E 1 R F 1]); % [Mk V E N R F W]
                dim = {'Mk' 'V' 'E' 'N' 'R' 'F' 'W' 'K'};
                prm = [ 2    4   3   6   5   7   8   1 ];
                dim = strjoin(dim(prm),' ');
                coh = permute(coh,prm); %[V N E F R W K Mk]
                res.gram.COH(:,:,:,:,:,:,wInd,:) = coh;
            end

            %%% Output window time
            res.gram.t(:,:,:,:,:,:,wInd,:) = tWin;
        end
        if K == 1
            res.gram.COH = [];
        end


        res.gram.f = f;
        res.gram.K = K;
        res.gram.T = N/Fs;
        res.gram.E = [];
        res.gram.win = [mean(reshape(diff(res.gram.t,[],1),[],1)) mean(reshape(diff(res.gram.t,[],7),[],1))];
        res.gram.param  = rmfield(param,{'complex' 'psd' 'svd' 'psdGram' 'svdGram' 'psdTrialGram' 'svdTrialGram'});
    else
        res.gram = [];
    end



    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Over each event-related timewindow %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %[time x trial x run x taper x freq x vox x window x mode]
    %[   7       2     1       3      5     8       20]
    if ~skip.trialGram
        %%% tapers
        tp = permute(TP.gram.tp,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
        [Nk,Ek,Rk,Kk,Fk,Vk,Wk,M] = size(tp);
        K = Kk;

        %%% windows and trials
        w = permute(param.psdTrialGram.winInd,[2 3 4 5 6 7 1 8]); % windows[time x trial x run x taper x freq x vox x window x mode]
        [Nw,Ew,Rw,Kw,Fw,Vw,Ww,Mw] = size(w);
        W = Ww;
        E = Ew;
        if Nw~=Nk; dbstack; error('X'); end
        N = Nk;


        %%% time
        t = permute(linspace(0,(N-1)/Fs,N),[2 1 3 4 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
        [Nt,Et,Rt,Kt,Ft,Vt,Wt,Mt] = size(t);

        %%% freq
        pad = 0;
        % pad = E-1;
        NFFT=max(2^(nextpow2(N*E)+pad),N*E);
        [f,fInd]=getfgrid(Fs,NFFT,[0 Fs/2]);
        f = permute(f,[1 3 4 5 2 6 7 8]); % frequencies[time x trial x run x taper x freq x vox x window x mode]
        [Nf,Ef,Rf,Kf,Ff,Vf,Wf] = size(f);
        F = Ff;


        %%% channels and runs
        [~,V,~,R] = size(funTs.vec);

        %%% modes
        M = min([V K]);
        Mek = min([V K*E]);

        %%% allocate
        res.trialGram.PSD    = zeros(1,E,1,1,F,V,W,1  ); % psd       [time x trial x run x taper x freq x vox x window x mode] at each trial
        res.trialGram.PSDeav = zeros(1,1,1,1,F,V,W,1  ); % psd       [time x trial x run x taper x freq x vox x window x mode] averaged across trials                      (eVENT avRAGED       )
        res.trialGram.PSDepc = zeros(1,1,1,1,F,V,W,1  ); % psd       [time x trial x run x taper x freq x vox x window x mode] phase-coherently averaged across trials     (eVENT pHASE cOHERENT)
        res.trialGram.COH    = zeros(1,E,1,1,F,1,W,M  ); % coherence [time x trial x run x taper x freq x vox x window x mode] at each trial
        res.trialGram.COHeav = zeros(1,1,1,1,F,1,W,M  ); % coherence [time x trial x run x taper x freq x vox x window x mode] averaged across trials                      (eVENT avRAGED       )
        res.trialGram.COHepc = zeros(1,1,1,1,F,1,W,M  ); % coherence [time x trial x run x taper x freq x vox x window x mode] phase-coherently averaged across trials     (eVENT pHASE cOHERENT)
        res.trialGram.COHek  = zeros(1,1,1,1,F,1,W,Mek); % coherence [time x trial x run x taper x freq x vox x window x mode] trials concatenated as extra sets of tapers (eVENT AS TAPERS k   )


        %%% loop over windows
        for wInd = 1:W
            %%% Compute J
            ind  = w(:,:,:,:,:,:,wInd);
            tWin = reshape(funTs.t(ind,:,:,:),size(ind)); tWin = tWin([1 end],:);
            d = funTs.vec(ind,:,:,:);   % [time  x vox x taper x run               ]
            d = permute(d,[3 2 4 1]);   % [taper x vox x run   x time*trial        ]
            d = reshape(d,[1 V R N E]); % [taper x vox x run   x time       x trial]
            d = permute(d,[4 5 3 1 6 2 7 8]); %[time x trial x run x taper x freq x vox x window x mode]
            if orderWin == 0
                d = d - mean(d,1);
            end
            tp = tp;
            t = t; % adjust t here for sub-tr stimulus onsets
            J = getJ2(d,tp,t,f)/Fs; % [time x trial x run x taper x freq x vox x window]
            %                         [N      E       R     K       F      V     W]

            %%% Compute psd at each trial
            res.trialGram.PSD(:,:,:,:,:,:,wInd)    = mean(  conj(J).*J  ,4);
            %%% Compute psd averaged across trials
            res.trialGram.PSDeav(:,:,:,:,:,:,wInd)   = mean(  res.trialGram.PSD(:,:,:,:,:,:,wInd)  ,2);
            %%% Compute psd phase-coherently averaged across trials
            res.trialGram.PSDepc(:,:,:,:,:,:,wInd) = mean(  conj(mean(J,2)).*mean(J,2)  ,4);


            if K > 1
                %%% Compute coherence at each trial
                dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
                prm = [ 6   4   2   1   3   5   7   8  ];
                dim = strjoin(dim(prm),' ');
                j = permute(J,prm); %[V K E N R F W Mk]
                j = reshape(j,[V K E*1*R*F*1*1]); %[V K E*N*R*F*W*Mk]
                [u,s,v] = pagesvd(j,'econ','vector'); % s[M V E*N*R*F*W]
                coh = s.^2./sum(s.^2,1); % coherence[Mk V E*N*R*F*W]
                coh = reshape(coh,[M 1 E 1 R F 1]); % [M V E N R F W]
                %       1   2   3   4   5   6   7   8
                dim = {'M' 'V' 'E' 'N' 'R' 'F' 'W' 'K'};
                prm = [4 3 5 8 6 2 7 1];
                     %[N E R K F V W M]
                dim = strjoin(dim(prm),' ');
                coh = permute(coh,prm); %[N E R K F V W M]
                res.trialGram.COH(:,:,:,:,:,:,wInd,:) = coh;

                %%% Average coherence across trials
                res.trialGram.COHeav(:,:,:,:,:,:,wInd,:) = mean(res.trialGram.COH(:,:,:,:,:,:,wInd,:),2);

                %%% Compute coherence after phase-coherently averaging across trials
                dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
                prm = [ 6   4   2   1   3   5   7   8  ];
                dim = strjoin(dim(prm),' ');
                j = permute(mean(J,2),prm); %[V K E N R F W Mk]
                j = reshape(j,[V K 1*1*R*F*1*1]); %[V K E*N*R*F*W*Mk]
                [u,s,v] = pagesvd(j,'econ','vector'); % s[Mk V E*N*R*F*W]
                coh = s.^2./sum(s.^2,1); % coherence[Mk V E*N*R*F*W]
                coh = reshape(coh,[M 1 1 1 R F 1]); % [Mk V E N R F W]
                dim = {'Mk' 'V' 'E' 'N' 'R' 'F' 'W' 'K'};
                prm = [ 2    4   3   6   5   7   8   1 ];
                dim = strjoin(dim(prm),' ');
                coh = permute(coh,prm); %[V N E F R W K Mk]
                res.trialGram.COHepc(:,:,:,:,:,:,wInd,:) = coh;
            end

            %%% Compute coherence with trials concatenated as extra sets of tapers (equivalent to averaging across trials)
            dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
            prm = [ 6   4   2   1   3   5   7   8  ];
            dim = strjoin(dim(prm),' ');
            j = permute(J,prm); %[V K E N R F W Mk]
            j = reshape(j,[V K E 1*R*F*1*1]); %[V K E N*R*F*W*Mk]
            dim = {'V' 'K' 'E' 'N*R*F*W*Mk'};
            prm = [ 1   4   2   3          ];
            dim = strjoin(dim(prm),' ');
            j = permute(j,prm); %[V N*R*F*W*Mk K E]
            j = reshape(j,[V 1*R*F*1*1 K*E]); %[V N*R*F*W*Mk K*E]
            j = permute(j,[1 3 2]); %[V K*E N*R*F*W*Mk];
            [u,s,v] = pagesvd(j,'econ','vector'); % s[Mke V N*R*F*W K E]
            coh = s.^2./sum(s.^2,1); % coherence[Mke V N*R*F*W K E]
            coh = reshape(coh,[Mek 1 1 1 F 1 1 1]); % [Mke V N R F W K E]
            %       1   2   3   4   5   6   7   8
            dim = {'Mke' 'V' 'N' 'R' 'F' 'W' 'K' 'E'};
            prm = [3 8 4 7 5 2 6 1];
                 %[N E R K F V W M]
            dim = strjoin(dim(prm),' ');
            coh = permute(coh,prm); %[V N E F R W K Mke]
            res.trialGram.COHek(:,:,:,:,:,:,wInd,:) = coh;

            %%% Output window time
            res.trialGram.t(:,:,:,:,:,:,wInd,:) = tWin;
        end
        if K == 1
            res.trialGram.COH    = [];
            res.trialGram.COHeav = [];
            res.trialGram.COHepc = [];
        end

        res.trialGram.f         = f;
        res.trialGram.K         = K;
        res.trialGram.T         = N/Fs;
        res.trialGram.E         = E;
        res.trialGram.win       = [mean(reshape(diff(res.trialGram.t,[],1),[],1)) mean(reshape(diff(res.trialGram.t,[],7),[],1))];
        res.trialGram.onsetList = permute(param.onsetList,[2 1 3 4 5 6 7 8]);
        res.trialGram.durList   = permute(param.durList  ,[2 1 3 4 5 6 7 8]);
        res.trialGram.param     = rmfield(param,{'complex' 'psd' 'svd' 'psdGram' 'svdGram' 'psdTrialGram' 'svdTrialGram'});
    else
        res.trialGram = [];
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Over each event-related timewindow   %
    % (with missing data seperating tials) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %[time x trial x run x taper x freq x vox x window x mode]
    %[   7       2     1       3      5     8       20]
    if ~skip.trialGram
        %%% tapers
        tp = permute(TP.trialGramMD.tp,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
        [Nk,Ek,Rk,Kk,Fk,Vk,Wk,M] = size(tp);
        K = Kk;

        %%% windows and trials
        w = permute(param.psdTrialGram.winInd,[2 3 4 5 6 7 1 8]); % windows[time x trial x run x taper x freq x vox x window x mode]
        [Nw,Ew,Rw,Kw,Fw,Vw,Ww,Mw] = size(w);
        W = Ww;
        E = Ew;
        if Nw*Ew~=Nk; dbstack; error('X'); end
        N = Nk;

        %%% time
        %%%% actual time
        t = TP.trialGramMD.t;
        [Nt,Et,Rt,Kt,Ft,Vt,Wt,Mt] = size(t);
        %%%% shifted to 0 at stim onset (phase coherent cross-trial avg)
        tPC = reshape(reshape(t,Nw,E) - funTs.dsgn.onsets',Nw*E,1);
        

        %%% freq
        pad = 0;
        NFFT=max(2^(nextpow2(N)+pad),N);
        [f,fInd]=getfgrid(Fs,NFFT,[0 Fs/2]);
        f = permute(f,[1 3 4 5 2 6 7 8]); % frequencies[time x trial x run x taper x freq x vox x window x mode]
        [Nf,Ef,Rf,Kf,Ff,Vf,Wf] = size(f);
        F = Ff;

        %%% channels and runs
        [~,V,~,R] = size(funTs.vec);

        %%% modes
        M = min([V K]);
        
        %%% allocate
        res.trialGramMD.PSD    = zeros(1,1,1,1,F,V,W,1); % psd       [time x trial x run x taper x freq x vox x window x mode] phase-coherently averaged across trials     (eVENT pHASE cOHERENT)
        res.trialGramMD.PSDepc = zeros(1,1,1,1,F,V,W,1); % psd       [time x trial x run x taper x freq x vox x window x mode] phase-coherently averaged across trials     (eVENT pHASE cOHERENT)
        res.trialGramMD.COH    = zeros(1,1,1,1,F,1,W,M); % coherence [time x trial x run x taper x freq x vox x window x mode] trials concatenated as extra sets of tapers (eVENT AS TAPERS k   )
        res.trialGramMD.COHepc = zeros(1,1,1,1,F,1,W,M); % coherence [time x trial x run x taper x freq x vox x window x mode] trials concatenated as extra sets of tapers (eVENT AS TAPERS k   )

        %%% loop over windows
        for wInd = 1:W
            %%% Compute J
            ind  = w(:,:,:,:,:,:,wInd);
            tWin = reshape(funTs.t(ind,:,:,:),size(ind)); tWin = tWin([1 end],:);
            d = permute(funTs.vec(ind,:,:,:),[2 3 4 1]); % [vox x taper x run x time]
            d = reshape(d,[V 1 R Nw E]); % [vox x taper x run x time x trial]
            if orderWin == 0
                d = d - mean(d,4);
            end
            d = reshape(d,[V 1 R Nw*E]); % [vox x taper x run x time]
            d = permute(d,[4 5 3 2 6 1 7 8]); % [time x trial x run x taper x freq x vox x window]
            % d = funTs.vec(ind,:,:,:);         % [time  x vox x taper x run               x window]
            % d = permute(d,[1 5 4 3 6 2 7 8]); % [time x trial x run x taper x freq x vox x window]
            tp = tp;

            
            %%% Using normal time (incoherent phase averaging across trials)
            t = t; % adjust t here phase coherent cross-trial averaging
            J = getJ2(d,tp,t,f)/Fs; % [time x trial x run x taper x freq x vox x window]
            %                         [N      E       R     K       F      V     W]

            %%%% Compute psd
            res.trialGramMD.PSD(:,:,:,:,:,:,wInd)    = mean(  conj(J).*J  ,4);
            
            %%%% Compute coherence
            if K > 1
                dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
                prm = [ 6   4   2   1   3   5   7   8  ];
                dim = strjoin(dim(prm),' ');
                j = permute(J,prm); %[V K E N R F W M]
                j = reshape(j,[V K 1*1*R*F*1*1]); %[V K E*N*R*F*W*Mk]
                [u,s,v] = pagesvd(j,'econ','vector'); % s[M V E*N*R*F*W]
                coh = s.^2./sum(s.^2,1); % coherence[M V E*N*R*F*W]
                coh = reshape(coh,[M 1 1 1 R F 1]); % [M V E N R F W]
                %       1   2   3   4   5   6   7   8
                dim = {'M' 'V' 'E' 'N' 'R' 'F' 'W' 'K'};
                prm = [4 3 5 8 6 2 7 1];
                %     [N E R K F V W M]
                dim = strjoin(dim(prm),' ');
                coh = permute(coh,prm); %[N E R K F V W M]
                res.trialGramMD.COH(:,:,:,:,:,:,wInd,:) = coh;
            end


            %%% Now adjusting time for phase coherent cross-trial averaging
            tPC = tPC;
            J = getJ2(d,tp,tPC,f)/Fs; % [time x trial x run x taper x freq x vox x window]
            %                         [N      E       R     K       F      V     W]

            %%%% Compute psd
            res.trialGramMD.PSDepc(:,:,:,:,:,:,wInd)    = mean(  conj(J).*J  ,4);
            
            %%%% Compute coherence
            if K > 1
                dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
                prm = [ 6   4   2   1   3   5   7   8  ];
                dim = strjoin(dim(prm),' ');
                j = permute(J,prm); %[V K E N R F W M]
                j = reshape(j,[V K 1*1*R*F*1*1]); %[V K E*N*R*F*W*Mk]
                [u,s,v] = pagesvd(j,'econ','vector'); % s[M V E*N*R*F*W]
                coh = s.^2./sum(s.^2,1); % coherence[M V E*N*R*F*W]
                coh = reshape(coh,[M 1 1 1 R F 1]); % [M V E N R F W]
                %       1   2   3   4   5   6   7   8
                dim = {'M' 'V' 'E' 'N' 'R' 'F' 'W' 'K'};
                prm = [4 3 5 8 6 2 7 1];
                %     [N E R K F V W M]
                dim = strjoin(dim(prm),' ');
                coh = permute(coh,prm); %[V N E F R W K M]
                res.trialGramMD.COHepc(:,:,:,:,:,:,wInd,:) = coh;
            end

            %%% Output window time
            res.trialGramMD.t(:,:,:,:,:,:,wInd,:) = tWin;
        end
        if K == 1
            res.trialGramMD.COH    = [];
            res.trialGramMD.COHepc = [];
        end

        res.trialGramMD.f         = f;
        res.trialGramMD.K         = K;
        res.trialGramMD.T         = N/Fs;
        res.trialGramMD.E         = E;
        res.trialGramMD.win       = [mean(reshape(diff(res.trialGram.t,[],1),[],1)) mean(reshape(diff(res.trialGram.t,[],7),[],1))];
        res.trialGramMD.onsetList = permute(param.onsetList,[2 1 3 4 5 6 7 8]);
        res.trialGramMD.durList   = permute(param.durList  ,[2 1 3 4 5 6 7 8]);
        res.trialGramMD.param     = rmfield(param,{'complex' 'psd' 'svd' 'psdGram' 'svdGram' 'psdTrialGram' 'svdTrialGram'});
    else
        res.trialGramMD = [];
    end




end

if skip.psd
    psd = [];
else
    fields = {'COH' 'COH_pVal' 'COH_fdr' 'spSV' 'spSV_pVal' 'spSV_fdr'}; fields = fields(ismember(fields,fieldnames(res.full))); 
    psd = rmfield(res.full,fields);
    psd.info = 'time x trial x run x taper x freq x vox x window x mode';
end

if skip.gram || skip.psd
    psdGram = [];
else
    fields = {'COH' 'COH_pVal' 'COH_fdr' 'spSV' 'spSV_pVal' 'spSV_fdr'}; fields = fields(ismember(fields,fieldnames(res.gram)));
    psdGram = rmfield(res.gram,fields);
    psdGram.info = 'time x trial x run x taper x freq x vox x window x mode';
end

if skip.trialGram || skip.psd
    psdTrialGram = [];
    psdTrialGramMD = [];
else
    fields = {'PSD' 'PSDeav' 'PSDepc' 'COH' 'COHeav' 'COHepc' 'COHek'}; fields = fields(ismember(fields,fieldnames(res.trialGram)));
    psdTrialGram = rmfield(res.trialGram,fields);
    psdTrialGram.vec.psd   = res.trialGram.PSDeav;
    psdTrialGram.vec.psdPC = res.trialGram.PSDepc;
    psdTrialGram.info = 'time x trial x run x taper x freq x vox x window x mode';

    fields = {'PSD' 'PSDepc' 'COH' 'COHepc'}; fields = fields(ismember(fields,fieldnames(res.trialGramMD)));
    psdTrialGramMD = rmfield(res.trialGramMD,fields);
    psdTrialGramMD.vec.psd   = res.trialGramMD.PSD;
    psdTrialGramMD.vec.psdPC = res.trialGramMD.PSDepc;
    psdTrialGramMD.info = 'time x trial x run x taper x freq x vox x window x mode';
end

if skip.svd
    svd = [];
else
    fields = {'PSD'}; fields = fields(ismember(fields,fieldnames(res.full)));
    svd = rmfield(res.full,fields);
    svd.info = 'time x trial x run x taper x freq x vox x window x mode';
end

if skip.svd || skip.gram
    svdGram = [];
else
    fields = {'PSD'}; fields = fields(ismember(fields,fieldnames(res.gram)));
    svdGram = rmfield(res.gram,fields);
    svdGram.info = 'time x trial x run x taper x freq x vox x window x mode';
end

if skip.trialGram || skip.svd
    svdTrialGram = [];
    svdTrialGramMD = [];
else
    fields = {'PSD' 'PSDeav' 'PSDepc' 'COH' 'COHeav' 'COHepc' 'COHek'}; fields = fields(ismember(fields,fieldnames(res.trialGram)));
    svdTrialGram = rmfield(res.trialGram,fields);
    svdTrialGram.vec.coh    = res.trialGram.COHeav;
    svdTrialGram.vec.cohEPC = res.trialGram.COHepc;
    svdTrialGram.vec.cohEK  = res.trialGram.COHek;
    svdTrialGram.info = 'time x trial x run x taper x freq x vox x window x mode';

    fields = {'PSD' 'PSDepc' 'COH' 'COHepc'}; fields = fields(ismember(fields,fieldnames(res.trialGramMD)));
    svdTrialGramMD = rmfield(res.trialGramMD,fields);
    svdTrialGramMD.vec.coh    = res.trialGramMD.COH;
    svdTrialGramMD.vec.cohEPC = res.trialGramMD.COHepc;
    svdTrialGramMD.info = 'time x trial x run x taper x freq x vox x window x mode';
end

function [COH_permAbove,spSVmag_permAbove,COH_perm,spSVmag_perm] = pagesvdPerm(j2,allPerm,allPermN,V,K,F,COH,spSVmag)
% for each voxel, randomly permute tapers, using the same
% permutation across frequencies
[uPerm,sPerm,~] = pagesvd(...
    permute(  reshape(  j2(:,allPerm(randi(allPermN,V,1),:)' + (0:K:V*K-1))  ,[F K V])  ,[3 2 1])...
    ,'econ','vector');
COH_perm          = sPerm.^2./sum(sPerm.^2,1);
COH_permAbove     = COH_perm > COH;
spSVmag_perm      = abs(uPerm);
spSVmag_permAbove = spSVmag_perm > spSVmag;
