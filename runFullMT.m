function funPsd = runFullMT(funTs,W,K,win,mask,memFlag,skipSVD,skipPSD,verbose,taperPerm,phaseRand)
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

if K==1; skipSVD = true; end

if iscell(funTs)
    for I = 1:numel(funTs)
        funPsd{I} = runFullMT(funTs{I},W,K,win,mask,memFlag,skipSVD,skipPSD,verbose,taperPerm,phaseRand);
    end
elseif isstruct(funTs)
    for I = 1:numel(funTs)
        cohFperm = [0 1.1];
        funPsd(I) = doIt(funTs(I),W,K,win,mask,memFlag,skipSVD,skipPSD,verbose,taperPerm,phaseRand,[],cohFperm);

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





function funPsd = doIt(funTs,W,K,win,mask,memFlag,skipSVD,skipPSD,verbose,taperPerm,phaseRand,cohF,cohFperm)
if length(K)>1; tp = K; K = size(tp,2); else tp = []; end
if ~exist('win','var') || isempty(win); param.win = [inf 0]; else, param.win = win; end; clear win
if funTs.nvoxels==1; if verbose; disp('only one timeseries, skipping SVD'); end; skipSVD = true; end

if param.win(1)==inf; skipGram = true; else skipGram = false; end
if length(param.win)<=2; skipTrialGram = 1; else; skipTrialGram = 0; end

if length(param.win)>2
    param.onsetList = param.win(3:end)';
    param.win(3:end) = [];
else
    param.onsetList = [];
end


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
if length(winSz2)~=1; error('badly defined trial-locked windows'); end
allWin4 = zeros(winSz2,size(allWin3,2),size(allWin3,3));
for winInd = 1:prod(size(allWin3,[2 3]))
    allWin4(:,winInd) = find(allWin3(:,winInd));
end
allWinTrialLock = permute(allWin4,[2 1 3]); % [win X timeIndex]
% reshape(allWinTrialLock,size(allWinTrialLock,1),prod(size(allWinTrialLock,[2 3])))
clear allWin2 allWin3 allWin4 winSz2



%% initiate stuff
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
    if skipPSD
        dbstack; error('code that')
    else

        %%% full timeseries
        if verbose; disp('getting taper for full timeseries'); end
        K   = param.tapers(2);
        tr  = funTs(sInd).tr/1000;
        N   = funTs(sInd).nframes;
        TP.full = getTapers(K,tr,N);
        % [~,funPsd(sInd).psd.f,funPsd(sInd).psd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramPsd,funTs(sInd).t(:,:,:,:,:,1),tp);
        
        %%% time-resolved
        if ~skipGram
            if verbose; disp('getting taper for time-resolved analysis'); end
            K   = param.tapers(2);
            tr  = funTs(sInd).tr/1000;
            N   = size(allWin,2);
            TP.gram = getTapers(K,tr,N);
        else
            TP.gram = [];
        end

        %%% trial-locked
        if ~skipTrialGram
            if verbose; disp('getting taper for trial-locked analysis'); end
            %%% Defining KE cross-trial tapers as K regular taper per E
            %%% events (trials)
            K   = param.tapers(2);
            tr  = funTs(sInd).tr/1000;
            N   = size(allWin,2);
            TP.trialGram = getTapers(K,tr,N);
            E = length(param.onsetList);
            % eval(['TP.trialGram = blkdiag(' strjoin(repmat({'TP.trialGram'},1,E),',') ');']);
            eval(['TP.trialGram = reshape(blkdiag(' strjoin(repmat({'TP.trialGram'},1,E),',') '),[N*E K E]);']);
            % TPx = nan(size(funTs(sInd).t,1),size(TP.trialGram,2));
            % TPx(allWinTrialLock(1,:),:) = TP.trialGram;
            % plot(funTs(sInd).t,TPx)


            %%% Defining missing data tapers
            % if verbose; disp('getting taper for trial-locked analysis'); end
            % K   = param.tapers(2);
            % tr  = funTs(sInd).tr/1000;
            % N   = false(1,funTs(sInd).nframes); N(allWinTrialLock(1,:)) = true; % here N is a logical vector specifying included data, producing tapers appropriate for timeseries with missing data
            % TP.trialGram = getTapers(K,tr,N);
        else
            TP.trialGram = [];
        end
        TP.info = 'time x taper x trial';
        

        % if verbose; disp('getting taper for psd'); end
        % 
        % %%% full spectrum
        % paramPsd = param;
        % [~,funPsd(sInd).psd.f,funPsd(sInd).psd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramPsd,funTs(sInd).t(:,:,:,:,:,1),tp);
        % funPsd(sInd).psd.tp = permute(funPsd(sInd).psd.tp,[3 4 5 1 2]);
        % szPsd = size(funTs(sInd).vec);
        % 
        % %%% time-resolved spectrum
        % if ~skipGram
        %     paramPsdGram = param;
        %     [~,funPsd(sInd).psdGram.f,funPsd(sInd).psdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramPsdGram.win(1)/tr,1,:,1), paramPsdGram,funTs(sInd).t(1:paramPsdGram.win(1)/tr,:,1),tp);
        %     funPsd(sInd).psdGram.tp = permute(funPsd(sInd).psdGram.tp,[3 4 5 1 2]);
        % else
        %     paramPsdGram = [];
        % end
        % szPsdGram = size(funTs(sInd).vec);
        % 
        % 
        % %%% time-resolved trial-locked spectrum
        % if ~skipTrialGram
        %     paramPsdTrialGram = param;
        %     t   = funTs(sInd).t(allWinTrialLock(1,:));
        %     N   = length(t);
        %     K   = paramPsdTrialGram.tapers(2);
        %     Fs  = paramPsdTrialGram.Fs;
        %     T   = N .* 1/Fs;
        %     [TW,W,K] = K2W(T,K);
        %     paramPsdTrialGram.tapers = [TW K];
        %     paramPsdTrialGram.win(1) = T;
        % 
        %     [~,funPsd(sInd).psdTrialGram.tp] = MDslepian(W,K,t,Fs);
        % 
        %     funPsd(sInd).psdTrialGram.tp = permute(funPsd(sInd).psdTrialGram.tp,[3 4 5 1 2]);
        %     funPsd(sInd).psdTrialGram.f = getfgrid2(Fs,max(2^(nextpow2(N)),N),[0 Fs/2]);
        % 
        %     % % [~,funPsd(sInd).psdTrialGram.f,funPsd(sInd).psdTrialGram.tp] = mtspectrumc4(funTs(sInd).vec(allWinTrialLock(1,:),1,:,1), paramPsdTrialGram,funTs(sInd).t(allWinTrialLock(1,:),:,1),tp);
        %     % [~,funPsd(sInd).psdTrialGram.f,funPsd(sInd).psdTrialGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramPsdTrialGram.win(1)/tr,1,:,1), paramPsdTrialGram,funTs(sInd).t(1:paramPsdTrialGram.win(1)/tr,:,1),tp);
        %     % funPsd(sInd).psdTrialGram.tp = permute(funPsd(sInd).psdTrialGram.tp,[3 4 5 1 2]);
        % else
        %     paramPsdTrialGram = [];
        % end
        % szPsdTrialGram = size(funTs(sInd).vec);
        % 
        % 
        % 
        % 
        % % paramPsdGram = param;
        % % if verbose && ~isempty(tp); disp('getting taper for psd'); end
        % % if param.win(2)==inf
        % %     %%% full spectrum
        % %     paramPsd = param;
        % %     [~,funPsd(sInd).psd.f,funPsd(sInd).psd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramPsd,funTs(sInd).t(:,:,:,:,:,1),tp);
        % %     funPsd(sInd).psd.tp = permute(funPsd(sInd).psd.tp,[3 4 5 1 2]);
        % % else
        % %     %%% time-window spectrum
        % %     paramPsd = param;
        % %     % [paramPsd.tapers(1),~,paramPsd.tapers(2)] = K2W(funTs(sInd).nframes.*tr,ceil(param.tapers(2)/2),verbose);
        % %     [~,funPsd(sInd).psd.f,funPsd(sInd).psd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramPsd,funTs(sInd).t(:,:,:,:,:,1),tp);
        % %     funPsd(sInd).psd.tp = permute(funPsd(sInd).psd.tp,[3 4 5 1 2]);
        % %     %%% time-resolved spectrum
        % %     paramPsdGram = param;
        % %     % [~,funPsd(sInd).psdGram.f,funPsd(sInd).psdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramPsdGram.win(1)/tr,1,:,1), paramPsdGram,funTs(sInd).t(:,:,:,1:paramPsdGram.win(1)/tr,:,1),tp);
        % %     [~,funPsd(sInd).psdGram.f,funPsd(sInd).psdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramPsdGram.win(1)/tr,1,:,1), paramPsdGram,funTs(sInd).t(1:paramPsdGram.win(1)/tr,:,1),tp);
        % %     funPsd(sInd).psdGram.tp = permute(funPsd(sInd).psdGram.tp,[3 4 5 1 2]);
        % % end
        % % % % use time windows
        % % % [~,funPsd(sInd).psdGram.f,funPsd(sInd).psdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramPsdGram.win(1)/tr,1,:,1), paramPsdGram,funTs(sInd).t(:,:,:,1:paramPsdGram.win(1)/tr,:,1),tp);
        % % % funPsd(sInd).psdGram.tp = permute(funPsd(sInd).psdGram.tp,[3 4 5 1 2]);
        % % % if paramPsdGram.win(2)~=inf
        % % %     % prepare for also performng psd over the full timeseries,
        % % %     % using half as many tapers
        % % %     paramPsd = paramPsdGram;
        % % %     [paramPsd.tapers(1),~,paramPsd.tapers(2)] = K2W(funTs(sInd).nframes.*tr,ceil(K/2),verbose);
        % % %     [~,funPsd(sInd).psd.f,funPsd(sInd).psd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramPsd,funTs(sInd).t(:,:,:,:,:,1),tp);
        % % %     funPsd(sInd).psd.tp = permute(funPsd(sInd).psd.tp,[3 4 5 1 2]);
        % % % else
        % % %     paramPsd = param;
        % % %     funPsd(sInd).psd.tp = funPsd(sInd).psdGram.tp;
        % % % end
    end


    % if skipSVD
    %     paramSvd = [];
    % else
    %     if verbose; disp('getting taper for svd'); end
    % 
    %     %%% full coherence spectrum
    %     paramSvd = param;
    %     [~,funPsd(sInd).svd.f,funPsd(sInd).svd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramSvd,funTs(sInd).t(1,1,1,:,1,1),tp);
    %     funPsd(sInd).svd.tp = permute(funPsd(sInd).svd.tp,[2 3 1]);
    %     szSvd = size(funTs(sInd).vec);
    % 
    %     %%% time-window coherence spectra
    %     if ~skipGram
    %         paramSvdGram = param;
    %         [~,funPsd(sInd).svdGram.f,funPsd(sInd).svdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramSvdGram.win(1)/tr,1,:,1),paramSvdGram,funTs(sInd).t(1:paramSvdGram.win(1)/tr,1,1),tp);
    %         funPsd(sInd).svdGram.tp = permute(funPsd(sInd).svdGram.tp,[2 3 1]);
    %     else
    %         paramSvdGram = [];
    %     end
    %     szSvdGram = size(funTs(sInd).vec);
    % 
    %     if ~skipTrialGram
    %         if ~skipPSD
    %             paramSvdTrialGram = paramPsdTrialGram;
    %             funPsd(sInd).svdTrialGram.f  = funPsd(sInd).psdTrialGram.f;
    %             funPsd(sInd).svdTrialGram.tp = funPsd(sInd).psdTrialGram.tp;
    %         else
    %             dbstack; error('code that');
    %         end
    %     else
    %         paramSvdTrialGram = [];
    %     end
    %     szSvdTrialGram = size(funTs(sInd).vec);
    % 
    % 
    % 
    %     % if param.win(2)==inf
    %     %     %%% full coherence spectrum
    %     %     paramSvd = param;
    %     %     % [paramSvd.tapers(1),~,paramSvd.tapers(2)] = K2W(funTs(sInd).nframes.*tr,ceil(param.tapers(2)*2),verbose);
    %     %     [~,funPsd(sInd).svd.f,funPsd(sInd).svd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramSvd,funTs(sInd).t(1,1,1,:,1,1),tp);
    %     %     funPsd(sInd).svd.tp = permute(funPsd(sInd).svd.tp,[2 3 1]);
    %     % else
    %     %     %%% time-window coherence spectrum
    %     %     paramSvd = param;
    %     %     [~,funPsd(sInd).svd.f,funPsd(sInd).svd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramSvd,funTs(sInd).t(1,1,1,:,1,1),tp);
    %     %     funPsd(sInd).svd.tp = permute(funPsd(sInd).svd.tp,[2 3 1]);
    %     %     %%% time-resolved coherence spectrum
    %     %     paramSvdGram = param;
    %     %     % [paramSvdGram.tapers(1),~,paramSvdGram.tapers(2)] = K2W(funTs(sInd).nframes.*tr,ceil(param.tapers(2)*2),verbose);
    %     %     % [~,funPsd(sInd).svdGram.f,funPsd(sInd).svdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramSvdGram.win(1)/tr,1,:,1),paramSvdGram,funTs(sInd).t(1,1,1,1:paramSvdGram.win(1)/tr,1,1),tp);
    %     %     [~,funPsd(sInd).svdGram.f,funPsd(sInd).svdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramSvdGram.win(1)/tr,1,:,1),paramSvdGram,funTs(sInd).t(1:paramSvdGram.win(1)/tr,1,1),tp);
    %     %     funPsd(sInd).svdGram.tp = permute(funPsd(sInd).svdGram.tp,[2 3 1]);
    %     % end
    %     % 
    %     % % paramSvdGram = paramPsdGram;
    %     % % paramSvdGram.tapers(2) = paramSvdGram.tapers(2)*2;
    %     % % paramSvdGram.tapers(1) = (paramSvdGram.tapers(2)+1)/2;
    %     % % % use time windows
    %     % % [~,funPsd(sInd).svdGram.f,funPsd(sInd).svdGram.tp] = mtspectrumc4(funTs(sInd).vec(1:paramSvdGram.win(1)/tr,1,:,1),paramSvdGram,funTs(sInd).t(1,1,1,1:paramSvdGram.win(1)/tr,1,1),tp);
    %     % % funPsd(sInd).svdGram.tp = permute(funPsd(sInd).svdGram.tp,[2 3 1]);
    %     % % if paramSvdGram.win(2)~=inf
    %     % %     % prepare for also performng svd over the full timeseries,
    %     % %     % using half as many tapers
    %     % %     paramSvd = paramSvdGram;
    %     % %     [paramSvd.tapers(1),~,paramSvd.tapers(2)] = K2W(funTs(sInd).nframes.*tr,ceil(paramSvdGram.tapers(2)/2),verbose);
    %     % %     [~,funPsd(sInd).svd.f,funPsd(sInd).svd.tp] = mtspectrumc4(funTs(sInd).vec(:,1,:,1), paramSvd,funTs(sInd).t(1,1,1,:,1,1),tp);
    %     % %     funPsd(sInd).svd.tp = permute(funPsd(sInd).svd.tp,[2 3 1]);
    %     % % end
    % end
    % 
    % 
    % 
    % if ~isempty(tp); tp = []; end
    % 
    % % szPsd = size(funTs(sInd).vec);%single-window over the full time series
    % % if ~skipSVD
    % %     szSvd = size(funTs(sInd).vec);
    % % end
    % % if ~skipGram
    % %     szPsdGram = size(funTs(sInd).vec);
    % %     if ~skipSVD
    % %         szSvdGram = size(funTs(sInd).vec);
    % %     end
    % % end
    % 
    % 
    % %%% Number of frequency points
    % szPsd(1) = length(funPsd(sInd).psd.f);
    % if ~skipSVD
    %     szSvd(1) = length(funPsd(sInd).svd.f);
    % end
    % if ~skipGram
    %     szPsdGram(1) = length(funPsd(sInd).psdGram.f);
    %     if ~skipSVD
    %         szSvdGram(1) = length(funPsd(sInd).svdGram.f);
    %     end
    % end
    % if ~skipTrialGram
    %     szPsdTrialGram(1) = length(funPsd(sInd).psdTrialGram.f);
    %     if ~skipSVD
    %         szSvdTrialGram(1) = length(funPsd(sInd).svdTrialGram.f);
    %     end
    % end


    % % % %%% Number of tapers to save



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


    % if ~skipSVD
    %     szSvd(3) = 1;
    % end
    % if ~skipGram
    %     if ~skipSVD
    %         szSvdGram(3) = 1;
    %     end
    %     szPsdGram(3) = 1;
    % end
    % if ~skipTrialGram
    %     szPsdTrialGram(3) = 1;
    %     if ~skipSVD
    %         szSvdTrialGram(3) = 1;
    %     end
    % end
    % 
    % 
    % %%% Number of runs in the same structure
    % if ~skipPSD
    %     szPsd(4) =1;
    % end
    % if ~skipSVD
    %     szSvd(4) =1;
    % end
    % if ~skipGram
    %     if ~skipPSD
    %         szPsdGram(4) =1;
    %     end
    %     if ~skipSVD
    %         szSvdGram(4) =1;
    %     end
    % end
    % if ~skipTrialGram
    %     if ~skipPSD
    %         szPsdTrialGram(4) =1;
    %     end
    %     if ~skipSVD
    %         szSvdTrialGram(4) =1;
    %     end
    % end
    % 
    % % szPsd(4) = size(funTs(sInd).vec,4);
    % % funPsd(sInd).psd.vec = nan(szPsd,class(funTs(sInd).vec));
    % % if ~skipSVD
    % %     szSvd = szPsd; szSvd(1) = length(funPsd(sInd).psd.f);
    % % end
    % 
    % 
    % %%% Number of time windows
    % if ~skipGram
    %     if ~skipPSD
    %         szPsdGram(5) = size(allWin,1); %number of windows
    %     end
    %     if ~skipSVD
    %         szSvdGram(5) = size(allWin,1); %number of windows
    %     end
    % end
    % if ~skipTrialGram
    %     if ~skipPSD
    %         szPsdTrialGram(5) = size(allWinTrialLock,1); %number of windows
    %     end
    %     if ~skipSVD
    %         szSvdTrialGram(5) = size(allWinTrialLock,1); %number of windows
    %     end
    % end



    % %%% allocate
    % tmp = ones(1,6); tmp(1:length(szPsd)) = szPsd; szPsd = tmp;
    % tmp = ones(1,6); tmp(1:length(szSvd)) = szSvd; szSvd = tmp;
    % tmp = ones(1,6); tmp(1:length(szPsdGram)) = szPsdGram; szPsdGram = tmp;
    % tmp = ones(1,6); tmp(1:length(szSvdGram)) = szSvdGram; szSvdGram = tmp;
    % tmp = ones(1,6); tmp(1:length(szPsdTrialGram)) = szPsdTrialGram; szPsdTrialGram = tmp;
    % tmp = ones(1,6); tmp(1:length(szSvdTrialGram)) = szSvdTrialGram; szSvdTrialGram = tmp;
    % 
    % szPsd = szPsd([1 3 2 4 5 6]);
    % szPsdGram = szPsdGram([1 3 2 4 5 6]);
    % szPsdTrialGram = szPsdTrialGram([1 3 2 4 5 6]);
    % 
    % szSvd = szSvd([1 3 2 4 5 6]);
    % szPsdGram = szPsdGram([1 3 2 4 5 6]);
    % szPsdTrialGram = szPsdTrialGram([1 3 2 4 5 6]);
    % 
    % sz = [length(funPsd(sInd).psd.f) paramPsd.tapers(2) size(funTs(sInd).vec,2) size(funTs(sInd).vec,4)];
    % if avTapers; sz(2) = 1; end
    % funPsd(sInd).psd.vec = nan(sz,class(funTs(sInd).vec)); % freq taper vox run
    % 
    % 
    % %PSD: 'freq/time' 'tapers' 'vox'       'run' 'timeWindow' 'mode'
    % %SVD: 'vox'       'mode'   'freq/time' 'run' 'timeWindow' 'tapers'
    % funPsd(sInd).psdGram.vec = nan(szPsdGram([]),class(funTs(sInd).vec)); % freq taper vox run
    % funPsd(sInd).psdTrialGram.vec = nan(szPsdTrialGram([]),class(funTs(sInd).vec)); % freq taper vox run
    % 
    % 
    % 
    % if ~skipGram
    %     sz = [length(funPsd(sInd).psdGram.f) paramPsdGram.tapers(2) size(funTs(sInd).vec,2) size(funTs(sInd).vec,4)];
    %     sz(2) = 1; sz(5) = size(allWin,1);
    %     funPsd(sInd).psdGram.vec = nan(sz,class(funTs(sInd).vec)); % freq taper vox run
    % end

    %% Compute (also with taper-level permutation)
    % TP.psd.tp = funPsd.psd.tp;
    % if isfield(funPsd,'psdGram')
    %     TP.psdGram.tp = funPsd.psdGram.tp;
    % else
    %     paramPsdGram = [];
    % end
    % if ~skipSVD
    %     TP.svd.tp = funPsd.svd.tp;
    % end
    % if isfield(funPsd,'svdGram') && ~skipSVD
    %     TP.svdGram.tp = funPsd.svdGram.tp;
    % else
    %     paramSvdGram = [];
    % end
    % if ~skipGram
    %     TP.psdTrialGram.tp = funPsd.psdTrialGram.tp;
    %     TP.svdTrialGram.tp = funPsd.svdTrialGram.tp;
    % end

    skip.psd       = skipPSD;
    skip.svd       = skipSVD;
    skip.gram      = skipGram;
    skip.trialGram = skipTrialGram;
    paramInit = param;
    param.psd          = rmfield(paramInit,'onsetList');
    param.svd          = rmfield(paramInit,'onsetList');
    param.psdGram      = rmfield(paramInit,'onsetList');
    param.svdGram      = rmfield(paramInit,'onsetList');
    param.psdTrialGram = paramInit;
    param.svdTrialGram = paramInit;
    clear paramInit
    param.psdGram.winInd      = allWin;
    param.svdGram.winInd      = allWin;
    param.psdTrialGram.winInd = allWinTrialLock;
    param.svdTrialGram.winInd = allWinTrialLock;
    
    
    tic
    [funPsd.psd,funPsd.psdGram,funPsd.svd,funPsd.svdGram,funPsd.psdTrialGram,funPsd.svdTrialGram] = computeAll(funTs,TP,[],skip,nBloc,avTapers,param,verbose,tr,taperPerm,[],cohFperm);
    disp('+++++')
    toc
    disp('+++++')
    %% Compute with fourrier domain phase randomization
    tic
    if phaseRand
        dbstack; error('double check that')
        % cohF = [0 0.2];
        nf = nnz(funPsd.svd.f>=cohFperm(1) & funPsd.svd.f<=cohFperm(2));
        funTsShuf = funTs;
        % sz = size(funPsd.svd.u,1:5); sz(3) = nf;
        % funPsd.svd.uPhaseRand = zeros([size(funPsd.svd.u,1:5) phaseRand],'single');
        sz = size(funPsd.svd.s,1:5); sz(3) = nf;
        funPsd.svd.sPhaseRand = zeros([sz phaseRand],'single');
        % sz = size(funPsd.svd.v,1:5); sz(3) = nf;
        % funPsd.svd.vPhaseRand = zeros([sz phaseRand],'single');
        sz = size(funPsd.svd.coh,1:5); sz(3) = nf;
        funPsd.svd.cohPhaseRand = zeros([sz phaseRand],'single');
        % sz = size(funPsd.svdGram.u,1:5); sz(3) = nf;
        % funPsd.svdGram.uPhaseRand = zeros([sz phaseRand],'single');
        % sz = size(funPsd.svdGram.s,1:5); sz(3) = nf;
        % funPsd.svdGram.sPhaseRand = zeros([sz phaseRand],'single');
        % sz = size(funPsd.svdGram.v,1:5); sz(3) = nf;
        % funPsd.svdGram.vPhaseRand = zeros([sz phaseRand],'single');
        % sz = size(funPsd.svdGram.coh,1:5); sz(3) = nf;
        % funPsd.svdGram.cohPhaseRand = zeros([sz phaseRand],'single');
        for phaseRandInd = 1:phaseRand
            disp(['phase randomization ' num2str(phaseRandInd) '/' num2str(phaseRand)])
            [X,Y] = pol2cart(rand(size(funTs(sInd).vec)).*(2*pi),abs(fft(funTs(sInd).vec,[],1)));
            funTsShuf(sInd).vec = ifft(complex(X,Y));
            paramTmp = param; paramTmp.win = [inf inf];
            [~,~,svd,svdGram] = computeAll(funTsShuf,TP,szPsd,1,skipSVD,nBloc,avTapers,paramPsd,paramSvd,0,tr,paramTmp,allWin,paramPsdGram,paramSvdGram,[],cohFperm,[]);
            % funPsdShuf = computeAll(funTsShuf,funPsd,szPsd,1,skipSVD,nBloc,avTapers,paramPsd,paramSvd,0,sInd,tr,param,allWin,paramPsdGram,paramSvdGram);
            % funPsd.svd.uPhaseRand(:,:,:,:,:,phaseRandInd) = svd.u;
            funPsd.svd.sPhaseRand(:,:,:,:,:,phaseRandInd) = svd.s;
            % funPsd.svd.vPhaseRand(:,:,:,:,:,phaseRandInd) = svd.v;
            funPsd.svd.cohPhaseRand(:,:,:,:,:,phaseRandInd) = svd.coh;
            % funPsd.svdGram.uPhaseRand(:,:,:,:,:,phaseRandInd) = svdGram.u;
            % funPsd.svdGram.sPhaseRand(:,:,:,:,:,phaseRandInd) = svdGram.s;
            % funPsd.svdGram.vPhaseRand(:,:,:,:,:,phaseRandInd) = svdGram.v;
            % funPsd.svdGram.cohPhaseRand(:,:,:,:,:,phaseRandInd) = svdGram.coh;
        end
    end
    disp('+++++')
    toc
    disp('+++++')

    

    %% Sort outputs
    funPsd(sInd).psd.f = permute(funPsd(sInd).psd.f,[2 1 3 4 5 6]);
    funPsd(sInd).psd.tp = TP.full;
    funPsd(sInd).psd.vec;
    funPsd(sInd).psd.K = param.psd.tapers(2);
    [~,funPsd(sInd).psd.W,~] = K2W(funPsd(sInd).psd.T,funPsd(sInd).psd.K,0);
    funPsd(sInd).psd.param = param.psd;
    funPsd(sInd).psd.info = strjoin({'freq/time' 'tapers' 'vox' 'run' 'timeWindow' 'mode'},' x ');
    funPsd(sInd).psd.mask = funTs(sInd).vol2vec;

    if ~skipGram
        funPsd(sInd).psdGram.f = permute(funPsd(sInd).psdGram.f,[2 1 3 4 5 6]);
        funPsd(sInd).psdGram.tp = TP.gram;
        funPsd(sInd).psdGram.vec;
        % funPsd(sInd).psdGram.t;
        % funPsd(sInd).psdGram.tWin = mean(funPsd(sInd).psdGram.t,1);
        funPsd(sInd).psdGram.lWin = param.psdGram.win(1);
        funPsd(sInd).psdGram.K = param.psdGram.tapers(2);
        [~,funPsd(sInd).psdGram.W,~] = K2W(funPsd(sInd).psdGram.T,funPsd(sInd).psdGram.K,0);
        funPsd(sInd).psdGram.param = param.psdGram;
        funPsd(sInd).psdGram.info = strjoin({'freq/time' 'tapers' 'vox' 'run' 'timeWindow' 'mode'},' x ');
        funPsd(sInd).psdGram.mask = funTs(sInd).vol2vec;
    end
    
    if ~skipSVD
        funPsd(sInd).svd.f;
        funPsd(sInd).svd.tp = permute(TP.full,[3 4 1 5 6 2]);
        funPsd(sInd).svd.u;
        funPsd(sInd).svd.s = permute(funPsd(sInd).svd.s,[2 1 3 4 5 6]);
        funPsd(sInd).svd.v = permute(funPsd(sInd).svd.v,[6 2 3 4 5 1]);
        funPsd(sInd).svd.coh = permute(funPsd(sInd).svd.coh,[2 1 3 4 5 6]);
        funPsd(sInd).svd.K = param.svd.tapers(2);
        [~,funPsd(sInd).svd.W,~] = K2W(funPsd(sInd).svd.T,funPsd(sInd).svd.K,0);
        funPsd(sInd).svd.param = param.svd;
        funPsd(sInd).svd.info = strjoin({'vox' 'mode' 'freq/time' 'run' 'timeWindow' 'tapers'},' x ');
        funPsd(sInd).svd.mask = funTs(sInd).vol2vec;
    else
        funPsd(sInd).svd = [];
    end

    if ~skipGram && ~skipSVD
        funPsd(sInd).svdGram.f;
        funPsd(sInd).svdGram.tp = permute(TP.gram,[3 4 1 5 6 2]);
        funPsd(sInd).svdGram.u;
        funPsd(sInd).svdGram.s = permute(funPsd(sInd).svdGram.s,[2 1 3 4 5 6]);
        funPsd(sInd).svdGram.v = permute(funPsd(sInd).svdGram.v,[6 2 3 4 5 1]);
        funPsd(sInd).svdGram.coh = permute(funPsd(sInd).svdGram.coh,[2 1 3 4 5 6]);
        % funPsd(sInd).svdGram.t;
        % funPsd(sInd).svdGram.tWin = mean(funPsd(sInd).svdGram.t,1);
        funPsd(sInd).svdGram.lWin = param.svdGram.win(1);
        funPsd(sInd).svdGram.K = param.svdGram.tapers(2);
        [~,funPsd(sInd).svdGram.W,~] = K2W(funPsd(sInd).svdGram.T,funPsd(sInd).svdGram.K,0);
        funPsd(sInd).svdGram.param = param.svdGram;
        funPsd(sInd).svdGram.info = strjoin({'vox' 'mode' 'freq/time' 'run' 'timeWindow' 'tapers'},' x ');
        funPsd(sInd).svdGram.mask = funTs(sInd).vol2vec;
    else
        funPsd(sInd).svdGram = [];
    end


    if ~skip.trialGram
        if ~skip.psd
            funPsd(sInd).psdTrialGram.f = permute(funPsd(sInd).psdTrialGram.f,[2 1 3 4 5 6]);
            funPsd(sInd).psdTrialGram.tp = TP.trialGram;
            funPsd(sInd).psdTrialGram.vec;
            funPsd(sInd).psdTrialGram.lWin = param.psdTrialGram.win(1);
            funPsd(sInd).psdTrialGram.K = param.psdTrialGram.tapers(2);
            [~,funPsd(sInd).psdTrialGram.W,~] = K2W(funPsd(sInd).psdTrialGram.T,funPsd(sInd).psdTrialGram.K,0);
            funPsd(sInd).psdTrialGram.param = param.psdTrialGram;
            funPsd(sInd).psdTrialGram.info = strjoin({'freq/time' 'tapers' 'vox' 'run' 'timeWindow' 'mode'},' x ');
            funPsd(sInd).psdTrialGram.mask = funTs(sInd).vol2vec;
        else
            funPsd(sInd).psdTrialGram = [];
        end
        if ~skip.svd
            funPsd(sInd).svdTrialGram.f;
            funPsd(sInd).svdTrialGram.tp = permute(TP.trialGram,[3 4 1 5 6 2]);
            % funPsd(sInd).svdTrialGram.u;
            % funPsd(sInd).svdTrialGram.s = permute(funPsd(sInd).svdTrialGram.s,[2 1 3 4 5 6]);
            % funPsd(sInd).svdTrialGram.v = permute(funPsd(sInd).svdTrialGram.v,[6 2 3 4 5 1]);
            funPsd(sInd).svdTrialGram.coh = permute(funPsd(sInd).svdTrialGram.coh,[2 1 3 4 5 6]);
            funPsd(sInd).svdTrialGram.lWin = param.svdTrialGram.win(1);
            funPsd(sInd).svdTrialGram.K = param.svdTrialGram.tapers(2);
            [~,funPsd(sInd).svdTrialGram.W,~] = K2W(funPsd(sInd).svdTrialGram.T,funPsd(sInd).svdTrialGram.K,0);
            funPsd(sInd).svdTrialGram.param = param.svdTrialGram;
            funPsd(sInd).svdTrialGram.info = strjoin({'vox' 'mode' 'freq/time' 'run' 'timeWindow' 'tapers'},' x ');
            funPsd(sInd).svdTrialGram.mask = funTs(sInd).vol2vec;
        else
            funPsd(sInd).svdTrialGram = [];
        end
    end
    
    %% For backward compatibility
    funPsd(sInd).vec = permute(funPsd(sInd).psd.vec,[1 3 2 4]); funPsd(sInd).psd.vec = [];
    funPsd(sInd).f = permute(funPsd(sInd).psd.f,[2 3 4 5 6 1]);
    funPsd(sInd).K = funPsd(sInd).psd.K;
    funPsd(sInd).T = funPsd(sInd).psd.T;
    funPsd(sInd).W = funPsd(sInd).psd.W;
    funPsd(sInd).tr = mode(diff(funPsd(sInd).f))*1000;
    funPsd(sInd).nfreq = size(funPsd(sInd).vec,1);
    funPsd(sInd).nframes=size(funTs(sInd).vec,1);
    funPsd(sInd).nruns=size(funPsd(sInd).vec,4);
    % if funPsd(sInd).nframes~=size(funTs(sInd).vec,1); dbstack; warning('something wrong with nframes'); end
    % if funPsd(sInd).nruns~=size(funPsd(sInd).vec,4); dbstack; warning('something wrong with nframes'); end
    funPsd(sInd).ntapers = size(funPsd(sInd).vec,3);
end


function [psd,psdGram,svd,svdGram,psdTrialGram,svdTrialGram] = computeAll(funTs,TP,nRun,skip,nBloc,avTapers,param,verbose,tr,nShuf,cohFrange,cohFperm)
if ~exist('nShuf','var');         nShuf = []; end
if ~exist('nRun','var');           nRun = []; end
if ~exist('cohFrange','var'); cohFrange = []; end
if isempty(nShuf);         nShuf = 0; end
if isempty(nRun);           nRun = 1; end
if isempty(cohFrange); cohFrange = [0 inf]; end
orderWin = 0;
% skip.gram = any(ismember(param.win(2),[0 inf nan]));
for runInd = 1:nRun
    if verbose && nRun>1; disp(['---Run ' num2str(runInd) '/' num2str(nRun) '---']); end
    
    
    %% Compute power spectra
    if ~skip.psd
        
        %%% psd over the full time series
        if verbose; disp('full power spectra: computing'); end
        tp = permute(TP.full,[1 3 2]);


% %%%%%%%%%%%
% avTapers = 1;
% 
% vec = funTs.vec(:,:,:,runInd);
% vec = vec - mean(vec,1);
% 
% [S,f]=mtspectrumc(vec,param);
% 
% curVerbose = 1;
% [vecFFT,fFFT] = fastMtPSD(tp,vec,param.psd.Fs,[],nBloc,avTapers,curVerbose);
% 
% NT = size(vec,1);
% Fs = param.psd.Fs;
% tvec=1/Fs *(0:NT-1)';
% [vecFSeries,fFSeries] = fastMtPSD(tp,vec,param.psd.Fs,[],nBloc,avTapers,curVerbose,tvec);
% 
% figure('WindowStyle','docked');
% k = 1; v = 1;
% h1 = plot(fFFT,vecFFT(:,k,v)); hold on
% h2 = plot(fFSeries,vecFSeries(:,k,v));
% h3 = plot(f,S(:,v,k));
% xline(1/mean(diff(param.psdTrialGram.onsetList)))
% legend([h1 h2 h3],{'fft' 'fSeries' 'Chronux'})
% %%%%%%%%%%%




        [psd.vec(:,:,:,runInd),psd.f] = fastMtPSD(tp,funTs.vec(:,:,:,runInd),param.psd.Fs,[],nBloc,avTapers);
        psd.T(:,:,:,runInd) = size(funTs.vec(:,:,:,runInd),1).*tr;
        if verbose; disp('full power spectra: done'); end
        
        %%% psd over each time window
        if ~skip.gram
            tp = permute(TP.gram,[1 3 2]);
            [N,E] = size(param.psdGram.winInd,[2 3]);
            K = param.psdGram.tapers(2);
            C = size(funTs.vec,2);
            for winInd = 1:size(param.psdGram.winInd,1)
                if verbose; disp(['time-resolved power spectra: computing window ' num2str(winInd) '/' num2str(size(param.psdGram.winInd,1))]); end
                if runInd==1 && winInd==1; curVerbose = 2; else curVerbose = 1; end

                %%%% detrend each window
                if E>1 || nRun>1; dbstack; error('double-check that'); end
                vec = funTs.vec(param.psdGram.winInd(winInd,:),:,:,runInd); % [time x vox x 1 x run] %[N*E C 1 nRun]
                if orderWin~=-1
                    if orderWin==0
                        vec = vec - mean(vec,1); %[N nRun*C*1*E]
                    else
                        vec = dtrnd2(vec,tr,[],orderWin); %[N nRun*C*1*E]
                    end
                end

                %%%% compute
                [psdGram.vec(:,:,:,runInd,winInd),psdGram.f] = fastMtPSD(tp,vec,param.psdGram.Fs,[],nBloc,avTapers,curVerbose);
                psdGram.T(:,:,:,runInd,winInd) = size(funTs.vec(param.psdGram.winInd(winInd,:),:,:,runInd),1).*tr;
                % psdGram.tWin(1,1,1,runInd,winInd,1) = mean(funTs.t(1,1,1,param.psdGram.winInd(winInd,[1 end]),1,runInd));
                psdGram.tWin(1,1,1,runInd,winInd,1) = funTs.t(param.psdGram.winInd(winInd,1),1,runInd) + param.psdGram.win(1)/2;
                % psdGram.tWin(1,1,1,runInd,winInd,1) = mean(funTs.t(param.psdGram.winInd(winInd,[1 end]),1,runInd));
            end
            if verbose; disp('time-resolved power spectra: done'); end
        else
            psdGram = [];
        end


        %%% psd over each multi-trial time window
        if ~skip.trialGram
            tp = permute(TP.trialGram,[1 3 2]);
            [N,E] = size(param.psdTrialGram.winInd,[2 3]);
            K = param.psdTrialGram.tapers(2);
            C = size(funTs.vec,2);
            nRun;
            for winInd = 1:size(param.psdTrialGram.winInd,1)
                if verbose; disp(['trial-locked power spectra: computing window ' num2str(winInd) '/' num2str(size(param.psdTrialGram.winInd,1))]); end
                if runInd==1 && winInd==1; curVerbose = 2; else curVerbose = 1; end
                
                %%%% detrend on a trial-by-trial basis
                vec = funTs.vec(param.psdTrialGram.winInd(winInd,:,:),:,:,runInd); % [time x vox x 1 x run] %[N*E C 1 nRun]
                vec = permute(vec,[4 2 3 1]); %[nRun C 1 N*E]
                vec = reshape(vec,[nRun C 1 N E]); %[nRun C 1 N E]
                if orderWin~=-1
                    vec = permute(vec,[4 1 2 3 5]); %[N nRun C 1 E]
                    vec = reshape(vec,[N nRun*C*1*E]); %[N nRun*C*1*E]
                    if orderWin==0
                        vec = vec - mean(vec,1); %[N nRun*C*1*E]
                    else
                        vec = dtrnd2(vec,tr,[],orderWin); %[N nRun*C*1*E]
                    end
                    vec = reshape(vec,[N nRun C 1 E]); %[N nRun C 1 E]
                    vec = permute(vec,[2 3 4 1 5]); %[nRun C 1 N E]
                end
                vec = permute(vec,[4 2 5 1 3]); %[N C E nRun 1]


                % vec = funTs.vec(param.psdTrialGram.winInd(winInd,:,:),:,:,runInd); % [time x vox x 1 x run] %[N*E C 1 nRun]
                % if orderWin~=-1
                %     vec = permute(vec,[4 2 3 1]); %[nRun C 1 N*E]
                %     vec = reshape(vec,[nRun C 1 N E]); %[nRun C 1 N E]
                %     vec = permute(vec,[4 1 2 3 5]); %[N nRun C 1 E]
                %     vec = reshape(vec,[N nRun*C*1*E]); %[N nRun*C*1*E]
                %     if orderWin==0
                %         vec = vec - mean(vec,1); %[N nRun*C*1*E]
                %     else
                %         vec = dtrnd2(vec,tr,[],orderWin); %[N nRun*C*1*E]
                %     end
                %     vec = reshape(vec,[N nRun C 1 E]); %[N nRun C 1 E]
                %     vec = permute(vec,[2 3 4 1 5]); %[nRun C 1 N E]
                %     vec = reshape(vec,[nRun C 1 N*E]); %[nRun C 1 N*E]
                %     vec = permute(vec,[4 2 3 1]); %[N*E C 1 nRun]
                % end
                
                %%%% precompute tvec to allow phase-reset at stimOnset
                tvec = zeros(N*E,E);
                for trialInd = 1:E
                    tvec(:,trialInd) = funTs.t(param.psdTrialGram.winInd(winInd,:),:,:,runInd) - param.onsetList(trialInd,1); % subtract stim onset from the tvec of each trial--this will effectively reset all phases to 0
                end
                tvec = repmat(tvec,[1 1 K]);
                % tvec = reshape(permute(repmat(tvec,[1 1 K]),[1 3 2]),[N*E,K*E]);

                %%%% compute psd
                curAvTaper = 2; % 0: don't average   1: average   2: phase-coherent cross-trial average
                [psdTrialGram.vec(:,:,:,runInd,winInd),psdTrialGram.f] = fastMtPSD(tp,vec,param.psdTrialGram.Fs,[],nBloc,curAvTaper,curVerbose,tvec);
                psdTrialGram.T(:,:,:,runInd,winInd) = size(funTs.vec(param.psdTrialGram.winInd(winInd,:),:,:,runInd),1).*tr;
                psdTrialGram.tWin(1,1,1,runInd,winInd,1) = funTs.t(param.psdTrialGram.winInd(winInd,1)) + param.psdTrialGram.win(1)/2;

                %%%% compute svd
                if ~skip.svd
                    if ~exist('cohF','var') || isempty(cohFrange); cohFrange = [0 inf]; end
                    nMode = inf;
                    Fs = param.Fs;
                    Fpass = param; if isfield(Fpass,'fpass'); Fpass = Fpass.fpass; else Fpass = [0 Fs/2]; end
                    NFFT=max(2^(nextpow2(N*E)+0),N*E);
                    f = permute(getfgrid(Fs,NFFT,Fpass),[1 3 2]);
                    F = length(f);
                    J = getJ(vec,tp,tvec,f)/Fs; %[freq x taper x vox x trial]
                    switch curAvTaper
                        case 0
                        case 1
                            J = permute(J,[1 3 2 4]); %[freq x vox         x taper       x trial]
                            J = reshape(J,[F C K*E]); %[freq x vox         x taper*trial        ]
                            J = permute(J,[1 3 2 4]); %[freq x taper*trial x vox                ]
                        case 2
                            J = mean(J,4);
                    end
                    [u,s,v] = pagesvd(  permute(J,[3 2 1 4])  ,"econ","vector");
                    
                    c = s.^2./sum(s.^2,1); % coherence
                    
                    svdTrialGram.coh(:,:,:,runInd,winInd) = c; %[mode x 1 x freq]
                    svdTrialGram.T(:,:,:,runInd,winInd) = size(funTs.vec(param.svdTrialGram.winInd(winInd,:),:,:,runInd),1).*tr;
                    svdTrialGram.tWin(1,1,1,runInd,winInd,1) = funTs.t(param.svdTrialGram.winInd(winInd,1),:,:,runInd) + param.svdTrialGram.win(1)/2;
                end

                % sz   = size(param.psdTrialGram.winInd(winInd,:,:)); sz([1 4]) = size(funTs.vec,[2 4]);
                % vec  = funTs.vec(param.psdTrialGram.winInd(winInd,:),:,:,runInd); % [time x vox x trial x run]
                % vec  = vec - mean(vec,1);
                % 
                % %%%% precompute tvec to allow phase-reset at stimOnset
                % tvec = zeros(size(param.psdTrialGram.winInd(winInd,:,:)));
                % for trialInd = 1:size(param.psdTrialGram.winInd,3)
                %     tvec(:,:,trialInd) = funTs.t(param.psdTrialGram.winInd(winInd,:,trialInd),:,:,runInd) - param.onsetList(trialInd,1); % subtract stim onset from the tvec of each trial--this will effectively reset all phases to 0
                % end
                % tvec = reshape(tvec,size(tvec,1),prod(size(tvec,[2 3])))';
                % 
                % %%%% compute
                % [psdTrialGram.vec(:,:,:,runInd,winInd),psdTrialGram.f] = fastMtPSD(tp,vec,param.psdTrialGram.Fs,[],nBloc,avTapers,curVerbose,tvec);
                % psdTrialGram.T(:,:,:,runInd,winInd) = size(funTs.vec(param.psdTrialGram.winInd(winInd,:),:,:,runInd),1).*tr;
                % psdTrialGram.tWin(1,1,1,runInd,winInd,1) = funTs.t(param.psdTrialGram.winInd(winInd,1)) + param.psdTrialGram.win(1)/2;
                

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % sz   = size(param.psdTrialGram.winInd(winInd,:,:)); sz([1 4]) = size(funTs.vec,[2 4]);
% % vec  = reshape(permute( funTs.vec(param.psdTrialGram.winInd(winInd,:,:),:,:,runInd) ,[2 3 1]),sz); % [time x vox x trial x run]
% % vec  = vec - mean(vec,1);
% % paramX = param.psdTrialGram;
% % paramX.trialave = 1;
% % paramX.err = [1 0.05];
% % 
% % K = paramX.tapers(2);
% % T = size(vec,1);
% % [TW,W,K] = K2W(T,K,verbose);
% % paramX.tapers = [TW K];
% % paramX = rmfield(paramX,'win');
% % 
% % [S,f,Serr]=mtspectrumc(permute(vec(:,1,:),[1 3 2]),paramX)
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                

                

                % % % % %%%% Precompute tvec to allow phase-reset at stimOnset
                % % % % tvec = zeros(size(param.psdTrialGram.winInd(winInd,:,:)));
                % % % % for trialInd = 1:size(param.psdTrialGram.winInd,3)
                % % % %     % param.svdTrialGram.winInd(winInd,:,trialInd)
                % % % %     % tvec(:,:,trialInd) = funTs.t(param.svdTrialGram.winInd(winInd,:,trialInd),:,:,runInd);
                % % % %     tvec(:,:,trialInd) = funTs.t(param.psdTrialGram.winInd(winInd,:,trialInd),:,:,runInd) - param.onsetList(trialInd,1); % subtract stim onset from the tvec of each trial--this will effectively reset all phases to 0
                % % % % end
                % % % % tvec = reshape(tvec,size(tvec,1),prod(size(tvec,[2 3])))';
                % % % 
                % % % [psdTrialGram.vec(:,:,:,runInd,winInd),psdTrialGram.f] = fastMtPSD(tp,vec,param.psdTrialGram.Fs,[],nBloc,avTapers,curVerbose,tvec);
                % % % 
                % % % % [psdTrialGram.vec(:,:,:,runInd,winInd),psdTrialGram.f] = fastMtPSD(tp,vec,param.psdTrialGram.Fs,[],nBloc,avTapers,curVerbose);
                % % % psdTrialGram.T(:,:,:,runInd,winInd) = size(funTs.vec(param.psdTrialGram.winInd(winInd,:),:,:,runInd),1).*tr;
                % % % psdTrialGram.tWin(1,1,1,runInd,winInd,1) = funTs.t(param.psdTrialGram.winInd(winInd,1)) + param.psdTrialGram.win(1)/2;
            end
            if skip.svd
                svdTrialGram = [];
            else
                svdTrialGram.f = f;
                svdTrialGram.info = 'taper/vox/mode x mode x freq x run x window';
            end
            if verbose; disp('time-resolved power spectra: done'); end
        else
            psdTrialGram = [];
        end

    else
        psd = [];
        psdGram = [];
        psdTrialGram = [];
    end

    
    %% Compute coherence spectra
    if ~skip.svd
        if ~exist('cohF','var') || isempty(cohFrange); cohFrange = [0 inf]; end
        nMode = inf;

        %%% svd over the full time series
        if verbose; disp('full coherence spectra: computing'); end
        tp = TP.full;
        [svd.u(:,:,:,runInd),...
            svd.s(:,:,:,runInd),...
            svd.v(:,:,:,runInd),...
            svd.coh(:,:,:,runInd),...
            svd.f]...
            = fastKleinMtSVD(tp,funTs.vec(:,:,:,runInd),param.svd.Fs,funTs.t(:,:,:,runInd),cohFrange,nMode);
        % = fastKleinMtSVD(svd.tp,funTs.vec(:,:,:,runInd),paramSvd.Fs,funTs.t(:,:,:,:,:,runInd),cohF,nMode);
        svd.T(:,:,:,runInd) = size(funTs.t,1).*tr;
        svd.info = 'taper/vox/mode x mode x freq';
        if verbose; disp('full coherence spectra: done'); end
        %%%% Random permutation in taper space
        if nShuf
            dbstack; error('double-check that')
            if nShuf<100; warning('nShuf<100, width of the confidence interval will be underestimated'); end
            if ~exist('cohFperm','var') || isempty(cohFperm); cohFperm = [0 0.2]; end
            disp('full coherence spectra with taper permutation: computing')
            [svd.uNullMean,...
                svd.sNullMean,...
                svd.vNullMean,...
                svd.cohNullMean,...
                svd.fNull,...
                svd.uNull90,...
                svd.sNull90,...
                svd.vNull90,...
                svd.cohNull90]...
                = fastKleinMtSVD(tp,funTs.vec(:,:,:,runInd),param.svd.Fs,funTs.t(:,:,:,runInd),cohFperm,nMode,nShuf);
            svd.infoTaperPerm = 'taper/vox/mode x mode x freq x perm';
            disp('full coherence spectra with taper permutation: done')
        end


        %%% svd over each time window
        if ~skip.gram
            tp = TP.gram;
            for winInd = 1:size(param.svdGram.winInd,1)
                if verbose; disp(['time-resolved coherence spectra: computing window ' num2str(winInd) '/' num2str(size(param.svdGram.winInd,1))]); end

                %%%% detrend each window
                vec = funTs.vec(param.svdGram.winInd(winInd,:),:,:,runInd); % [time x vox x 1 x run] %[N*E C 1 nRun]
                if orderWin~=-1
                    if orderWin==0
                        vec = vec - mean(vec,1); %[N nRun*C*1*E]
                    else
                        vec = dtrnd2(vec,tr,[],orderWin); %[N nRun*C*1*E]
                    end
                end

                %%%% compute
                [svdGram.u(:,:,:,runInd,winInd),...
                    svdGram.s(:,:,:,runInd,winInd),...
                    svdGram.v(:,:,:,runInd,winInd),...
                    svdGram.coh(:,:,:,runInd,winInd),...
                    svdGram.f]...
                    = fastKleinMtSVD(tp,vec,param.svdGram.Fs,funTs.t(param.svdGram.winInd(winInd,:),:,:,runInd),cohFrange,nMode);
                % = fastKleinMtSVD(svdGram.tp,funTs.vec(param.svdGram.winInd(winInd,:),:,:,runInd)-mean(funTs.vec(param.svdGram.winInd(winInd,:),:,:,runInd),1),paramSvdGram.Fs,funTs.t(:,:,:,param.svdGram.winInd(winInd,:),:,runInd),cohF,nMode);
                svdGram.T(:,:,:,runInd,winInd) = size(funTs.vec(param.svdGram.winInd(winInd,:),:,:,runInd),1).*tr;
                svdGram.tWin(1,1,1,runInd,winInd,1) = funTs.t(param.svdGram.winInd(winInd,1),:,:,runInd) + param.svdGram.win(1)/2;
                svdGram.info = 'taper/vox/mode x mode x freq x run x window';

                if nShuf
                    dbstack; error('code that')
                    disp(['time-resolved coherence spectra: computing window ' num2str(winInd) '/' num2str(size(param.svdGram.winInd,1)) '; taper permutation: computing'])
                    [svdGram.uTaperPerm(:,:,:,:,runInd,winInd),...
                        svdGram.sTaperPerm(:,:,:,:,runInd,winInd),...
                        svdGram.vTaperPerm(:,:,:,:,runInd,winInd),...
                        svdGram.cohTaperPerm(:,:,:,:,runInd,winInd),...
                        ~]...
                        = fastKleinMtSVD(svdGram.tp,funTs.vec(param.svdGram.winInd(winInd,:),:,:,runInd)-mean(funTs.vec(param.svdGram.winInd(winInd,:),:,:,runInd),1),paramSvdGram.Fs,funTs.t(param.svdGram.winInd(winInd,:),:,:,runInd),cohF,nMode,nShuf,0);
                    svdGram.infoTaperPerm = 'taper/vox/mode x mode x freq x perm x run x window';
                    disp(['time-resolved coherence spectra: computing window ' num2str(winInd) '/' num2str(size(param.svdGram.winInd,1)) '; taper permutation: done'])
                end
            end
            if verbose; disp('time-resolved coherence spectra: done'); end
        else
            svdGram = [];
        end



        %%% svd over each multi-trial time window
        % if ~skip.trialGram
        %     tp = TP.trialGram;
        %     [N,E] = size(param.svdTrialGram.winInd,[2 3]);
        %     K = param.svdTrialGram.tapers(2);
        %     C = size(funTs.vec,2);
        %     nRun;
        %     for winInd = 1:size(param.svdTrialGram.winInd,1)
        %         if verbose; disp(['trial-locked coherence spectra: computing window ' num2str(winInd) '/' num2str(size(param.svdTrialGram.winInd,1))]); end
        % 
        %         %%%% detrend on a trial-by-trial basis
        %         vec = funTs.vec(param.svdTrialGram.winInd(winInd,:,:),:,:,runInd); % [time x vox x 1 x run] %[N*E C 1 nRun]
        %         if orderWin~=-1
        %             vec = permute(vec,[4 2 3 1]); %[nRun C 1 N*E]
        %             vec = reshape(vec,[nRun C 1 N E]); %[nRun C 1 N E]
        %             vec = permute(vec,[4 1 2 3 5]); %[N nRun C 1 E]
        %             vec = reshape(vec,[N nRun*C*1*E]); %[N nRun*C*1*E]
        %             if orderWin==0
        %                 vec = vec - mean(vec,1); %[N nRun*C*1*E]
        %             else
        %                 vec = dtrnd2(vec,tr,[],orderWin); %[N nRun*C*1*E]
        %             end
        %             vec = reshape(vec,[N nRun C 1 E]); %[N nRun C 1 E]
        %             vec = permute(vec,[2 3 4 1 5]); %[nRun C 1 N E]
        %             vec = reshape(vec,[nRun C 1 N*E]); %[nRun C 1 N*E]
        %             vec = permute(vec,[4 2 3 1]); %[N*E C 1 nRun]
        %         end
        % 
        %         %%%%% precompute tvec to allow phase-reset at stimOnset
        %         tvec = zeros(N*E,E);
        %         for trialInd = 1:E
        %             tvec(:,trialInd) = funTs.t(param.svdTrialGram.winInd(winInd,:),:,:,runInd) - param.onsetList(trialInd,1); % subtract stim onset from the tvec of each trial--this will effectively reset all phases to 0
        %         end
        %         tvec = reshape(permute(repmat(tvec,[1 1 K]),[1 3 2]),[N*E,K*E]);
        % 
        %         %%%%% compute
        %         [svdTrialGram.u(:,:,:,runInd,winInd),...
        %             svdTrialGram.s(:,:,:,runInd,winInd),...
        %             svdTrialGram.v(:,:,:,runInd,winInd),...
        %             svdTrialGram.coh(:,:,:,runInd,winInd),...
        %             svdTrialGram.f]...
        %             = fastKleinMtSVD(tp,vec,param.svdTrialGram.Fs,tvec,cohFrange,nMode);
        %         svdTrialGram.T(:,:,:,runInd,winInd) = size(funTs.vec(param.svdTrialGram.winInd(winInd,:),:,:,runInd),1).*tr;
        %         svdTrialGram.tWin(1,1,1,runInd,winInd,1) = funTs.t(param.svdTrialGram.winInd(winInd,1),:,:,runInd) + param.svdTrialGram.win(1)/2;
        % 
        % 
        % 
        %         % %%%% Missing data taper aproach (not working)
        %         % if 0
        %         %     %%%% precompute tvec to allow phase-reset at stimOnset
        %         %     tvec = zeros(size(param.svdTrialGram.winInd(winInd,:,:)));
        %         %     for trialInd = 1:size(param.svdTrialGram.winInd,3)
        %         %         tvec(:,:,trialInd) = funTs.t(param.svdTrialGram.winInd(winInd,:,trialInd),:,:,runInd) - param.onsetList(trialInd,1); % subtract stim onset from the tvec of each trial--this will effectively reset all phases to 0
        %         %     end
        %         %     tvec = reshape(tvec,size(tvec,1),prod(size(tvec,[2 3])))';
        %         % 
        %         %     %%%% compute
        %         %     [svdTrialGram.u(:,:,:,runInd,winInd),...
        %         %         svdTrialGram.s(:,:,:,runInd,winInd),...
        %         %         svdTrialGram.v(:,:,:,runInd,winInd),...
        %         %         svdTrialGram.coh(:,:,:,runInd,winInd),...
        %         %         svdTrialGram.f]...
        %         %         = fastKleinMtSVD(tp,funTs.vec(param.svdTrialGram.winInd(winInd,:),:,:,runInd)-mean(funTs.vec(param.svdTrialGram.winInd(winInd,:),:,:,runInd),1),param.svdTrialGram.Fs,tvec,cohFrange,nMode);
        %         % 
        %         %     svdTrialGram.T(:,:,:,runInd,winInd) = size(funTs.vec(param.svdTrialGram.winInd(winInd,:),:,:,runInd),1).*tr;
        %         %     svdTrialGram.tWin(1,1,1,runInd,winInd,1) = funTs.t(param.svdTrialGram.winInd(winInd,1),:,:,runInd) + param.svdTrialGram.win(1)/2;
        %         %     svdTrialGram.info = 'taper/vox/mode x mode x freq x run x window';
        %         % end
        %     end
        %     svdTrialGram.info = 'taper/vox/mode x mode x freq x run x window';
        %     if verbose; disp('time-trial coherence spectra: done'); end
        % else
        %     svdTrialGram = [];
        % end


    else
        svd = [];
        svdGram = [];
        svdTrialGram = [];
    end
end

