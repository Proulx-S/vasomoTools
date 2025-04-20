function funPsd = runFullMT6(rCond,W,K,winSec,dsgn,mask,skipSVD,skipPSD,force,verbose)
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
if ~exist('extra','var'); extra = []; end
if ~exist('W','var'); W = []; end
if ~exist('K','var'); K = []; end
if isempty(K) && isempty(W); K = 1; end
if ~exist('skipSVD','var') || isempty(skipSVD); skipSVD = false; end
if ~exist('skipPSD','var') || isempty(skipPSD); skipPSD = false; end
if ~exist('mask','var'); mask = []; end
if ~exist('onsets','var'); onsets = []; end
if ~exist('ondurs','var'); ondurs = []; end
if ~exist('taperperm','var'); taperPerm = []; end
if ~exist('phaseRand','var'); phaseRand = []; end
if ~exist('testFlag','var');   testFlag = []; end
if ~isfield(extra,'catFlag'); extra.catFlag = []; end
if ~isfield(extra,'padTo');     extra.padTo = []; end
if isempty(testFlag); testFlag = 0; end
% catFlag = extra.catFlag;
% padTo   = extra.padTo;



% if isempty(onsets)
%     if length(volTs)>1; warning('using first dsgn for all instances of volTs'); end
%     if isfield(volTs(1),'dsgn')
%         if isfield(volTs(1).dsgn,'onsets')
%             onsets = volTs(1).dsgn.onsets;
%         else
%             onsets = volTs(1).dsgn.onsetList;
%         end
%     elseif isfield(volTs(1),'ts') && isfield(volTs(1).ts,'dsgn')
%         onsets = volTs(1).ts.dsgn.onsetList;
%     end
% end
% if size(onsets,1)==1; onsets = onsets'; end
% 
% if isempty(ondurs)
%     if isfield(volTs(1),'dsgn')
%         if isfield(volTs(1).dsgn,'ondurs')
%             ondurs = volTs(1).dsgn.ondurs;
%         else
%             ondurs = volTs(1).dsgn.ondurList;
%         end
%     elseif isfield(volTs(1),'ts') && isfield(volTs(1).ts,'dsgn')
%         ondurs = volTs(1).ts.dsgn.ondurList;
%     end
% end
% if size(ondurs,1)==1; ondurs = ondurs'; end


%% Compute on each run
for I = 1:size(rCond.fPreprocList,1)
    rCond.r = I;
    rCond.R = size(rCond.fPreprocList,1);
    if isempty(rCond.volTs)
        rCond.volTs = vol2vec(MRIread(rCond.fPreprocList{I,1,1}));
    else
        rCond.volTs(I,1,1) = vol2vec(MRIread(rCond.fPreprocList{I,1,1}));
    end
    volMt(I) = doIt(rCond,W,K,winSec,dsgn,mask,extra,skipSVD,skipPSD,verbose,taperPerm,phaseRand,[],[],testFlag);
end

%% Compute on the average of all runs
rCond.volTs(1).vol = mean(cat(5,rCond.volTs.vol),5);
rCond.volTs(2:end) = [];
rCond.r = 1;
rCond.R = 1;
volMt(end+1) = doIt(rCond,W,K,winSec,dsgn,mask,extra,skipSVD,skipPSD,verbose,taperPerm,phaseRand,[],[],testFlag);





% %% Quick inspection of results
% close all
% ax   = {};
% axAv = {};
% for I = 0:size(rCond.fPreprocList,1)+1
%     figure('WindowStyle','docked');
%     if I==0 % cat across runs
%         f   = volMt(1).psd.f;
%         psd = cat(1,volMt(1:size(rCond.fPreprocList,1)).psd);
%         psd = cat(3,psd.PSD);
%     else
%         f   = volMt(I).psd.f;
%         psd = volMt(I).psd.PSD;
%     end
%     psd   = mean(psd   ,6); % average across voxels
%     psdAv = mean(psd   ,3); % average across runs
%     psdEr = std( psd,[],3); % error across runs

%     plot(squeeze(f),squeeze(psdAv),'k');

    
%     if I==0 % cat across runs
%         f   = volMt(1).psdTrialGramMD.f;
%         t   = volMt(1).psdTrialGramMD.t;
%         psd = cat(1,volMt(1:size(rCond.fPreprocList,1)).psdTrialGramMD);
%         psd = cat(6,psd.vec);
%         psd = cat(3,psd.psdPC);
%         psd = psd(:,:,:,:,:,:,end);
%     else
%         f   = volMt(I).psdTrialGramMD.f;
%         t   = volMt(I).psdTrialGramMD.t;
%         psd = volMt(I).psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,end);
%     end
%     psd   = mean(psd   ,6); % average across voxels
%     psdAv = mean(psd   ,3); % average across runs
%     psdEr = std( psd,[],3); % error across runs

    
%     hold on
%     plot(squeeze(f),squeeze(psdAv),'r');

%     if I==0
%         xline(mean(1./diff(volMt(1).param.dsgn.onsetList)).*(1:5),'g');

%         winSpan = mean(t - volMt(1).psdTrialGramMD.param.dsgn.onsetList,2); %sec
%         ws = winSpan(:,:,:,:,:,:,end);

%         axAv{1} = gca;
%         title('results averaged across runs');
%     elseif I==size(rCond.fPreprocList,1)+1
%         xline(mean(1./diff(volMt(I).param.dsgn.onsetList)).*(1:5),'g');

%         winSpan = mean(t - volMt(1).psdTrialGramMD.param.dsgn.onsetList,2); %sec
%         ws = winSpan(:,:,:,:,:,:,end);

%         axAv{end+1} = gca;
%         title(['timeseries averaged across runs']);
%     else
%         xline(mean(1./diff(volMt(I).param.dsgn.onsetList)).*(1:5),'g');

%         winSpan = mean(t - volMt(1).psdTrialGramMD.param.dsgn.onsetList,2); %sec
%         ws = winSpan(:,:,:,:,:,:,end);

%         ax{end+1} = gca;
%         title(['run ' num2str(I)]);
%     end
%     xlabel('Frequency (Hz)');
%     ylabel('PSD');
%     legend('full timeseries',[num2str(ws(1),'%0.1f') 's to ' num2str(ws(end),'%0.1f') 's post-stim onset (' num2str(mean(volMt(1).psdTrialGramMD.param.dsgn.ondurList),'%0.1f') 's stim dur)'],'stim fundamental and harmonics');
%     yscale('log')
%     grid on
% end
% yLim = get([ax{:}],'YLim'); yLim = cat(1,yLim{:}); yLim = [min(yLim(:,1)) max(yLim(:,2))];
% set([ax{:}],'YLim',yLim);
% yLim = get([axAv{:}],'YLim'); yLim = cat(1,yLim{:}); yLim = [min(yLim(:,1)) max(yLim(:,2))];
% set([axAv{:}],'YLim',yLim);

% axSave = [ax,axAv];
% outDir = fullfile(rCond.dirs.bidsDeriv,'..','..','forDavid');
% if ~exist(outDir,'dir'); mkdir(outDir); end
% for i = 1:length(axSave)
%     outFile = fullfile(outDir,replace(axSave{i}.Title.String,' ','_'));
%     saveas(axSave{i},[outFile '.fig']);
%     saveas(axSave{i},[outFile '.png']);
% end
% %% Output data for David





disp('!!!!!!!')
disp('!!!!!!!')
disp('!!!!!!!')
disp('!!!!!!!')
disp('!!!!!!!')
disp('continue the work here')
dbstack;
disp('!!!!!!!')
disp('!!!!!!!')
disp('!!!!!!!')
disp('!!!!!!!')
keyboard








funPsd = reshape(funPsd,size(volTs));


for I = 1:numel(volTs)
    fieldList = {'psd' 'psdGram' 'psdTrialGram' 'psdTrialGramMD' 'svd' 'svdGram' 'svdTrialGram' 'svdTrialGramMD'};
    for i = 1:length(fieldList)
        if isfield(funPsd(I),fieldList{i}) && ~isempty(funPsd(I).(fieldList{i}))
            if isfield(funPsd(I),'mask')
                funPsd(I).(fieldList{i}).param.mask = funPsd(I).mask;
            else
                funPsd(I).(fieldList{i}).param.mask = [];
            end
        end
    end
end





function volMt = doIt(rCond,W,K,winSec,dsgn,mask,extra,skipSVD,skipPSD,verbose,taperPerm,phaseRand,cohF,cohFperm,testFlag)
tpFlag = false;
if ~exist('winSec','var');        winSec = []; end
if ~exist('dsgn','var');            dsgn = []; end
onsetList = dsgn.onsetList';
durList   = dsgn.ondurList';
% if ~exist('onsetList','var'); onsetList = []; end
% if ~exist('durList','var');     durList = []; end
if ~exist('testFlag','var');   testFlag = []; end
if ~isfield(extra,'padTo'); extra.padTo = []; end
padTo   = extra.padTo;

if isempty(testFlag); testFlag = 0; end

if isempty(winSec); winSec = [inf 0]; end
windFlag = 0;

% try
if ~isempty(rCond.volTs(rCond.r).vol)
    nVox = prod(size(rCond.volTs(rCond.r).vol,[1 2 3]));
elseif isfield(rCond.mri,'vec') && ~isempty(rCond.mri.vec)
    nVox = size(rCond.mri.vec,2);
else
    nVox = rCond.mri.nvoxels;    
end
% catch
%     keyboard
% end
if nVox==1
    if verbose; disp('only one timeseries, skipping SVD'); end
    skipSVD = true;
end

param.dsgn = dsgn;
param.dsgn.winSec = winSec;
% param.win = win; % win always in seconds; param.win in seconds here, but will be converted to frames later
% param.onsetList = onsetList; % alwaysin seconds
% param.durList = durList; % alwaysin seconds
if ~isfield(param,'pad') || isempty(param.pad); param.pad = 0; end

if param.dsgn.winSec(1)==inf; skipGram = true; skipTrialGram = true; else skipGram = false; skipTrialGram = false; end
if isempty(param.dsgn.onsetList)
    skipTrialGram = true;
else
    if all((size(param.dsgn.onsetList)==1)==[0 1])
        dbstack; error('onsetList must be a row vector')
    end
end
skipTrialGramMD = skipTrialGram;


% if length(param.dsgn.win)>2
%     param.onsetList = param.dsgn.win(3:end)';
%     param.dsgn.win(3:end) = [];
% else
%     param.onsetList = [];
% end

% %% Load data
% if (~isfield(funTs,'vec') || isempty(funTs.vec)) && (~isfield(funTs,'vol') || isempty(funTs.vol))
%     funTs = MRIload3(funTs,mask,[],1);
% else

    %% Mask
    if ~isempty(mask)
        rCond.volTs(rCond.r) = vol2vec(rCond.volTs(rCond.r),mask,1);
        % funTs.mri = applyMask(funTs.mri,mask);
        % funTs.mask = mask;
    end
% end

% %% Assert
% if isfield(funTs,'vec') && ~isempty(funTs.vec); tmp = all(funTs.vec==0,1); else tmp = all(funTs.vol==0,4); end
% if any(tmp(:)); warning('Some voxels are all 0s. Adjust your mask to avoid later problems'); end

%% Complete some stuff
% if isfield(funTs,'tr') && length(funTs.tr)>1
%     funTs.tr = mean(funTs.tr);
% end
% if ~isfield(funTs,'tr') || isempty(funTs.tr) || isnan(funTs.tr)
%     if isfield(funTs,'Fs')
%         funTs.tr = 1/funTs.Fs(1) *1000;
%     else
%         dbstack; error('code that');
%     end
% end

tr = rCond.tr(rCond.r);
nDummyRemoved = rCond.nFrameOrig(rCond.r) - rCond.nFrame(rCond.r);
t = ((nDummyRemoved+1):rCond.nFrameOrig(rCond.r))-1; t = (t.*tr)';

% if ~isfield(funTs,'nDummyRemoved')
%     funTs.nDummyRemoved = 0;
% end
% if ~isfield(funTs,'t') || isempty(funTs.t) || any(isnan(funTs.t))
%     n = funTs.nframes;
%     s = tr.*(0+funTs.nDummyRemoved  );
%     e = tr.*(n+funTs.nDummyRemoved-1);
%     funTs.t = linspace(s,e,n)';
% end



%% Set parameters
Wflag = ~isempty(W);
Kflag = ~isempty(K);
% % tpFlag = ~isempty(tp); if tpFlag; Wflag = false; Kflag = false; end
% if ~isfield(funTs,'nruns'); funTs.nruns = 1; end
%%% Window size (defined in seconds up to here, then in number of frames)
if param.dsgn.winSec(1)==inf
    %%%% single-window over the full time series
    T = tr.*rCond.nFrame(rCond.r);
    param.dsgn.win(1) = rCond.nFrame(rCond.r);
    param.dsgn.win(2) = 0;
    % param.dsgn.winSec = [T 0];
else
    %%%% multiple time windows
    % Seconds to frames
    param.dsgn.win = round(param.dsgn.winSec./tr);
    if length(param.dsgn.win)==1
        param.dsgn.win(2) = 1;
    end
    T = param.dsgn.win(1)*tr;
end
param.dsgn.winSec = param.dsgn.win.*tr;
% if Wflag && Kflag
%     error('Cannot specify both W and K');
% elseif Wflag
%     [TW,W,K] = W2K(T,W);
% elseif Kflag || tpFlag
%     if tpFlag && verbose
%         disp('using precomputed tapers');
%     end
%     [TW,W,K] = K2W(T,K,verbose);
% end
% param.tapers = [TW K];




%%%%%%%%%%%%%%%%%%%
%% xgram parameters
%%%%%%%%%%%%%%%%%%%
%%%%%%% Hack: detect when there is a gap in the time vector, indicating
%%%%%%% multiple runs where concatenated. Extract that time vector of the
%%%%%%% first runs to compute windows
if isfield(rCond,'t') && ~isempty(rCond.t) && all(~isnan(rCond.t)) && max(abs(diff(diff(rCond.t))))/mode(diff(rCond.t))>1.1 && length(rCond.nFrame)>1
    hackFlag = 1;
    nFrame = rCond.nFrame(1);
    onsetList = onsetList(onsetList < rCond.t(nFrame));
else
    hackFlag = 0;
    % nFrame = funTs.nframes;
    nFrame = rCond.nFrame(rCond.r);
end
    
    

if ~skipGram
    param.dsgn.win(3) = ceil(nFrame / (param.dsgn.win(2))); % define in number of volume
    allWin = repmat(1:param.dsgn.win(1),[param.dsgn.win(3) 1]);
    allWin = allWin + (((1:param.dsgn.win(3))-1)*param.dsgn.win(2))'; % win x t
    allWin(any(allWin>nFrame,2),:) = [];
    allWin = allWin - allWin(end,end) + nFrame;
    % allWin(end+1,:) = (nFrame-param.dsgn.win(1)+1:nFrame)';
    param.dsgn.win(3) = [];
    % param.dsgn.win = param.dsgn.win.*tr;
    if verbose && param.dsgn.winSec(2)~=inf
        disp(['win(1) (window width): ' num2str(param.dsgn.winSec(1),'%0.3f') 'sec or ' num2str(param.dsgn.winSec(1)/tr) 'vol'])
        disp(['win(2) (step size)   : ' num2str(param.dsgn.winSec(2),'%0.3f') 'sec or ' num2str(param.dsgn.winSec(2)/tr) 'vol'])
    end
else
    allWin = [];
end

allWin = unique(allWin,'rows');
%% %%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% trial-locked xgram parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~skipTrialGram
    allWin; % [win X timeIndex]
    % onsetList = param.onsetList;
    if windFlag; onsetList(1) = []; end
    winSz = param.dsgn.win(1);
    n = nFrame+winSz;
    % n = max(allWin(:));
    nWin = size(allWin,1);
    nTrial = size(onsetList,1);
    try
        allWin2 = repmat({zeros(n,nWin)},[nTrial 1]);
    catch ME
        warning(ME.message)
        keyboard
    end
    for winInd = 1:nWin
        for trialInd = 1:nTrial
            if trialInd == 1
                allWin2{trialInd}(allWin(winInd,:),winInd) = 1;
            else
                % offsetInd = floor((onsetList(trialInd) - onsetList(1)) ./ tr);
                offsetInd = floor((onsetList(trialInd) - onsetList(1)) ./ tr) + 1;
                tInd = allWin(winInd,:) + offsetInd;
                tInd(tInd>n) = [];
                allWin2{trialInd}(tInd,winInd) = 1;
            end
        end
    end

    
    
    %%% remove windows exceeding timeseries
    allWin3 = any(cat(3,allWin2{:}),3); % [win X time in trial]
    % endInd = find(allWin3(end,:),1)+1;
    endInd = find(allWin3(nFrame,:),1)+1;
    % endInd = find(allWin3(end,:),1);
    for trialInd = 1:nTrial
        allWin2{trialInd}(:,endInd:end) = [];
    end
    % allWin3(:,endInd:end) = [];

    % %%% remove completely overlapping windows
    % endInd = find(sum((allWin2{1}+allWin2{2}(:,1))==2,1)==winSz);
    % for trialInd = 1:nTrial
    %     allWin2{trialInd}(:,endInd:end) = [];
    % end
    % % allWin3(:,endInd:end) = [];


    %%% remove completely overlapping windows
    [a,endInd] = max(sum(allWin2{end}==allWin2{end-1}(:,end),1));
    if a<n; endInd = endInd - 1; end
    for trialInd = 1:nTrial
        allWin2{trialInd}(:,1:endInd) = [];
    end
    % imagesc(any(cat(3,allWin2{:}),3))
    
    
    %%% allign windows to the end of the timeseries
    timeShiftVol = n-find(allWin2{end}(:,end),1,'last');
    % imagesc(any(cat(3,allWin2{:}),3))
    % find(allWin2{1}(:,end),1,'last')
    % find(allWin2{2}(:,end),1,'last')
    % find(allWin2{3}(:,end),1,'last')
    % find(allWin2{4}(:,end),1,'last')
    % find(allWin2{5}(:,end),1,'last')
    for trialInd = 1:nTrial
        allWin2{trialInd}(end-timeShiftVol+1:end,:) = false;
        allWin2{trialInd} = circshift(allWin2{trialInd},timeShiftVol,1);
    end

    %%% trim to data length
    for trialInd = 1:nTrial
        ind = size(allWin2{trialInd},1) - nFrame;
        allWin2{trialInd}(1:ind,:) = [];
    end

    %!!!!!!!!!
    %!!!!!!!!!
    %!!!!!!!!!
    %!!!!!!!!!
    %%% trim incomplete windows
    allWin3 = any(cat(3,allWin2{:}),3); % [win X time in trial]
    ind = strfind(diff(sum(allWin3,1)<mode(sum(allWin3,1))),[-1 0]);
    if ~isempty(ind)
        ind = ind(1);
        for trialInd = 1:nTrial
            allWin2{trialInd}(:,1:ind) = [];
        end
    end
    %!!!!!!!!!
    %!!!!!!!!!
    %!!!!!!!!!
    %!!!!!!!!!
    

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

    % figure('WindowStyle','docked');
    % imagesc(any(allWin3,3))
    % yline(param.onsetList./(funTs.tr/1000),'r')

    winSz2 = unique(sum(allWin3,1));
    if length(winSz2)~=1
        dbstack; warning([newline 'badly defined trial-locked windows' newline 'win too long for ISI?' newline 'anyway, skipping trialGram and trialGramMD'])
        skipTrialGram   = 1;
        skipTrialGramMD = 1;
        allWinTrialLock = [];
        % dbstack; error('badly defined trial-locked windows');
    else
        allWin4 = zeros(winSz2,size(allWin3,2),size(allWin3,3));
        for winInd = 1:prod(size(allWin3,[2 3]))
            allWin4(:,winInd) = find(allWin3(:,winInd));
        end
        allWinTrialLock = permute(allWin4,[2 1 3]); % [win X timeIndex]
        % allWinTrialLock = reshape(allWinTrialLock,size(allWinTrialLock,1),prod(size(allWinTrialLock,[2 3])));
    end
    clear allWin2 allWin3 allWin4 winSz2
else
    allWinTrialLock = [];
end

if hackFlag
    % restore to before the hack
    nFrame = rCond.nFrame;
    onsetList = param.onsetList;

    % repeat windows for subsequent runs
    nSum = cumsum(nFrame);
    nSum(end) = [];
    allWin2 = allWin;
    allWinTrialLock2 = allWinTrialLock;
    while ~isempty(nSum)
        allWin2 = cat(1,allWin2,allWin+nSum(1));
        allWinTrialLock2 = cat(3,allWinTrialLock2,allWinTrialLock+nSum(1));
        nSum(1) = [];
    end
    allWin          = allWin2         ; clear allWin2
    allWinTrialLock = allWinTrialLock2; clear allWinTrialLock2
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Detect impossible parameter combination and turn analysis off
% for i = 1:length(K)
%     [~,Wtmp,~] = K2W(param.win(1),K(i),0);
%     if Wtmp>0.5
%         warning(['cannot do gram with this parameter combination' newline 'skipping gram and trialGram'])
%         skipGram = 1;
%         skipTrialGram = 1;
%     end
% end


%% initiate stuff
% volMt = rmfield(rCond.volTs(rCond.r),{'vol' 'vec'});
% volMt.vol = []; volMt.vec = []; setNiceFieldOrder(volMt,{'volInfo' 'vecInfo' 'vol' 'vol2vec' 'vec' });

% funPsd.vecInfo = strsplit(funPsd.vecInfo,' x ');
% funPsd.vecInfo{ismember(funPsd.vecInfo,'time/freq')} = 'freq/time';
% funPsd.vecInfo = strjoin(funPsd.vecInfo,' x ');
% funPsd.volInfo = strsplit(funPsd.volInfo,' x ');
% funPsd.volInfo{ismember(funPsd.volInfo,'time/freq')} = 'freq/time';
% funPsd.volInfo = strjoin(funPsd.volInfo,' x ');

% [funPsd.vol] = deal([]);
% [funPsd.vec] = deal([]);
% if ~isfield(funTs.volTs,'volInfo'); [funTs.volTs.volInfo] = deal(strjoin({'X' 'Y' 'Z' 'time/freq' 'taper/mode' 'run'},' x ')); end
% if ~isfield(funTs.volTs,'vecInfo'); [funTs.volTs.vecInfo] = deal(strjoin({'time/freq' 'vox' 'taper/mode' 'run'},' x ')); end
% tmp = strsplit(funTs(1).vecInfo,' x '); tmp{1} = 'freq/time'; tmp = strjoin(tmp,' x ');
% strsplit(funPsd(1).vecInfo,' x ')
% funPsd.vecInfo
% [funPsd.vecInfo] = deal(tmp);
% tmp = strsplit(funTs(1).volInfo,' x '); tmp{4} = 'freq/time'; tmp = strjoin(tmp,' x ');
% [funPsd.volInfo] = deal(tmp);

param.Fs = 1/tr;
param.complex = 1;


for sInd = 1:length(rCond)
    %% Get tapers
    % if ~isfield(funTs(sInd),'t') || isempty(funTs(sInd).t)
    %     funTs(sInd).t = linspace(0,(funTs(sInd).nframes-1)*funTs(sInd).tr/1000,funTs(sInd).nframes)';
    %     funTs(sInd).t = funTs(sInd).t + funTs(sInd).nDummyRemoved*funTs(sInd).tr/1000;
    % end
    % if ~isfield(funPsd(sInd),'t') || isempty(funPsd(sInd).t)
    %     funPsd(sInd).t = funTs(sInd).t;
    % end
    if skipPSD
        dbstack; error('code that')
    else
        %%% full timeseries
        if verbose; disp('getting taper for full timeseries'); end
        Kcur = K(end);
        % tr  = funTs(sInd).tr/1000;
        N   = nFrame;
        pad = param.pad;
        if ~isempty(padTo) && ~isnan(padTo(3))
            pad = pad-1;
            NFFTtarg = (padTo(3)-1)*2;
            NFFT = max(2^(nextpow2(N)+pad),N);
            while NFFT<NFFTtarg
                pad = pad+1;
                NFFT = max(2^(nextpow2(N)+pad),N);
            end
            clear NFFT NFFTtarg
        end

        % !!!!!!!!! does not account for actual missing data (e.g. manually labeled artefact timepoints)
        [TP.full.tp,TP.full.eigs,TP.full.tpDC,TP.full.N,TP.full.pad] = getTapers(Kcur,tr,N,t,pad);
        TP.full.t = t;



        %%% time-resolved
        if ~skipGram
            if verbose; disp('getting taper for time-resolved analysis'); end
            Kcur = K(1);
            % tr  = funTs(sInd).tr/1000;
            N   = size(allWin,2);
            pad = param.pad;
            if ~isempty(padTo) && ~isnan(padTo(1))
                pad = pad-1;
                NFFTtarg = (padTo(1)-1)*2;
                NFFT = max(2^(nextpow2(N)+pad),N);
                while NFFT<NFFTtarg
                    pad = pad+1;
                    NFFT = max(2^(nextpow2(N)+pad),N);
                end
                clear NFFT NFFTtarg
            end
            % !!!!!!!!! does not account for actual missing data (e.g. manually labeled artefact timepoints)
            [TP.gram.tp,TP.gram.eigs,TP.gram.tpDC,TP.gram.N,TP.gram.pad] = getTapers(Kcur,tr,N,[],pad);
            TP.gram.t = t(allWin(1,:));
        else
            TP.gram.tp   = [];
            TP.gram.eigs = [];
            TP.gram.t    = [];
            TP.gram.N    = [];
            TP.gram.pad  = [];
        end

        %%% trial-locked
        if ~skipTrialGram
            %%%% Regular tapers repeated at each trial
            E = length(param.dsgn.onsetList);
            TP.trialGram.tp = repmat(TP.gram.tp,[1 1 E]);
            TP.trialGram.eigs = TP.gram.eigs;
            TP.trialGram.t = TP.gram.t;
            TP.trialGram.N = TP.gram.N;
            TP.trialGram.pad = TP.gram.pad;
        else
            TP.trialGram.tp     = [];
            TP.trialGram.eigs   = [];
            TP.trialGram.t      = [];
            TP.trialGram.N = [];
            TP.trialGram.pad = [];
        end
        if ~skipTrialGramMD
            %%%% Missing data tapers (experimental)
            if length(K)>1
                Kcur = K(2);
            else
                Kcur = K;
            end
            % tr  = funTs(sInd).tr/1000;
            % t   = funTs(sInd).t;
            % t   = t(allWinTrialLock(1,:));
            N   = length(t(allWinTrialLock(1,:)));
            if ~isfield(param,'pad') || isempty(param.pad); param.pad = 0; end
            pad = param.pad;
            if ~isempty(padTo) && ~isnan(padTo(2))
                pad = pad-1;
                NFFTtarg = (padTo(2)-1)*2;
                NFFT = max(2^(nextpow2(N)+pad),N);
                while NFFT<NFFTtarg
                    pad = pad+1;
                    NFFT = max(2^(nextpow2(N)+pad),N);
                end
                clear NFFT NFFTtarg
            end
            [TP.trialGramMD.tp,TP.trialGramMD.eigs,TP.trialGramMD.tpDC,TP.trialGramMD.N,TP.trialGramMD.pad] = getTapers(Kcur,tr,N,t(allWinTrialLock(1,:)),pad);
            TP.trialGramMD.t = t(allWinTrialLock(1,:));
        else
            TP.trialGramMD.tp   = [];
            TP.trialGramMD.eigs = [];
            TP.trialGramMD.t    = [];
            TP.trialGramMD.N    = [];
            TP.trialGramMD.pad  = [];
        end
        TP.full.info = 'time x taper x trial';
        TP.gram.info = 'time x taper x trial';
        TP.trialGram.info = 'time x taper x trial';
        TP.trialGramMD.info = 'time x taper x trial';
    end


    % figure('WindowStyle','docked')
    % ht = tiledlayout(1,2);
    % ax1 = nexttile;
    % plot(TP.full.t,TP.full.tp)
    % ylim([-0.15 0.15])
    % title(legend(num2str(TP.full.eigs')),'eigenvalues'); grid on; xlabel('time (s)')
    % ax2 = nexttile;
    % sz = [size(TP.full.t,1) size(TP.trialGramMD.tp,2)];
    % tp = nan(sz);
    % tp(allWinTrialLock(1,:),:) = TP.trialGramMD.tp;
    % plot(t,tp)
    % title(legend(num2str(TP.trialGramMD.eigs')),'eigenvalues'); grid on; xlabel('time (s)')
    % linkaxes([ax1 ax2])
    % title(ax1,'full timeseries')
    % title(ax2,['time-resolved leveraging missing data' newline '(showing a single post-stimulus window)'])
    % title(ht,'slepian tapers')



    %% Compute (also with taper-level permutation
    skip.psd         = skipPSD;
    skip.svd         = skipSVD;
    skip.gram        = skipGram;
    skip.trialGram   = skipTrialGram;
    skip.trialGramMD = skipTrialGramMD;
    paramInit = param;
    param.psd          = rmfield(paramInit,'dsgn');
    param.svd          = rmfield(paramInit,'dsgn');
    param.psdGram      = rmfield(paramInit,'dsgn');
    param.svdGram      = rmfield(paramInit,'dsgn');
    param.psdTrialGram = paramInit;
    param.svdTrialGram = paramInit;
    clear paramInit
    param.psdGram.winInd      = allWin;
    param.svdGram.winInd      = allWin;
    param.psdTrialGram.winInd = allWinTrialLock;
    param.svdTrialGram.winInd = allWinTrialLock;

    param.perm = taperPerm;

    if isfield(extra,'keepJ')
        param.keepJ = extra.keepJ;
    else
        param.keepJ = 0;
    end
    if isfield(extra,'Kf')
        param.Kf = extra.Kf;
    else
        param.Kf = [];
    end
    
    [...
        volMt.psd ,volMt.psdGram ,volMt.psdTrialGram ,volMt.psdTrialGramMD ,...
        volMt.svd ,volMt.svdGram ,volMt.svdTrialGram ,volMt.svdTrialGramMD ,...
        volMt.harm,volMt.harmGram,volMt.harmTrialGram,volMt.harmTrialGramMD,...
        volMt.svdXfreq]...
        = computeAll(rCond.volTs(rCond.r),TP,param,windFlag,skip,verbose,testFlag);


    %% Sort outputs
    volMt(sInd).param = param;


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


function [psd,psdGram,psdTrialGram,psdTrialGramMD,svd,svdGram,svdTrialGram,svdTrialGramMD,harm,harmGram,harmTrialGram,harmTrialGramMD,svdXfreq] = computeAll(volTs,TP,param,windFlag,skip,verbose,testFlag)
if ~exist('nShuf','var');         nShuf = []; end
if ~exist('nRun','var');           nRun = []; end
if ~exist('cohFrange','var'); cohFrange = []; end
if ~exist('testFlag','var');   testFlag = []; end

if isempty(nShuf);         nShuf = 0; end
if isempty(nRun);           nRun = 1; end
if isempty(cohFrange); cohFrange = [0 inf]; end
if isempty(testFlag);   testFlag = 0; end
dtrndTsFlag  = 1;
dtrndWinFlag = 0;
dtrndWinOrd  = 0;




for runInd = 1:nRun
    if verbose && nRun>1; disp(['---Run ' num2str(runInd) '/' num2str(nRun) '---']); end

    Fs = param.Fs;
    if ~isfield(volTs,'vec') && isfield(volTs,'vol') && isempty(volTs.vol) && isfield(volTs,'vol2vec') && ~isempty(volTs.vol2vec)
        volTs = MRIload2(volTs);
    end
    volTs = vol2vec(volTs);


    %%%%%%%%%%
    %% Detrend
    %%%%%%%%%%
    if dtrndTsFlag
        volTs = dtrnd2(volTs);
    end
    %% %%%%%%%
    
    

    %%%%%%%%%%%%%%%%%%%%%%%
    %% Over full timeseires
    %%%%%%%%%%%%%%%%%%%%%%%
    %[time x trial x run x taper x freq x vox x window x mode]
    % tic
    disp('full timeseries analysis')
    %%% tapers
    tp   = permute(TP.full.tp  ,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
    tpDC = permute(TP.full.tpDC,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
    [Nk,Ek,Rk,Kk,Fk,Vk,Wk,M] = size(tp);
    K = Kk;
    N = Nk;

    %%% windows and trials
    W = 1;
    E = 1;

    %%% time
    if isfield(volTs,'t') && ~isempty(volTs.t)
        t = volTs.t; % tapers[time x trial x run x taper x freq x vox x window x mode]
    else
        dbstack; error('X');
        t = permute(linspace(0,(N-1)/Fs,N),[2 1 3 4 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
    end
    [Nt,Et,Rt,Kt,Ft,Vt,Wt,Mt] = size(t);

    %%% freq
    if isfield(param,'Kf') && ~isempty(param.Kf)
        f = permute(param.Kf,[1 3 4 5 2]);
    else
        if ~isfield(param,'pad') || isempty(param.pad)
            pad = 0;
        else
            pad = param.pad;
        end
        if isfield(TP.full,'pad') || ~isempty(TP.full.pad)
            if pad~=TP.full.pad; disp('!!!'); warning('overriding param.pad with TP.full.pad'); end
            pad = TP.full.pad;
        end
        NFFT=max(2^(nextpow2(N)+pad),N);
        [f,~]=getfgrid(Fs,NFFT,[0 Fs/2]);
        f = permute(f,[1 3 4 5 2 6 7 8]); % frequencies[time x trial x run x taper x freq x vox x window x mode]
    end
    [Nf,Ef,Rf,Kf,Ff,Vf,Wf] = size(f);
    F = Ff;

    %%% channels and runs
    [~,V,~,R] = size(volTs.vec);

    %%% modes
    M = min([V K]);
    Mxfreq = min([V K*F]);

    %%% allocate
    if param.keepJ
        res.full.J = zeros(1,E,1,1,F,V,W,1  ); % psd       [time x trial x run x taper x freq x vox x window x mode] at each trial
    end
    res.full.PSD   = zeros(1,E,1,1,F,V,W,1  ); % psd       [time x trial x run x taper x freq x vox x window x mode] at each trial
    res.full.COH   = zeros(1,E,1,1,F,1,W,M  ); % coherence [time x trial x run x taper x freq x vox x window x mode] at each trial

    %%% Compute J
    d = volTs.vec;              % [time  x vox x taper x run               ]
    d = permute(d,[3 2 4 1]);   % [taper x vox x run   x time*trial        ]
    d = reshape(d,[1 V R N E]); % [taper x vox x run   x time       x trial]
    d = permute(d,[4 5 3 1 6 2 7 8]);
    tp = tp;
    tt = t;
    J = getJ4(d,tp,tt,f)/Fs;% [N E R K F V W]
    % J2 = getJ4(d,tp,tt,f,1)/Fs; figure('WindowStyle','docked'); scatter(J2(:),J(:)); ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; grid on; max(abs(real(J(:) - J2(:))))./median(abs(real(J(:))))
    %%% Compute psd
    if param.keepJ
        res.full.J = J;
    end
    res.full.PSD   = mean(  conj(J).*J  ,4);


    %%% Compute line power (adapted from Chronux)
    % %%%% at all frequencies
    % [linePwr,lineF,lineP] = getLinePwr(J,tp,Fs);
    %%%% at specified frequencies
    if isfield(param,'onsetList') && ~isempty(param.dsgn.onsetList)
        if max(abs(diff(diff(param.dsgn.onsetList))))>0.000001; keyboard; end
        fStim = 1/mean(diff(param.dsgn.onsetList));
        hInd = 1:min(floor(Fs/2/fStim),10);
        fStim = permute(fStim.*hInd,[1 3 4 5 2 6 7]);
        Jstim = getJ4(d,tp,tt,fStim,1)/Fs; % [time x trial x run x taper x freq x vox x window] % [N E R K F V W]
        
        % JstimX = getJ4(d,tp,tt,fStim)/Fs; % [time x trial x run x taper x freq x vox x window] % [N E R K F V W]
        % Jx = getJ4(d,tp,tt,f,1)/Fs;% [N E R K F V W]
        % J1 = getJ4(d,tp,tt,fStim,1)/Fs; % [time x trial x run x taper x freq x vox x window] % [N E R K F V W]
        % J0 = getJ4(d,tp,tt,fStim,0)/Fs; % [time x trial x run x taper x freq x vox x window] % [N E R K F V W]
        
        [linePwrStim,lineFstim,linePstim] = getLinePwr(Jstim,tp,Fs);

        res.full.harm.label   = 'stimFreq';
        res.full.harm.linePwr = linePwrStim; clear linePwrStim
        res.full.harm.lineF   = lineFstim;   clear lineFstim
        res.full.harm.lineP   = linePstim;   clear linePstim
        res.full.harm.f       = fStim;       clear fStim
    else
        res.full.harm.label   = '';
        res.full.harm.linePwr = [];
        res.full.harm.lineF   = [];
        res.full.harm.lineP   = [];
        res.full.harm.f       = [];
    end
    

    % %%% Reconstruct harmonic signal fit
    % %%%% at prespecified frequencies (same for all voxels)
    % fInd = [];
    % [~,fInd(end+1)] = min(abs(f-0.1));
    % [~,fInd(end+1)] = min(abs(f-0.23));
    % ts = getLineTs(A,f,fInd,[],N,Fs);
    % %%%% at voxel-specific frequencies (based on p-values and threshold alpha)
    % ts = getLineTs(A,f,0.05,p,N,Fs);
    

    %%% Compute coherence
    if K > 1
        dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
        prm = [ 6   4   2   1   3   5   7   8  ];
        dim = strjoin(dim(prm),' ');
        % j = permute(J,prm); %[V K E N R F W Mk]
        % j = reshape(j,[V K E*1*R*F*1*1]); %[V K E*N*R*F*W*Mk]
        % [u,s,v] = pagesvd(j,'econ','vector'); % s[Mk V E*N*R*F*W]
        [u,s,~] = pagesvd(reshape(permute(J,prm),[V K E*1*R*F*1*1]),'econ','vector'); % s[Mk V E*N*R*F*W]
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


        
        % Cross-frequency mt-svd filter
        if isfield(param,'Kf') && ~isempty(param.Kf)
            % catenate different frequencies as extra columns for the svd (instead of different svd for each frequencies)
            dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
            prm = [ 6   4   5   2   1   3   7   8  ];
            dim = strjoin(dim(prm),' ');
            [u,s,~] = pagesvd(reshape(permute(J,prm),[V K*F E*1*R*1*1]),'econ','vector'); % s[Mk V E*N*R*F*W]
            coh     = s.^2./sum(s.^2,1); % coherence[Mk V E*N*R*W]
            dim = {'Mk' 'V' 'E' 'N' 'R' 'F' 'W' 'K'};
            prm = [4 3 5 8 6 2 7 1];
            res.xfreq.COH  = permute(reshape(coh,[Mxfreq 1 E 1 R 1 1]),prm); %[time x trial x run x taper x freq x vox x window x mode]
            res.xfreq.spSV = permute(u,[4 5 6 7 3 1 8 2]); %[time x trial x run x taper x freq x vox x window x mode]
        end


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

    if isfield(param,'Kf') && ~isempty(param.Kf)
        res.xfreq.f     = permute(param.Kf,[1 3 4 5 2]);
        res.xfreq.K     = K;
        res.xfreq.T     = N/Fs;
        res.xfreq.E      = [];
        res.xfreq.win    = [];
        res.xfreq.param  = rmfield(param,{'complex' 'psd' 'svd' 'psdGram' 'svdGram' 'psdTrialGram' 'svdTrialGram'});
        res.xfreq.dim     = [N E R K F V W Mxfreq];
        res.xfreq.dimInfo = '[N E R K F V W M]';    
    end
    %% %%%%%%%%%%%%%%%%%%%%

    



    %%%%%%%%%%%%%%%%%%%%%%%
    %% Over each timewindow
    %%%%%%%%%%%%%%%%%%%%%%%
    % tic
    disp('time-resolved analysis')
    if ~skip.gram
        %[time x trial x run x taper x freq x vox x window x mode]
        %[   7       2     1       3      5     8       20]
        tp   = permute(TP.gram.tp  ,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
        tpDC = permute(TP.gram.tpDC,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
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
        if ~isfield(param,'pad') || isempty(param.pad)
            pad = 0;
        else
            pad = param.pad;
        end
        if isfield(TP.gram,'pad') || ~isempty(TP.gram.pad)
            if pad~=TP.gram.pad; disp('!!!'); warning('overriding param.pad with TP.gram.pad'); end
            pad = TP.gram.pad;
        end
        NFFT=max(2^(nextpow2(N)+pad),N);
        [f,fInd]=getfgrid(Fs,NFFT,[0 Fs/2]);
        f = permute(f,[1 3 4 5 2 6 7 8]); % frequencies[time x trial x run x taper x freq x vox x window x mode]
        [Nf,Ef,Rf,Kf,Ff,Vf,Wf] = size(f);
        F = Ff;

        %%% channels and runs
        [~,V,~,R] = size(volTs.vec);

        %%% modes
        M = min([V K]);

        %%% allocate
        res.gram.PSD    = zeros(1,E,1,1,F,V,W,1  ); % psd       [time x trial x run x taper x freq x vox x window x mode] at each trial
        res.gram.COH    = zeros(1,E,1,1,F,1,W,M  ); % coherence [time x trial x run x taper x freq x vox x window x mode] at each trial

        %%% loop over windows
        if verbose>1
            fprintf([repmat('|',1,W) '\n\n']);
        end
        for wInd = 1:W
            %%% Compute J
            ind  = w(:,:,:,:,:,:,wInd);
            tWin = reshape(volTs.t(ind,:,:,:),size(ind)); tWin = tWin([1 end],:);
            d = volTs.vec(ind,:,:,:);   % [time  x vox x taper x run               ]
            d = permute(d,[3 2 4 1]);   % [taper x vox x run   x time*trial        ]
            d = reshape(d,[1 V R N E]); % [taper x vox x run   x time       x trial]
            d = permute(d,[4 5 3 1 6 2 7 8]); %[time x trial x run x taper x freq x vox x window x mode]
            if dtrndWinFlag
                d = detrend(d,dtrndWinOrd);
                % d = d - mean(d,1);
            end
            tp = tp;
            tt = t; % adjust t here for sub-tr stimulus onsets

            J = getJ4(d,tp,tt,f,[],testFlag)/Fs; % [N E R K F V W]
            
            %%% Compute psd at each trial
            res.gram.PSD(:,:,:,:,:,:,wInd)    = mean(  conj(J).*J  ,4);

            %%% Compute coherence at each trial
            if K > 1
                dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
                prm = [ 6   4   2   1   3   5   7   8  ];
                dim = strjoin(dim(prm),' ');
                % j = permute(J,prm); %[V K E N R F W Mk]
                % j = reshape(j,[V K E*1*R*F*1*1]); %[V K E*N*R*F*W*Mk]
                % [u,s,v] = pagesvd(j,'econ','vector'); % s[Mk V E*N*R*F*W]
                [~,s,~] = pagesvd(reshape(permute(J,prm),[V K E*1*R*F*1*1]),'econ','vector'); % s[Mk V E*N*R*F*W]
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
            if verbose>1
                fprintf('\b''\n');
            end
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

    % toc
    %% %%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Over each event-related timewindow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Much time could be saved by reusing the output of the time-resolved
    % analysis, as long as the step sizes allow to construct every
    % event-related window from available time-resolved windows. However,
    % the missing-data event-related analysis seems superior to the
    % event-related analysis so the latter might not be needed.
    % tic
    disp('trial-locked time-resolved analysis')
    %[time x trial x run x taper x freq x vox x window x mode]
    %[   7       2     1       3      5     8       20]
    if ~skip.trialGram
        %%% tapers
        tp   = permute(TP.gram.tp  ,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
        tpDC = permute(TP.gram.tpDC,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
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
        if ~isfield(param,'pad') || isempty(param.pad)
            pad = 0;
        else
            pad = param.pad;
        end
        if isfield(TP.trialGram,'pad') || ~isempty(TP.trialGram.pad)
            if pad~=TP.trialGram.pad; disp('!!!'); warning('overriding param.pad with TP.trialGram.pad'); end
            pad = TP.trialGram.pad;
        end
        NFFT=max(2^(nextpow2(N)+pad),N);
        [f,fInd]=getfgrid(Fs,NFFT,[0 Fs/2]);
        f = permute(f,[1 3 4 5 2 6 7 8]); % frequencies[time x trial x run x taper x freq x vox x window x mode]
        [Nf,Ef,Rf,Kf,Ff,Vf,Wf,Mf] = size(f);
        F = Ff;


        %%% channels and runs
        [~,V,~,R] = size(volTs.vec);

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
        if verbose>1
            fprintf([repmat('|',1,W) '\n\n']);
end
        
        for wInd = 1:W
            %%% Compute J
            ind  = w(:,:,:,:,:,:,wInd,:);
            tWin = reshape(volTs.t(ind,:,:,:),size(ind)); tWin = tWin([1 end],:);
            d = volTs.vec(ind,:,:,:);   % [time  x vox x taper x run               ]
            d = permute(d,[3 2 4 1]);   % [taper x vox x run   x time*trial        ]
            d = reshape(d,[1 V R N E]); % [taper x vox x run   x time       x trial]
            d = permute(d,[4 5 3 1 6 2 7 8]); %[time x trial x run x taper x freq x vox x window x mode]
            if dtrndWinFlag
                d = detrend(d,dtrndWinOrd);
            end
            tp = tp;

            % adjust t here for sub-tr stimulus onsets
            onsetList = param.dsgn.onsetList;
            if windFlag; onsetList(1) = []; end
            tt = t + tWin(1,:) - onsetList;

            J = getJ4(d,tp,tt,f,[],testFlag)/Fs; % [N E R K F V W]

            
            %%% Compute psd at each trial
            res.trialGram.PSD(:,:,:,:,:,:,wInd,:)    = mean(  conj(J).*J  ,4); % Compute power then average across tapers...
            %%% Compute psd averaged across trials
            res.trialGram.PSDeav(:,:,:,:,:,:,wInd,:)   = mean(  res.trialGram.PSD(:,:,:,:,:,:,wInd,:)  ,2); % then average across trials.
            %%% Compute psd phase-coherently averaged across trials
            res.trialGram.PSDepc(:,:,:,:,:,:,wInd,:) = mean(  conj(mean(J,2)).*mean(J,2)  ,4); % Average across trials then compute power then average across tapers.


            if K > 1
                %%% Compute coherence at each trial
                dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'M'};
                prm = [ 6   4   2   1   3   5   7   8  ];
                dim = strjoin(dim(prm),' ');
                % j = permute(J,prm); %[V K E N R F W M]
                % j = reshape(j,[V K E*1*R*F*1*1]); %[V K E*N*R*F*W*M]
                % [u,s,v] = pagesvd(j,'econ','vector'); % s[M V E*N*R*F*W]
                [~,s,~] = pagesvd(reshape(permute(J,prm),[V K E*1*R*F*1*1]),'econ','vector'); % s[M V E*N*R*F*W]
                coh = s.^2./sum(s.^2,1); % coherence[M V E*N*R*F*W]
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
                % j = permute(mean(J,2),prm); %[V K E N R F W Mk]
                % j = reshape(j,[V K 1*1*R*F*1*1]); %[V K E*N*R*F*W*M]
                % [u,s,v] = pagesvd(j,'econ','vector'); % s[M V E*N*R*F*W]
                [~,s,~] = pagesvd(reshape(permute(mean(J,2),prm),[V K 1*1*R*F*1*1]),'econ','vector'); % s[M V E*N*R*F*W]
                coh = s.^2./sum(s.^2,1); % coherence[M V E*N*R*F*W]
                coh = reshape(coh,[M 1 1 1 R F 1]); % [M V E N R F W]
                dim = {'Mk' 'V' 'E' 'N' 'R' 'F' 'W' 'K'};
                prm = [ 2    4   3   6   5   7   8   1 ];
                dim = strjoin(dim(prm),' ');
                coh = permute(coh,prm); %[V N E F R W K M]
                res.trialGram.COHepc(:,:,:,:,:,:,wInd,:) = coh;
            end

            %%% Compute coherence with trials concatenated as extra sets of tapers (equivalent to averaging across trials)
            % dim1 = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
            % prm1 = [ 6   4   2   1   3   5   7   8  ];
            % dim1 = strjoin(dim1(prm1),' ');
            % j = permute(J,prm1); %[V K E N R F W Mk]
            % j = reshape(j,[V K E 1*R*F*1*1]); %[V K E N*R*F*W*Mk]
            % dim2 = {'V' 'K' 'E' 'N*R*F*W*Mk'};
            % prm2 = [ 1   4   2   3          ];
            % dim2 = strjoin(dim2(prm2),' ');
            % j = permute(j,prm2); %[V N*R*F*W*Mk K E]
            % j = reshape(j,[V 1*R*F*1*1 K*E]); %[V N*R*F*W*Mk K*E]
            % j = permute(j,[1 3 2]); %[V K*E N*R*F*W*Mk];
            dim1 = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
            prm1 = [ 6   4   2   1   3   5   7   8  ];
            dim1 = strjoin(dim1(prm1),' ');
            dim2 = {'V' 'K' 'E' 'N*R*F*W*Mk'};
            prm2 = [ 1   4   2   3          ];
            dim2 = strjoin(dim2(prm2),' ');
            [~,s,~] = pagesvd(permute(reshape(permute(reshape(permute(J,prm1),[V K E 1*R*F*1*1]),prm2),[V 1*R*F*1*1 K*E]),[1 3 2]),'econ','vector'); % s[Mke V N*R*F*W K E]
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
            if verbose>1
                fprintf('\b''\n');
            end
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
        res.trialGram.onsetList = permute(param.dsgn.onsetList,[2 1 3 4 5 6 7 8]);
        res.trialGram.ondurList = permute(param.dsgn.ondurList  ,[2 1 3 4 5 6 7 8]);
        res.trialGram.param     = rmfield(param,{'complex' 'psd' 'svd' 'psdGram' 'svdGram' 'psdTrialGram' 'svdTrialGram'});
    else
        res.trialGram = [];
    end

    % toc
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Over each event-related timewindow
    % (with missing data seperating tials)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE that time vector is manipulated only for the "PC" version of psd
    % and coh, so for the "nonPC", this is really just a full-run analysis
    % with missing data at regular intervals
    
    % tic
    disp('trial-locked time-resolved analysis (using missing data tapers)')
    %[time x trial x run x taper x freq x vox x window x mode]
    if ~skip.trialGramMD
        %%% tapers
        tp   = permute(TP.trialGramMD.tp  ,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
        tpDC = permute(TP.trialGramMD.tpDC,[1 3 4 2 5 6 7 8]); % tapers[time x trial x run x taper x freq x vox x window x mode]
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
        onsets = param.psdTrialGram.dsgn.onsetList;
        if windFlag
            onsets(1) = [];
        end
        


        %%% freq
        if ~isfield(param,'pad') || isempty(param.pad)
            pad = 0;
        else
            pad = param.pad;
        end
        if isfield(TP.trialGramMD,'pad') || ~isempty(TP.trialGramMD.pad)
            if pad~=TP.trialGramMD.pad; disp('!!!'); warning('overriding param.pad with TP.trialGramMD.pad'); end
            pad = TP.trialGramMD.pad;
        end
        NFFT=max(2^(nextpow2(N)+pad),N);
        [f,fInd]=getfgrid(Fs,NFFT,[0 Fs/2]);
        f = permute(f,[1 3 4 5 2 6 7 8]); % frequencies[time x trial x run x taper x freq x vox x window x mode]
        [Nf,Ef,Rf,Kf,Ff,Vf,Wf] = size(f);
        F = Ff;

        %%% channels and runs
        [~,V,~,R] = size(volTs.vec);

        %%% modes
        M = min([V K]);
        
        %%% allocate
        res.trialGramMD.PSD    = zeros(1,1,1,1,F,V,W,1); % psd       [time x trial x run x taper x freq x vox x window x mode] phase-coherently averaged across trials     (eVENT pHASE cOHERENT)
        res.trialGramMD.PSDepc = zeros(1,1,1,1,F,V,W,1); % psd       [time x trial x run x taper x freq x vox x window x mode] phase-coherently averaged across trials     (eVENT pHASE cOHERENT)
        res.trialGramMD.COH    = zeros(1,1,1,1,F,1,W,M); % coherence [time x trial x run x taper x freq x vox x window x mode] trials concatenated as extra sets of tapers (eVENT AS TAPERS k   )
        res.trialGramMD.COHepc = zeros(1,1,1,1,F,1,W,M); % coherence [time x trial x run x taper x freq x vox x window x mode] trials concatenated as extra sets of tapers (eVENT AS TAPERS k   )

        %%% loop over windows
        if verbose>1
            fprintf([repmat('|',1,W) '\n\n']);
end
        
        for wInd = 1:W
            %%% Compute J
            ind  = w(:,:,:,:,:,:,wInd);
            tWin = reshape(volTs.t(ind,:,:,:),size(ind)); tWin = tWin([1 end],:);
            d = permute(volTs.vec(ind,:,:,:),[2 3 4 1]); % [vox x taper x run x time]
            d = reshape(d,[V 1 R Nw E]); % [vox x taper x run x time x trial]
            if dtrndWinFlag
                d = permute(detrend(permute(d,[4 1 2 3 5 6 7 8]),dtrndWinOrd),[2 3 4 1 5 6 7 8]); % detrend on a trial-by-trial basis
            end
            d = reshape(d,[V 1 R Nw*E]); % [vox x taper x run x time]
            d = permute(d,[4 5 3 2 6 1 7 8]); % [time x trial x run x taper x freq x vox x window]
            tp = tp;

            
            %%% Using normal time (incoherent phase averaging across trials, but phase coherency can be generated by stimulus design)
            tt = t;
            
            ttTmp = permute(reshape(permute(tt ,[3 4 5 6 7 8 1 2]),[R 1 1 1 1 1 N/E E]),[7 8 1 2 3 4 5 6]);
            dTmp  = permute(reshape(permute(d ,[3 4 5 6 7 8 1 2]),[R 1 1 V 1 1 N/E E]),[7 8 1 2 3 4 5 6]);
            tpTmp = permute(reshape(permute(tp,[3 4 5 6 7 8 1 2]),[1 K 1 1 1 1 N/E E]),[7 8 1 2 3 4 5 6]);
            fTmp  = permute(reshape(permute(f ,[3 4 5 6 7 8 1 2]),[1 1 F 1 1 1 1   1]),[7 8 1 2 3 4 5 6]);
            J = getJ4(dTmp,tpTmp,ttTmp,fTmp,[],testFlag)/Fs; % [N E R K F V W]
            J = sum(J,2);


            %%%% Compute psd
            res.trialGramMD.PSD(:,:,:,:,:,:,wInd)    = mean(  conj(J).*J  ,4);
            
            %%%% Compute coherence
            if K > 1
                dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
                prm = [ 6   4   2   1   3   5   7   8  ];
                dim = strjoin(dim(prm),' ');
                j = permute(J,prm); %[V K E N R F W M]
                j = reshape(j,[V K 1*1*R*F*1*1]); %[V K E*N*R*F*W*Mk]
                [~,s,~] = pagesvd(j,'econ','vector'); % s[M V E*N*R*F*W]
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
            tt = reshape(reshape(t,Nw,E) - onsets,Nw*E,1); % aligning the time vector to 0 at stimulus onset effectively enforces coherent phase averaging across trials

            ttTmp = permute(reshape(permute(tt ,[3 4 5 6 7 8 1 2]),[R 1 1 1 1 1 N/E E]),[7 8 1 2 3 4 5 6]);
            dTmp  = permute(reshape(permute(d ,[3 4 5 6 7 8 1 2]),[R 1 1 V 1 1 N/E E]),[7 8 1 2 3 4 5 6]);
            tpTmp = permute(reshape(permute(tp,[3 4 5 6 7 8 1 2]),[1 K 1 1 1 1 N/E E]),[7 8 1 2 3 4 5 6]);
            fTmp  = permute(reshape(permute(f ,[3 4 5 6 7 8 1 2]),[1 1 F 1 1 1 1   1]),[7 8 1 2 3 4 5 6]);
            J = getJ4(dTmp,tpTmp,ttTmp,fTmp,[],testFlag)/Fs; % [N E R K F V W]
            J = sum(J,2);

            %%%% Compute psd
            res.trialGramMD.PSDepc(:,:,:,:,:,:,wInd)    = mean(  conj(J).*J  ,4);
            
            %%%% Compute coherence
            if K > 1
                dim = {'N' 'E' 'R' 'K' 'F' 'V' 'W' 'Mk'};
                prm = [ 6   4   2   1   3   5   7   8  ];
                dim = strjoin(dim(prm),' ');
                j = permute(J,prm); %[V K E N R F W M]
                j = reshape(j,[V K 1*1*R*F*1*1]); %[V K E*N*R*F*W*Mk]
                [~,s,~] = pagesvd(j,'econ','vector'); % s[M V E*N*R*F*W]
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
            if verbose>1
                fprintf('\b''\n');
            end
        end
        if K == 1
            res.trialGramMD.COH    = [];
            res.trialGramMD.COHepc = [];
        end

        res.trialGramMD.f         = f;
        res.trialGramMD.K         = K;
        res.trialGramMD.T         = N/Fs;
        res.trialGramMD.E         = E;
        res.trialGramMD.win       = [mean(reshape(diff(res.trialGramMD.t,[],1),[],1)) mean(reshape(diff(res.trialGramMD.t,[],7),[],1))];
        res.trialGramMD.onsetList = permute(param.dsgn.onsetList,[2 1 3 4 5 6 7 8]);
        res.trialGramMD.ondurList = permute(param.dsgn.ondurList  ,[2 1 3 4 5 6 7 8]);
        res.trialGramMD.param     = rmfield(param,{'complex' 'psd' 'svd' 'psdGram' 'svdGram' 'psdTrialGram' 'svdTrialGram'});
    else
        res.trialGramMD = [];
    end


    % toc
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end


%% Refactor (for reason that don't seem to make sense anymore)

%%% psd
if skip.psd
    psd  = [];
    harm = [];
else
    fields = {'COH' 'COH_pVal' 'COH_fdr' 'spSV' 'spSV_pVal' 'spSV_fdr'}; fields = fields(ismember(fields,fieldnames(res.full))); 
    psd = rmfield(res.full,fields);
    psd.info = 'time x trial x run x taper x freq x vox x window x mode';
    if isfield(res.full,'harm') && ~isempty(res.full.harm)
        harm = res.full.harm;
    else
        harm = [];
    end
end

%%% psdGram
if skip.gram || skip.psd
    psdGram  = [];
    harmGram = [];
else
    fields = {'COH' 'COH_pVal' 'COH_fdr' 'spSV' 'spSV_pVal' 'spSV_fdr'}; fields = fields(ismember(fields,fieldnames(res.gram)));
    psdGram = rmfield(res.gram,fields);
    psdGram.info = 'time x trial x run x taper x freq x vox x window x mode';
    if isfield(res.gram,'harm') && ~isempty(res.gram.harm)
        harmGram = res.gram.harm;
    else
        harmGram = [];
    end
end

%%% psdTrialGram
if skip.trialGram || skip.psd
    psdTrialGram  = [];
    harmTrialGram = [];
else
    fields = {'PSD' 'PSDeav' 'PSDepc' 'COH' 'COHeav' 'COHepc' 'COHek'}; fields = fields(ismember(fields,fieldnames(res.trialGram)));
    psdTrialGram = rmfield(res.trialGram,fields);
    psdTrialGram.vec.psd   = res.trialGram.PSDeav;
    psdTrialGram.vec.psdPC = res.trialGram.PSDepc;
    psdTrialGram.info = 'time x trial x run x taper x freq x vox x window x mode';
    
    if isfield(res.trialGram,'harm') && ~isempty(res.trialGram.harm)
        harmTrialGram = res.trialGram.harm;
    else
        harmTrialGram = [];
    end
end

%%% psdTrialGramMD
if skip.trialGramMD || skip.psd
    psdTrialGramMD  = [];
    harmTrialGramMD = [];
else
    fields = {'PSD' 'PSDepc' 'COH' 'COHepc'}; fields = fields(ismember(fields,fieldnames(res.trialGramMD)));
    psdTrialGramMD = rmfield(res.trialGramMD,fields);
    psdTrialGramMD.vec.psd   = res.trialGramMD.PSD;
    psdTrialGramMD.vec.psdPC = res.trialGramMD.PSDepc;
    psdTrialGramMD.info = 'time x trial x run x taper x freq x vox x window x mode';
    if isfield(res.trialGramMD,'harm') && ~isempty(res.trialGramMD.harm)
        harmTrialGramMD = res.trialGramMD.harm;
    else
        harmTrialGramMD = [];
    end
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
else
    fields = {'PSD' 'PSDeav' 'PSDepc' 'COH' 'COHeav' 'COHepc' 'COHek'}; fields = fields(ismember(fields,fieldnames(res.trialGram)));
    svdTrialGram = rmfield(res.trialGram,fields);
    svdTrialGram.vec.coh    = res.trialGram.COHeav;
    svdTrialGram.vec.cohEPC = res.trialGram.COHepc;
    svdTrialGram.vec.cohEK  = res.trialGram.COHek;
    svdTrialGram.info = 'time x trial x run x taper x freq x vox x window x mode';
end

if skip.trialGramMD || skip.svd
    svdTrialGramMD = [];
else
    fields = {'PSD' 'PSDepc' 'COH' 'COHepc'}; fields = fields(ismember(fields,fieldnames(res.trialGramMD)));
    svdTrialGramMD = rmfield(res.trialGramMD,fields);
    svdTrialGramMD.vec.coh    = res.trialGramMD.COH;
    svdTrialGramMD.vec.cohEPC = res.trialGramMD.COHepc;
    svdTrialGramMD.info = 'time x trial x run x taper x freq x vox x window x mode';
end

if isfield(param,'Kf') && ~isempty(param.Kf)
    svdXfreq = res.xfreq;
    svdXfreq.info = 'time x trial x run x taper x freq x vox x window x mode';
else
    svdXfreq = [];
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
