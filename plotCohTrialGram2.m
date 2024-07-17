function [ax,F] = plotCohTrialGram2(volPsd,H,mFlag)
if ~exist('H','var');         H = [];                             end
if isempty(H);                H = figure('WindowStyle','docked'); end
if ~exist('mFlag','var'); mFlag = [];                             end
if isempty(mFlag);        mFlag = 'av'; end

switch class(H)
    case 'matlab.graphics.layout.TiledChartLayout'
        F = H.Parent;
    case 'matlab.ui.Figure'
        F = H;
    otherwise
end

figure(F);
ax = {};
ax{end+1} = nexttile;

%% spatial average
switch mFlag
    case 'av'
        if isempty(volPsd.svdTrialGram.vec.coh); return; end
        c = permute(volPsd.svdTrialGram.vec.coh(:,:,:,:,:,:,:,1),[5 7 2 8 1 3 4 6]);
    case 'pc'
        if isempty(volPsd.svdTrialGram.vec.cohEPC); return; end
        c = permute(volPsd.svdTrialGram.vec.cohEPC(:,:,:,:,:,:,:,1),[5 7 2 8 1 3 4 6]);
    case 'ek'
        c = permute(volPsd.svdTrialGram.vec.cohEK(:,:,:,:,:,:,:,1),[5 7 2 8 1 3 4 6]);
    case 'MD'
        c = permute(volPsd.svdTrialGramMD.vec.coh(:,:,:,:,:,:,:,1),[5 7 2 8 1 3 4 6]);
    case 'MDpc'
        c = permute(volPsd.svdTrialGramMD.vec.cohEPC(:,:,:,:,:,:,:,1),[5 7 2 8 1 3 4 6]);
end
Fs  = volPsd.param.Fs;
switch mFlag
    case {'av' 'pc' 'ek'}
        f   = permute(volPsd.svdTrialGram.f              ,[5 7 2 8 1 3 4 6]);
        t   = permute(mean(volPsd.svdTrialGram.t(:,1,:,:,:,:,:),1),[5 7 2 8 1 3 4 6]);
        K   = volPsd.psdTrialGram.K;
        T   = volPsd.psdTrialGram.T;
        E   = volPsd.psdTrialGram.E;
        win = volPsd.psdTrialGram.win;
    case {'MD' 'MDpc'}
        f   = permute(volPsd.svdTrialGramMD.f              ,[5 7 2 8 1 3 4 6]);
        t   = permute(mean(volPsd.svdTrialGramMD.t(:,1,:,:,:,:,:),1),[5 7 2 8 1 3 4 6]);
        K   = volPsd.psdTrialGramMD.K;
        T   = volPsd.psdTrialGramMD.T;
        E   = volPsd.psdTrialGramMD.E;
        win = volPsd.psdTrialGramMD.win;
end
[TW,W] = K2W(T,K,0);

% coh = permute(abs(volPsd.svdTrialGram.coh(1,1,:,:,:)),[3 5 1 2 4 6 7 8]);
% f   = permute(volPsd.svdTrialGram.f                  ,[3 5 1 2 4 6 7 8]);
% t   = permute(volPsd.svdTrialGram.tWin               ,[3 5 1 2 4 6 7 8]);

imagesc(t,f,c)

xlabel('time (s)')
ylabel('freq (Hz)')
ylabel(colorbar,'coherence');



hold on
% W = permute(volPsd.svdTrialGram.W,[3 5 1 2 4 6 7 8]); if length(unique(W))==1; W = unique(W); end
% K = permute(volPsd.svdTrialGram.K,[3 5 1 2 4 6 7 8]); if length(unique(K))==1; K = unique(K); end
% T = permute(volPsd.svdTrialGram.T,[3 5 1 2 4 6 7 8]); if length(unique(T))==1; T = unique(T); end
% E = length(volPsd.psdTrialGram.param.onsetList);
% KE = K*E;

runDur = volPsd.nframes/volPsd.psd.param.Fs;
xlim([0 runDur])

switch mFlag
    case 'ek'
        paramStr = ['(KE=' num2str(K*E) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec; E=' num2str(E) 'trials)'];
    otherwise
        paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec; E=' num2str(E) 'trials)'];
end
switch mFlag
    case 'av'
        title(['evoked coherogram ' paramStr])
    case 'pc'
        title(['phase-locked evoked coherogram ' paramStr])
    case 'ek'
        title(['evoked cross-trial coherogram ' paramStr])
    case 'MD'
        title(['discontinuous-taper evoked coherogram ' paramStr])
    case 'MDpc'
        title(['discontinuous-taper phase-locked evoked coherogram ' paramStr])
end



cLim = c(f>0.01,:);
cLim = [1/(K*E) max(cLim(:))];
clim(cLim)



addWin([],volPsd.svdTrialGram)
addW([],volPsd.svdTrialGram)
addOnset([],volPsd.svdTrialGram.param.onsetList)



ax = [ax{:}];


%% Add to timeseries
axTs = findobj(allchild(F.Children),'type','axes');
ttl = get(axTs,'Title'); if ~iscell(ttl); ttl = {ttl}; end
ttl = get([ttl{:}],'String');
axTs = axTs(contains(ttl,'timeseries'));
if ~isempty(axTs)
    axes(axTs)

    %%% delete previous window size visual elements (magentat lines)
    hLine = findobj(axTs.Children,'type','Line'); %hLine = {hLine(:)};
    mLine = get(hLine,'Color'); if iscell(mLine); mLine = cat(1,mLine{:}); end
    mLine = all(cat(1,mLine)==[1 0 1],2);
    delete(hLine(mLine));

    for e = 1:E
        %%% add window size (first window)
        x = mean(volPsd.svdTrialGram.t(:,e,:,:,:,:,1,:),1);
        y = -inf;
        addWin([],volPsd.svdTrialGram,x,y)
        %%% add window size (last window)
        x = mean(volPsd.svdTrialGram.t(:,e,:,:,:,:,end,:),1);
        y = +inf;
        addWin([],volPsd.svdTrialGram,x,y)
    end
end


%% Add to coherogram
axSvdGram = findobj(allchild(F.Children),'type','axes');
ttl = get(axSvdGram,'Title'); if ~iscell(ttl); ttl = {ttl}; end
ttl = get([ttl{:}],'String');
axSvdGram = axSvdGram(contains(ttl,'coherogram') & ~contains(ttl,'trial-locked'));
if ~isempty(axSvdGram)
    axes(axSvdGram);
    %%% add trial onset
    addOnset([],volPsd.svdTrialGram.param.onsetList)
end
