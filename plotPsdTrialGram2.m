function [ax,F] = plotPsdTrialGram2(volPsd,H,mFlag)
if ~exist('H','var');         H = []; end
if ~exist('mFlag','var'); mFlag = []; end
if isempty(H);                H = figure('WindowStyle','docked'); end
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


%% 
switch mFlag
    case 'av'
        psd = permute(mean(volPsd.psdTrialGram.vec.psd,6),[5 7 2 8 1 3 4 6]);
    case 'pc'
        psd = permute(mean(volPsd.psdTrialGram.vec.psdPC,6),[5 7 2 8 1 3 4 6]);
    case 'MD'
        psd = permute(mean(volPsd.psdTrialGramMD.vec.psd,6),[5 7 2 8 1 3 4 6]);
    case 'MDpc'
        psd = permute(mean(volPsd.psdTrialGramMD.vec.psdPC,6),[5 7 2 8 1 3 4 6]);
end
switch mFlag
    case {'av' 'pc'}
        f = permute(volPsd.psdTrialGram.f              ,[5 7 2 8 1 3 4 6]);
        t = permute(mean(volPsd.psdTrialGram.t(:,1,:,:,:,:,:),1),[5 7 2 8 1 3 4 6]);
        K   = volPsd.psdTrialGram.K;
        T   = volPsd.psdTrialGram.T;
        E   = volPsd.psdTrialGram.E;
        win = volPsd.psdTrialGram.win;
    case {'MD' 'MDpc'}
        f = permute(volPsd.psdTrialGramMD.f              ,[5 7 2 8 1 3 4 6]);
        t = permute(mean(volPsd.psdTrialGramMD.t(:,1,:,:,:,:,:),1),[5 7 2 8 1 3 4 6]);
        K   = volPsd.psdTrialGramMD.K;
        T   = volPsd.psdTrialGramMD.T;
        E   = volPsd.psdTrialGramMD.E;
        win = volPsd.psdTrialGramMD.win;
end
Fs  = volPsd.param.Fs;
[TW,W] = K2W(T,K,0);


% %% spatial average
% psd = permute(mean(abs(volPsd.psdTrialGram.vec),3) ,[1 5 2 3 4 6 7 8]);
% f   = permute(volPsd.psdTrialGram.f   ,[1 5 2 3 4 6 7 8]);
% t   = permute(volPsd.psdTrialGram.tWin,[1 5 2 3 4 6 7 8]);

imagesc(t,f,psd)

xlabel('time (s)')
ylabel('freq (Hz)')
ylabel(colorbar,'psd');

ax{end}.ColorScale = 'log';

hold on
% W = permute(volPsd.psdTrialGram.W,[1 5 2 3 4 6 7 8]); if length(unique(W))==1; W = unique(W); end
% K = permute(volPsd.psdTrialGram.K,[1 5 2 3 4 6 7 8]); if length(unique(K))==1; K = unique(K); end
% T = permute(volPsd.psdTrialGram.T,[1 5 2 3 4 6 7 8]); if length(unique(T))==1; T = unique(T); end
% E = length(volPsd.psdTrialGram.param.onsetList);
% winSz = volPsd.psdTrialGram.lWin;
% winStep = volPsd.psdTrialGram.param.win(2);
% if range(diff(f)) / mean(f) > 1e-15; error('variable dx'); end
% df = mean(diff(f));

runDur = volPsd.nframes/volPsd.param.Fs;
xlim([0 runDur])

paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec TW=' num2str(TW) '; E=' num2str(E) 'trials)'];
switch mFlag
    case 'av'
        title(['evoked spectrogram ' paramStr])
    case 'pc'
        title(['phase-locked evoked spectrogram ' paramStr])
    case 'MD'
        title(['discontinuous-taper evoked spectrogram ' paramStr])
    case 'MDpc'
        title(['discontinuous-taper phase-locked evoked spectrogram ' paramStr])
end



cLim = psd(f>0.01,:);
cLim = [min(cLim(:)) max(cLim(:))];
clim(cLim)





addWin([],volPsd.psdTrialGram)
addW([],volPsd.psdTrialGram)
addOnset([],volPsd.psdTrialGram.onsetList)
addFreq([],volPsd.psdTrialGram.onsetList,volPsd.psdTrialGram.durList)






ax = [ax{:}];



%% Add to time series plot
axTs = findobj(allchild(F.Children),'type','axes'); ttl = get(axTs,'Title'); if ~iscell(ttl); ttl = {ttl}; end; ttl = get([ttl{:}],'String');
axTs = axTs(contains(ttl,'timeseries'));
if ~isempty(axTs)
    axes(axTs);
    yLim = ylim;

    %%% delete previous window size visual elements (magentat lines)
    hLine = findobj(axTs.Children,'type','Line'); %hLine = {hLine(:)};
    mLine = get(hLine,'Color'); if iscell(mLine); mLine = cat(1,mLine{:}); end
    mLine = all(cat(1,mLine)==[1 0 1],2);
    delete(hLine(mLine));

    for e = 1:E
        %%% add window size (first window)
        x = mean(volPsd.psdTrialGram.t(:,e,:,:,:,:,1,:),1);
        y = -inf;
        addWin([],volPsd.psdTrialGram,x,y)
        %%% add window size (last window)
        x = mean(volPsd.psdTrialGram.t(:,e,:,:,:,:,end,:),1);
        y = +inf;
        addWin([],volPsd.psdTrialGram,x,y)
    end
end

% %%% add trial onset
% addOnset([],volPsd.psdTrialGram.param.onsetList)



% %% Add to spectrogram
% axPsdGram = findobj(allchild(F.Children),'type','axes');
% ttl = get(axPsdGram,'Title'); if ~iscell(ttl); ttl = {ttl}; end
% ttl = get([ttl{:}],'String');
% axPsdGram = axPsdGram(contains(ttl,'spectrogram') & ~contains(ttl,'evoked'));
% if ~isempty(axPsdGram)
%     axes(axPsdGram);
% 
%     %%% add trial onset
%     addOnset([],volPsd.psdTrialGram.param.onsetList)
% 
% 
%     %%% add frequencies
%     onsetList = volPsd.psdTrialGram.onsetList;
%     ondurList   = volPsd.psdTrialGram.durList;
%     fStim = 1/mean(diff(onsetList));
%     yline(fStim,'Color','r','linestyle',':')
%     xLim = xlim;
%     text(xLim(2),fStim,'fStim','HorizontalAlignment','right','VerticalAlignment','baseline','Color','r')
%     if length(unique(ondurList))==1
%         fOn = 1/mean(ondurList(1));
%         yline(fOn,'Color','g','linestyle',':')
%         text(xLim(2),fOn,'fOn','HorizontalAlignment','right','VerticalAlignment','top','Color','g')
%         offsetList = onsetList + ondurList;
%         offdurList = onsetList(2:end) - offsetList(1:end-1);
%         fOff = 1/mean(offdurList);
%         yline(fOff,'Color','b','linestyle',':')
%         text(xLim(2),fOff,'fOff','HorizontalAlignment','right','VerticalAlignment','top','Color','b')
%     end
% end

