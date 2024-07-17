function [ax,F] = plotPsdTrialGram(volPsd,H)
if ~exist('H','var'); H = [];                             end
if isempty(H);        H = figure('WindowStyle','docked'); end

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
psd = permute(mean(abs(volPsd.psdTrialGram.vec),3) ,[1 5 2 3 4 6 7 8]);
f   = permute(volPsd.psdTrialGram.f   ,[1 5 2 3 4 6 7 8]);
t   = permute(volPsd.psdTrialGram.tWin,[1 5 2 3 4 6 7 8]);

imagesc(t,f,psd)

xlabel('time (s)')
ylabel('freq (Hz)')
ylabel(colorbar,'psd');

ax{end}.ColorScale = 'log';

hold on
W = permute(volPsd.psdTrialGram.W,[1 5 2 3 4 6 7 8]); if length(unique(W))==1; W = unique(W); end
K = permute(volPsd.psdTrialGram.K,[1 5 2 3 4 6 7 8]); if length(unique(K))==1; K = unique(K); end
T = permute(volPsd.psdTrialGram.T,[1 5 2 3 4 6 7 8]); if length(unique(T))==1; T = unique(T); end
E = length(volPsd.psdTrialGram.param.onsetList);
winSz = volPsd.psdTrialGram.lWin;
winStep = volPsd.psdTrialGram.param.win(2);
if range(diff(f)) / mean(f) > 1e-15; error('variable dx'); end
df = mean(diff(f));

runDur = volPsd.nframes/volPsd.psdTrialGram.param.Fs;
xlim([0 runDur])

paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec; E=' num2str(E) 'trials)'];
title(['spatially trial-locked spectrogram ' paramStr])




cLim = psd(f>0.01,:);
cLim = [min(cLim(:)) max(cLim(:))];
clim(cLim)





addWin([],volPsd.psdTrialGram)
addW([],volPsd.psdTrialGram)
addOnset([],volPsd.psdTrialGram.param.onsetList)





ax = [ax{:}];



%% Add to time series plot
axTs = findobj(allchild(F.Children),'type','axes'); ttl = get(axTs,'Title'); ttl = get([ttl{:}],'String');
axTs = axTs(contains(ttl,'timeseries'));
axes(axTs);
yLim = ylim;

%%% delete previous window size visual elements (magentat lines)
if length(axTs.Children)>1
    mLine = get([axTs.Children(:)],'Color');
else
    mLine = {get([axTs.Children(:)],'Color')};
end
mLine = axTs.Children(all(cat(1,mLine{:})==[1 0 1],2));
delete(mLine);

firstWin  = mean(volPsd.t(volPsd.psdTrialGram.param.winInd(1  ,[1 end]  ,:)),2);
lastWin   = mean(volPsd.t(volPsd.psdTrialGram.param.winInd(end,[1 end]  ,:)),2);
for e = 1:length(firstWin)
    %%% add window size (first window)
    x = firstWin(e);
    y = -inf;
    addWin([],volPsd.psdTrialGram,x,y)
    %%% add window size (last window)
    x = lastWin(e);
    y = +inf;
    addWin([],volPsd.psdTrialGram,x,y)
end

%%% add trial onset
addOnset([],volPsd.psdTrialGram.param.onsetList)



%% Add to spectrogram
axPsdGram = findobj(allchild(F.Children),'type','axes');
ttl = get(axPsdGram,'Title'); ttl = get([ttl{:}],'String');
axPsdGram = axPsdGram(contains(ttl,'spectrogram') & ~contains(ttl,'trial-locked'));
axes(axPsdGram);

%%% add trial onset
addOnset([],volPsd.psdTrialGram.param.onsetList)


