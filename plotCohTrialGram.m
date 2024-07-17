function [ax,F] = plotCohTrialGram(volPsd,H)
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
coh = permute(abs(volPsd.svdTrialGram.coh(1,1,:,:,:)),[3 5 1 2 4 6 7 8]);
f   = permute(volPsd.svdTrialGram.f                  ,[3 5 1 2 4 6 7 8]);
t   = permute(volPsd.svdTrialGram.tWin               ,[3 5 1 2 4 6 7 8]);

imagesc(t,f,coh)

xlabel('time (s)')
ylabel('freq (Hz)')
ylabel(colorbar,'coherence');



hold on
W = permute(volPsd.svdTrialGram.W,[3 5 1 2 4 6 7 8]); if length(unique(W))==1; W = unique(W); end
K = permute(volPsd.svdTrialGram.K,[3 5 1 2 4 6 7 8]); if length(unique(K))==1; K = unique(K); end
T = permute(volPsd.svdTrialGram.T,[3 5 1 2 4 6 7 8]); if length(unique(T))==1; T = unique(T); end
E = length(volPsd.psdTrialGram.param.onsetList);
KE = K*E;

runDur = volPsd.nframes/volPsd.psd.param.Fs;
xlim([0 runDur])

paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec; E=' num2str(E) 'trials)'];
title(['trial-locked coherogram ' paramStr])




cLim = coh(f>0.01,:);
cLim = [1/(K*E) max(cLim(:))];
clim(cLim)



addWin([],volPsd.svdTrialGram)
addW([],volPsd.svdTrialGram)
addOnset([],volPsd.svdTrialGram.param.onsetList)



ax = [ax{:}];


%% Add to timeseries
axTs = findobj(allchild(F.Children),'type','axes'); ttl = get(axTs,'Title'); ttl = get([ttl{:}],'String');
axTs = axTs(contains(ttl,'timeseries'));
axes(axTs)

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



%% Add to coherogram
axPsdGram = findobj(allchild(F.Children),'type','axes');
ttl = get(axPsdGram,'Title'); ttl = get([ttl{:}],'String');
axPsdGram = axPsdGram(contains(ttl,'coherogram') & ~contains(ttl,'trial-locked'));
axes(axPsdGram);

%%% add trial onset
addOnset([],volPsd.svdTrialGram.param.onsetList)

