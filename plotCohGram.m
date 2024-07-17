function [ax,F] = plotCohGram(volPsd,H)
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
coh = permute(abs(volPsd.svdGram.coh(1,1,:,:,:)),[3 5 1 2 4 6 7 8]);
f   = permute(volPsd.svdGram.f                  ,[3 5 1 2 4 6 7 8]);
t   = permute(volPsd.svdGram.tWin               ,[3 5 1 2 4 6 7 8]);

imagesc(t,f,coh)

xlabel('time (s)')
ylabel('freq (Hz)')
ylabel(colorbar,'coherence');



hold on
W = permute(volPsd.svdGram.W,[3 5 1 2 4 6 7 8]); if length(unique(W))==1; W = unique(W); end
K = permute(volPsd.svdGram.K,[3 5 1 2 4 6 7 8]); if length(unique(K))==1; K = unique(K); end
T = permute(volPsd.svdGram.T,[3 5 1 2 4 6 7 8]); if length(unique(T))==1; T = unique(T); end

runDur = volPsd.nframes/volPsd.psd.param.Fs;
xlim([0 runDur])


paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec)'];
title(['coherogram ' paramStr])



cLim = coh(f>0.01,:);
cLim = [1/K max(cLim(:))];
clim(cLim)



addWin([],volPsd.svdGram)
addW([],volPsd.svdGram)




ax = [ax{:}];


%% Add to timeseries
axTs = findobj(allchild(F.Children),'type','axes'); ttl = get(axTs,'Title'); ttl = get([ttl{:}],'String');
axTs = axTs(contains(ttl,'timeseries'));

%%% delete previous window size visual elements (magentat lines)
if length(axTs.Children)>1
    mLine = get([axTs.Children(:)],'Color');
else
    mLine = {get([axTs.Children(:)],'Color')};
end
mLine = axTs.Children(all(cat(1,mLine{:})==[1 0 1],2));
delete(mLine);

%%% add window
addWin(axTs,volPsd.svdGram)


