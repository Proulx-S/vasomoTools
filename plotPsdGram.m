function [ax,F] = plotPsdGram(volPsd,H)
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
psd = permute(mean(abs(volPsd.psdGram.vec),3) ,[1 5 2 3 4 6 7 8]);
f   = permute(volPsd.psdGram.f   ,[1 5 2 3 4 6 7 8]);
t   = permute(volPsd.psdGram.tWin,[1 5 2 3 4 6 7 8]);
K   = volPsd.psdGram.K;
W   = unique(volPsd.psdGram.W);
T   = unique(volPsd.psdGram.T);

imagesc(t,f,psd)

xlabel('time (s)')
ylabel('freq (Hz)')
ylabel(colorbar,'psd');

ax{end}.ColorScale = 'log';

runDur = volPsd.nframes/volPsd.psd.param.Fs;
xlim([0 runDur])

paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec)'];
title(['spatially averaged spectrogram ' paramStr])


cLim = psd(f>0.01,:);
cLim = [min(cLim(:)) max(cLim(:))];
clim(cLim)




addWin([],volPsd.psdGram)
addW([],volPsd.psdGram)






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

%%% add window size
addWin(axTs,volPsd.psdGram)

