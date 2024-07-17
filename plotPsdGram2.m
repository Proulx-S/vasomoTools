function [ax,F] = plotPsdGram2(volPsd,H)
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
if isfield(volPsd.psdGram,'vec')
    psd = permute(mean(volPsd.psdGram.vec,6),[5 7 2 8 1 3 4 6]);
else
    psd = permute(mean(volPsd.psdGram.PSD,6),[5 7 2 8 1 3 4 6]);
end
f   = permute(volPsd.psdGram.f               ,[5 7 2 8 1 3 4 6]);
Fs  = volPsd.psdGram.param.Fs;
t   = permute(mean(volPsd.psdGram.t,1),[5 7 2 8 1 3 4 6]);
K   = volPsd.psdGram.K;
T   = volPsd.psdGram.T;
[TW,W] = K2W(T,K,0);

imagesc(t,f,psd)

xlabel('time (s)')
ylabel('freq (Hz)')
ylabel(colorbar,'psd');

ax{end}.ColorScale = 'log';


runDur = volPsd.nframes/volPsd.psd.param.Fs;
xlim([0 runDur])

paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec; TW=' num2str(TW) ')'];
title(['spectrogram ' paramStr])


cLim = psd(f>0.01,:);
cLim = [min(cLim(:)) max(cLim(:))];
clim(cLim)




addWin([],volPsd.psdGram)
addW([],volPsd.psdGram)




if isfield(volPsd,'psdTrialGram') && ~isempty(volPsd.psdTrialGram)
    addFreq([],volPsd.psdTrialGram.onsetList,volPsd.psdTrialGram.durList)
end






ax = [ax{:}];



%% Add to timeseries
axTs = findobj(allchild(F.Children),'type','axes'); ttl = get(axTs,'Title'); ttl = get([ttl{:}],'String');
axTs = axTs(contains(ttl,'timeseries'));

%%% delete previous window size visual elements (magentat lines)
hLine = findobj(axTs.Children,'type','Line'); %hLine = {hLine(:)};
mLine = get(hLine,'Color'); if iscell(mLine); mLine = cat(1,mLine{:}); end
mLine = all(cat(1,mLine)==[1 0 1],2);
delete(hLine(mLine));

% % if length(hLine)==1
% %     hLine = {hLine};
% % end
% % 
% %     mLine = get(hLine,'Color');
% % else
% %     mLine = {get([axTs.Children(:)],'Color')};
% % end
% mLine = get(hLine,'Color'); if ~iscell(mLine); mLine = {mLine}; end
% ind = all(cat(1,mLine{:})==[1 0 1],2);
% mLine(ind)
% 
% mLine = hLine(all(cat(1,mLine{:})==[1 0 1],2));
% delete(mLine);
% % if length(axTs.Children)>1
% %     mLine = get(hLine,'Color');
% % else
% %     mLine = {get([axTs.Children(:)],'Color')};
% % end
% % mLine = axTs.Children(all(cat(1,mLine{:})==[1 0 1],2));
% % delete(mLine);

%%% add window size
addWin(axTs,volPsd.psdGram)

