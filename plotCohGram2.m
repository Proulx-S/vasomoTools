function [ax,F] = plotCohGram2(volPsd,H)
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

if isfield(volPsd.svdGram,'vec')
    if isempty(volPsd.svdGram.vec)
        ax = []; F = []; return
    end
else
    if isempty(volPsd.svdGram.COH)
        ax = []; F = []; return
    end
end


%% spatial average
if isfield(volPsd.svdGram,'vec')
    c = permute(mean(volPsd.svdGram.vec(:,:,:,:,:,:,:,1),6),[5 7 2 8 1 3 4 6]);
else
    c = permute(mean(volPsd.svdGram.COH(:,:,:,:,:,:,:,1),6),[5 7 2 8 1 3 4 6]);
end
f   = permute(volPsd.svdGram.f                         ,[5 7 2 8 1 3 4 6]);
t   = permute(mean(volPsd.svdGram.t,1)                 ,[5 7 2 8 1 3 4 6]);
Fs  = volPsd.svdGram.param.Fs;
K   = volPsd.svdGram.K;
T   = volPsd.svdGram.T;
[TW,W] = K2W(T,K,0);


% coh = permute(abs(volPsd.svdGram.coh(1,1,:,:,:)),[3 5 1 2 4 6 7 8]);
% f   = permute(volPsd.svdGram.f                  ,[3 5 1 2 4 6 7 8]);
% t   = permute(volPsd.svdGram.tWin               ,[3 5 1 2 4 6 7 8]);

imagesc(t,f,c)

xlabel('time (s)')
ylabel('freq (Hz)')
ylabel(colorbar,'coherence');



hold on
% W = permute(volPsd.svdGram.W,[3 5 1 2 4 6 7 8]); if length(unique(W))==1; W = unique(W); end
% K = permute(volPsd.svdGram.K,[3 5 1 2 4 6 7 8]); if length(unique(K))==1; K = unique(K); end
% T = permute(volPsd.svdGram.T,[3 5 1 2 4 6 7 8]); if length(unique(T))==1; T = unique(T); end

runDur = volPsd.nframes/volPsd.psd.param.Fs;
xlim([0 runDur])


paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec)'];
title(['coherogram ' paramStr])



cLim = c(f>0.01,:);
cLim = [1/K max(cLim(:))];
clim(cLim)



addWin([],volPsd.svdGram)
addW([],volPsd.svdGram)





if isfield(volPsd,'svdTrialGram')
    onsetList = volPsd.svdTrialGram.onsetList;
    ondurList   = volPsd.svdTrialGram.durList;
    fStim = 1/mean(diff(onsetList));
    yline(fStim,'Color','r','linestyle',':')
    if length(unique(ondurList))==1
        yline(1/mean(ondurList(1)),'Color','g','linestyle',':')
        offsetList = onsetList + ondurList;
        offdurList = onsetList(2:end) - offsetList(1:end-1);
        yline(1/mean(offdurList),'Color','b','linestyle',':')
    end
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

%%% add window
addWin(axTs,volPsd.svdGram)


