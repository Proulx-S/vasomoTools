function [ax,F] = plotCoh2(volPsd,H)
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

if isfield(volPsd.svd,'vec')
    if isempty(volPsd.svd.vec)
        ax = []; F = []; return
    end
else
    if isempty(volPsd.svd.COH)
        ax = []; F = []; return
    end
end

%%% coherence spectrum
f = squeeze(volPsd.svd.f);
if isfield(volPsd.svd,'vec')
    c = squeeze(volPsd.svd.vec(:,:,:,:,:,:,:,1));
else
    c = squeeze(volPsd.svd.COH(:,:,:,:,:,:,:,1));
end
plot(f,c,'k')
grid on
axis tight
xlabel('f (Hz)')
ylabel('coherence')

K = volPsd.svd.K;
T = volPsd.svd.T;
[TW,W] = K2W(T,K,0);



paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec)'];
title(['coherence spectrum ' paramStr])




yLim = [1/K max(c(f>0.01))];
ylim(yLim);
xlim([0 f(end)])




addW([],volPsd.svd)



% if isfield(volPsd,'svdTrialGram')
%     fStim = 1/mean(diff(volPsd.svdTrialGram.param.onsetList));
%     xline(fStim,'Color','r','linestyle','--')
% end

if isfield(volPsd,'svdTrialGram')
    onsetList = volPsd.svdTrialGram.onsetList;
    ondurList   = volPsd.svdTrialGram.durList;
    fStim = 1/mean(diff(onsetList));
    xline(fStim,'Color','r','linestyle','--')
    if length(unique(ondurList))==1
        xline(1/mean(ondurList(1)),'Color','g','linestyle','--')
        offsetList = onsetList + ondurList;
        offdurList = onsetList(2:end) - offsetList(1:end-1);
        xline(1/mean(offdurList),'Color','b','linestyle','--')
    end
end





ax = [ax{:}];




