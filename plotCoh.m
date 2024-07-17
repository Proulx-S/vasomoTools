function [ax,F] = plotCoh(volPsd,H)
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

%%% coherence spectrum
f = squeeze(volPsd.svd.f);
c = squeeze(volPsd.svd.coh(1,1,:));
plot(f,c,'k')
grid on
axis tight
xlabel('f (Hz)')
ylabel('coherence')

K = volPsd.svd.K;
W = volPsd.svd.W;
T = volPsd.svd.T;



paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec)'];
title(['coherence spectrum ' paramStr])




yLim = [1/K max(c(f>0.01))];
ylim(yLim);





addW([],volPsd.svd)



if isfield(volPsd,'svdTrialGram')
    fStim = 1/mean(diff(volPsd.svdTrialGram.param.onsetList));
    xline(fStim,'Color','r','linestyle','--')
end




ax = [ax{:}];




