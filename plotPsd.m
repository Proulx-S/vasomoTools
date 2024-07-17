function [ax,F] = plotPsd(volPsd,H)
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

%%% average spectrum
volPsd = vol2vec(volPsd);
f   = squeeze(volPsd.f);
psd = squeeze(mean(volPsd.vec,2));
plot(f,psd,'k')
grid on
axis tight
xlabel('f (Hz)')
ylabel('psd')
ax{end}.YScale = 'log';

W = volPsd.W;
K = volPsd.K;
T = volPsd.T;
paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec)'];
title(['spatially averaged spectrum ' paramStr])

yLim = psd(f>0.01);
yLim = [min(yLim) max(yLim)];
ylim(yLim)

addW([],volPsd.psd)



if isfield(volPsd,'psdTrialGram')
    fStim = 1/mean(diff(volPsd.psdTrialGram.param.onsetList));
    xline(fStim,'Color','r','linestyle','--')
end


ax = [ax{:}];




