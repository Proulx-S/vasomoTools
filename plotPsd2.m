function [ax,F] = plotPsd2(volPsd,H,tWin)
if ~exist('H','var'); H = [];                             end
if isempty(H);        H = figure('WindowStyle','docked'); end
if ~exist('tWin','var'); tWin = [];                             end

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

%%% average across space
% volPsd = vol2vec(volPsd);
if isfield(volPsd.psd,'vec')
    psd = squeeze(mean(volPsd.psd.vec,6));
else
    psd = squeeze(mean(volPsd.psd.PSD,6));
end
f = squeeze(volPsd.psd.f);
% f   = squeeze(volPsd.f);
% psd = squeeze(mean(volPsd.vec,2));
plot(f,psd,'k')
grid on
axis tight
xlabel('f (Hz)')
ylabel('psd')
ax{end}.YScale = 'log';

K = volPsd.psd.K;
T = volPsd.psd.T;
[TW,W] = K2W(T,K,0);

paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec); TW=' num2str(TW)];
title(['spatially averaged spectrum ' paramStr])

yLim = psd(f>0.01);
yLim = [min(yLim) max(yLim)];
ylim(yLim)
xlim([0 f(end)])

addW([],volPsd.psd)



if isfield(volPsd,'psdTrialGram') && ~isempty(volPsd.psdTrialGram)
    addFreq([],volPsd.psdTrialGram.onsetList,volPsd.psdTrialGram.durList)
end

if isfield(volPsd,'psdTrialGramMD') && ~isempty(volPsd.psdTrialGramMD) && ~isempty(tWin)
    t   = volPsd.psdTrialGramMD.t(1,1,1,1,1,1,:,1);
    f   = volPsd.psdTrialGramMD.f(1,1,1,1,:,1,1,1);
    if tWin==inf
        psd = mean(mean(volPsd.psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,:,:),6),7);
    else
        [~,wInd] = min(abs(t-tWin));
        psd = mean(volPsd.psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,wInd,:),6);
    end
    plot(squeeze(f),squeeze(psd),'--k')
end


ax = [ax{:}];




