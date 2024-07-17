function hF = viewFullMT_cohgram(funPsd)


xLim = [0 3800];
fLim = [0 0.4];

hF = figure('WindowStyle','docked');
modeN = 1;
t = squeeze(funPsd.svdGram.tWin);
f = squeeze(funPsd.svdGram.f);
gram = squeeze(funPsd.svdGram.coh(:,modeN,:,:,:));
for runInd = 1:funPsd.nruns
    hSurf = surf(t(runInd,:),f,squeeze(gram(:,runInd,:)),squeeze(gram(:,runInd,:)),'EdgeColor','none','FaceColor','interp');
    hold on
end
ax = gca;
ax.View = [0 90];
ax.Color = 'none';
ax.Colormap = jet;
ax.ColorScale = 'log';
grid off
xlabel('time (s)')
ylabel('Hz')
ylabel(colorbar,'coherence')
ylim(fLim)
ax.CLim(1) = 1/funPsd.svdGram.K;
title(['Coherogram; W=' num2str(funPsd.svdGram.W(1)) 'Hz; K=' num2str(funPsd.svdGram.K) '; winSz=' num2str(funPsd.svdGram.lWin) 's'])
xlim(xLim);
x = mean([funPsd.t(1,1,1,end,1,1) funPsd.t(1,1,1,1,1,2)]);
y = mean(ylim);
z = ax.CLim;
x = x + 0.5.*[-1 1;0 0].*funPsd.svdGram.lWin;
y = y + 1.*[0 0;-1 1].*funPsd.svdGram.W(1);
% z = z(2) + 0.1.*[1 1;1 1].*diff(z);
z = z(2) + 0.1.*diff(z);
hPlot = plot(x',y','g','LineWidth',2); uistack(hPlot,'top');
hPlot(1).ZData = ones(size(hPlot(1).XData)).*z; hPlot(2).ZData = ones(size(hPlot(2).XData)).*z;
% hPlot = plot3(x',y',z','g','LineWidth',2); uistack(hPlot,'top');
y = squeeze(funPsd.f(:,:,:,:,:,[1 end]));
for runInd = 1:funPsd.nruns
    x = squeeze(funPsd.t(:,:,:,[1 end],:,runInd));
    patch(x([1 2 2 1 1]),y([1 1 2 2 1]),[0.5 0.5 0.5],'EdgeColor','none');
end

% nexttile
% modeN = 1;
% t = squeeze(funPsd.psdGram.tWin);
% f = squeeze(funPsd.psdGram.f);
% gram = squeeze(mean(funPsd.psdGram.vec(:,:,:,:,:,modeN),3));
% for runInd = 1:funPsd.nruns
%     % hSurf = surf(t(runInd,:),f,squeeze(gram(:,runInd,:)),'EdgeColor','w','EdgeAlpha',0.2,'FaceColor','interp');
%     hSurf = surf(t(runInd,:),f,squeeze(gram(:,runInd,:)),'EdgeColor','none','FaceColor','interp');
%     hold on
% end
% ax = gca;
% ax.View = [0 90];
% ax.Color = 'none';
% ax.ColorScale = 'log';
% ax.Colormap = jet;
% grid off
% xlabel('time (s)')
% ylabel('Hz')
% ylabel(colorbar,'psd')
% ylim([0 0.5])
% if isfield(funPsd,'scale')
%     ax.CLim(1) = 1;
%     title(['Spectrogram; W=' num2str(funPsd.psdGram.W(1)) 'Hz; K=' num2str(funPsd.psdGram.K) '; winSz=' num2str(funPsd.psdGram.lWin) 's; normalized'])
% else
%     title(['Spectrogram; W=' num2str(funPsd.psdGram.W(1)) 'Hz; K=' num2str(funPsd.psdGram.K) '; winSz=' num2str(funPsd.psdGram.lWin) 's'])
% end
% sub = subList{I};
% title(hTile,['sub-' sub]);
% xlim(xLim);
% x = mean([funPsd.t(1,1,1,end,1,1) funPsd.t(1,1,1,1,1,2)]);
% y = mean(ylim);
% z = ax.CLim;
% x = x + 0.5.*[-1 1;0 0].*funPsd.psdGram.lWin;
% y = y + 1.*[0 0;-1 1].*funPsd.psdGram.W(1);
% % z = z(2) + 0.1.*[1 1;1 1].*diff(z);
% z = z(2) + 0.1.*diff(z);
% hPlot = plot(x',y','g','LineWidth',2); uistack(hPlot,'top');
% hPlot(1).ZData = ones(size(hPlot(1).XData)).*z; hPlot(2).ZData = ones(size(hPlot(2).XData)).*z;
% % hPlot = plot3(x',y',z','g','LineWidth',2); uistack(hPlot,'top');
% y = squeeze(funPsd.f(:,:,:,:,:,[1 end]));
% for runInd = 1:funPsd.nruns
%     x = squeeze(funPsd.t(:,:,:,[1 end],:,runInd));
%     patch(x([1 2 2 1 1]),y([1 1 2 2 1]),[0.5 0.5 0.5],'EdgeColor','none')
% end
% hFig = gcf;
% axGram = findobj(hFig.Children.Children,'Type','Axes');
