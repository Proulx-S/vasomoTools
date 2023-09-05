function axIm = viewMTsvdKlein(svdStruct,funPsd,peakFreq,imType,axIm)
if ~exist("imType",'var') || isempty(imType)
    imType = 'sv'; % sv or psd
end

tmp = strsplit(funPsd.fspec,filesep); tmp = strsplit(tmp{end-1},'_');
sub = replace(tmp{contains(tmp,'sub-')},'sub-','');
ses = replace(tmp{contains(tmp,'ses-')},'ses-','');
run = replace(tmp{contains(tmp,'run-')},'run-','');
fpass = svdStruct.param.fpass;
K = svdStruct.param.tapers(2);

%%% plot power and eigen spectra
figure('WindowStyle','docked');
hTile = tiledlayout(3,3);
hTile.Padding = 'tight'; hTile.TileSpacing = 'tight';
nexttile([1 3])
yyaxis left
plot(funPsd.roi.f,funPsd.roi.psd); hold on
axis tight
ylabel('psd')
yyaxis right
plot(svdStruct.f,svdStruct.c); hold on
ylim([1/K 1])
ylabel('coherence')
ax1 = gca;
ax1.YAxis(2).Scale = 'linear';
ax1.YAxis(1).Scale = 'log';
xlim(fpass)
grid minor
title(['sub-' sub '; ses-' ses '; run-' run ' K=' num2str(K)])
drawnow

% figure('WindowStyle','docked');
% hTile = tiledlayout(3,3);
% hTile.Padding = 'tight'; hTile.TileSpacing = 'tight';
% nexttile([1 3])
% yyaxis left
% plot(funPsd.roi.f,funPsd.roi.psd); hold on
% axis tight
% ylabel('psd')
% ax = gca;
% ax.YScale = 'log';
% xlim(fpass)
% grid minor
% title(['sub-' sub '; ses-' ses '; run-' run ' K=' num2str(K)])
% drawnow
% yyaxis right
% ylabel(' ')

plot(([1 1].*peakFreq')',repmat(ylim,length(peakFreq),1)','-k')
yLim = ylim;
plot((peakFreq'+[-1 1].*svdStruct.w)',yLim(1).*ones(size(peakFreq,2),2)','-g','LineWidth',3)
text(peakFreq+0.005,ones(size(peakFreq)).*yLim(2),cellstr(num2str(peakFreq','%.3fHz')),'Rotation',-90,'VerticalAlignment','baseline')
xlabel('Hz')


if ~exist('axIm','var') || isempty(axIm)
    axIm = {};
end
axIm{end+1} = nexttile;
imagesc(funPsd.tMean); hold on
maskC = getMaskOutline(funPsd.roi.mask,5);
hMask = plot(maskC); hMask.FaceColor = 'none'; hMask.EdgeColor = 'r';
axIm{end}.PlotBoxAspectRatio = [1 1 1]; axIm{end}.DataAspectRatio = [1 1 1]; axIm{end}.XAxis.Visible = 'off'; axIm{end}.YAxis.Visible = 'off'; axIm{end}.Colormap = gray;
title('timeseries mean')
switch imType
    case 'psd'
        %%% power spectra maps
        for fInd = 1:length(peakFreq)
            if fInd>5; break; end
            axIm{end+1} = nexttile;
            [~,b] = min(abs(svdStruct.f - peakFreq(fInd)));
            pwrSpace = vec2vol(funPsd);
            cLim = permute(pwrSpace.vol,[4 3 1 2]);
            cLim = cLim(b,1,logical(svdStruct.mask));
            cLim = [1 max(cLim(:))];
            hIm = imagesc(pwrSpace.vol(:,:,:,b),cLim);
            axIm{end}.ColorScale = 'log';
            axIm{end}.PlotBoxAspectRatio = [1 1 1]; axIm{end}.DataAspectRatio = [1 1 1]; axIm{end}.XAxis.Visible = 'off'; axIm{end}.YAxis.Visible = 'off'; axIm{end}.Colormap = jet;
            title([num2str(peakFreq(fInd),'%.3f') 'Hz'])
        end
        linkaxes([axIm{:}])
        cb = colorbar;
        cbTicks = cb.Ticks; if cbTicks(1) ~= cLim(1); cbTicks = [cLim(1) cbTicks]; end; if cbTicks(end) ~= cLim(2); cbTicks = [cbTicks cLim(2)]; end; cb.Ticks = cbTicks; cb.TickLabels(2:end-1) = {''};
        cb.TickLabelInterpreter = 'latex';
        cb.TickLabels{1} = ['$$\begin{array}{c}' 'noise' '\\' 'floor' '\\' '\end{array}$$'];
        cb.TickLabels{end} = 'max';
        hYlabel = ylabel(cb,'psd');
        hYlabel.Units = 'normalized';
        hYlabel.Position(1) = 1;
    case 'sv'
        %%% spatial sv
        for fInd = 1:length(peakFreq)
            if fInd>5; break; end
            axIm{end+1} = nexttile;
            [~,b] = min(abs(svdStruct.f - peakFreq(fInd)));
            svSpace = nan(size(svdStruct.mask));
            svSpace(logical(svdStruct.mask)) = svdStruct.sp(:,b,1);
            imagesc(abs(svSpace));
            axIm{end}.PlotBoxAspectRatio = [1 1 1]; axIm{end}.DataAspectRatio = [1 1 1]; axIm{end}.XAxis.Visible = 'off'; axIm{end}.YAxis.Visible = 'off'; axIm{end}.Colormap = jet;
            title([num2str(peakFreq(fInd),'%.3f') 'Hz'])
            cLim = abs(svSpace(logical(svdStruct.mask)));
            cLim = [min(cLim(:)) max(cLim(:))];
            axIm{end}.CLim = cLim;
        end
        cb = colorbar;
        cbTicks = cb.Ticks; if cbTicks(1) ~= cLim(1); cbTicks = [cLim(1) cbTicks]; end; if cbTicks(end) ~= cLim(2); cbTicks = [cbTicks cLim(2)]; end; cb.Ticks = cbTicks; cb.TickLabels(2:end-1) = {''};
        cb.TickLabels{1} = 'min'; cb.TickLabels{end} = 'max';
        hYlabel = ylabel(cb,'mag of sv weigths');
        hYlabel.Units = 'normalized';
        hYlabel.Position(1) = 1;
end

%% Adjust some stuff
cb.Location = 'manual';
drawnow
cb.Position(1) = sum(axIm{end}.Position([1 3]));
linkaxes([axIm{:}])
