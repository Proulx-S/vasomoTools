function axIm = viewMTsvdKlein2(svdStruct,funPsd,peakFreq,imType,axIm)

if ~exist("imType",'var') || isempty(imType)
    imType = 'psd'; % svMag, svPhase or psd
end

tmp = strsplit(funPsd.fspec,filesep); tmp = strsplit(tmp{end-1},'_');
sub = replace(tmp{contains(tmp,'sub-')},'sub-','');
ses = replace(tmp{contains(tmp,'ses-')},'ses-','');
run = replace(tmp{contains(tmp,'run-')},'run-','');
fpass = svdStruct.param.fpass;
fLim = [0 2.5];
K = svdStruct.param.tapers(2);
f = svdStruct.f;

%% Plot timeseries mean (underlay) only of not already plotted in axIm
if ~exist('axIm','var') || isempty(axIm)
    axIm = {};
    tmp = []; tmp.String = '';
else
    tmp = [axIm{:}]; tmp = [tmp.Title];
end
if ~any(contains({tmp.String},'timeseries mean'))
    figure('WindowStyle','docked');
    hTile = tiledlayout(6,1); hTile.Padding = 'tight'; hTile.TileSpacing = "tight";
    nexttile([5 1]);
    imagesc(funPsd.tMean); hold on
    maskC = getMaskOutline(funPsd.roi.mask,5);
    hMask = plot(maskC); hMask.FaceColor = 'none'; hMask.EdgeColor = 'r';
    axIm{end+1} = gca;
    %%% Crop image
    cropPrc = 75;
    xLim = [find(any(funPsd.tMean>prctile(funPsd.tMean(:),cropPrc),1),1,'first') find(any(funPsd.tMean>prctile(funPsd.tMean(:),cropPrc),1),1,'last')];
    yLim = [find(any(funPsd.tMean>prctile(funPsd.tMean(:),cropPrc),2),1,'first') find(any(funPsd.tMean>prctile(funPsd.tMean(:),cropPrc),2),1,'last')];
    w = max([diff(xLim) diff(yLim)]);
    xLim = mean(xLim)+[-1 1].*w/2;
    yLim = mean(yLim)+[-1 1].*w/2;
    axis(axIm{end},[xLim yLim]);
    axIm{end}.PlotBoxAspectRatio = [1 1 1]; axIm{end}.DataAspectRatio = [1 1 1];
    axIm{end}.XTick = []; axIm{end}.YTick = [];
    axIm{end}.Colormap = gray;
    cb = colorbar;
    ylabel(cb,'MR intensity')
    title(['sub-' sub '; ses-' ses '; run-' run ' K=' num2str(K) '; timeseries mean'])

    %%% Plot spectrum
    nexttile([1 1]);
    yyaxis left
    plot(funPsd.roi.f,funPsd.roi.psd); hold on
    axis tight
    ylabel('psd')
    yyaxis right
    plot(f,svdStruct.c); hold on
    ylim([1/K 1])
    ylabel('coherence')
    ax1 = gca;
    ax1.YAxis(2).Scale = 'linear';
    ax1.YAxis(1).Scale = 'log';
    xlim(fLim)
    grid minor

    %%% Plot spectrum peaks
    if ~isnan(peakFreq)
        plot(([1 1].*peakFreq')',repmat(ylim,length(peakFreq),1)','-k')
        yLim = ylim;
        plot((peakFreq'+[-1 1].*svdStruct.w)',yLim(1).*ones(size(peakFreq,2),2)','-g','LineWidth',3)
        xlabel('Hz')
    end

    %%% Adjust
    cb.Location = 'manual';
    drawnow
    cb.Position = [sum(axIm{end}.Position([1 3])) axIm{end}.Position(2) 0.02 axIm{end}.Position(4)];
    cb.Visible = 'off';
end

peakFreq = peakFreq(~isnan(peakFreq));
%% Plot each frequency
for fInd = 1:length(peakFreq)
    [~,b] = min(abs(f - peakFreq(fInd)));
    figure('WindowStyle','docked');
    %%% Plot map
    axIm{end+1} = axes(gcf,'Units',axIm{end}.Units,'Position',axIm{end}.Position);
    switch imType
        case 'psd'
            curIm = vec2vol(funPsd); curIm = curIm.vol(:,:,:,b);
            hIm = imagesc(curIm);
            axIm{end}.Colormap = jet;
            cb = colorbar;
            ylabel(cb,'psd')
            axIm{end}.ColorScale = 'log';
            cLim = curIm(logical(funPsd.roi.mask)); cLim = [1 max(cLim)];
            axIm{end}.CLim = cLim;
        case 'svMag'
            curIm = nan(size(svdStruct.mask));
            curIm(logical(svdStruct.mask)) = svdStruct.sp(:,b,1);
            hIm = imagesc(abs(curIm));
            axIm{end}.Colormap = jet;
            cb = colorbar;
            ylabel(cb,'spatial sv weigth mag')
            cLim = abs(curIm(logical(funPsd.roi.mask))); cLim = [min(cLim) max(cLim)];
            axIm{end}.CLim = cLim;
            cbTicks = cb.Ticks; if cbTicks(1) ~= cLim(1); cbTicks = [cLim(1) cbTicks]; end; if cbTicks(end) ~= cLim(2); cbTicks = [cbTicks cLim(2)]; end; cb.Ticks = cbTicks; cb.TickLabels(2:end-1) = {''};
            cb.TickLabels{1} = 'min'; cb.TickLabels{end} = 'max';
        case 'svPhase'
            curIm = nan(size(svdStruct.mask));
            curIm(logical(svdStruct.mask)) = svdStruct.sp(:,b,1);
            hIm = imagesc(wrapToPi(angle(curIm) - angle(mean(curIm(:),'omitnan'))));
            hIm.AlphaData = ( abs(curIm) - min(abs(curIm(:))) ) ./ max(abs(curIm(:)));
            axIm{end}.Color = [1 1 1].*0.5;
            axIm{end}.Colormap = hsv;
            cb = colorbar;
            hYl = ylabel(cb,'spatial sv weigth phase');
            cLim = [-pi pi];
            axIm{end}.CLim = cLim;
            cb.Ticks = [-pi -pi/2 0 pi/2 pi];
            cb.TickLabels = {'-pi' '-pi/2' '0' 'pi/2' 'pi'};     
    end
    axIm{end}.PlotBoxAspectRatio = [1 1 1]; axIm{end}.DataAspectRatio = [1 1 1];
    axIm{end}.XTick = []; axIm{end}.YTick = [];
    title(['sub-' sub '; ses-' ses '; run-' run ' K=' num2str(K) '; ' num2str(peakFreq(fInd),'%.3f') 'Hz'])
    axis([axIm{end-1}.XLim axIm{end-1}.YLim])

    %% Plot spectra
    axes(gcf,'Units',axIm{end-1}.Parent.Children(1).Units,'Position',axIm{end-1}.Parent.Children(1).Position);
    yyaxis left
    plot(funPsd.roi.f,funPsd.roi.psd); hold on
    axis tight
    ylabel('psd')
    ytickformat('%2.0f')
    yyaxis right
    plot(f,svdStruct.c); hold on
    ylim([1/K 1])
    ylabel('coherence')
    ax = gca; ax.YAxis(2).Scale = 'linear'; ax.YAxis(1).Scale = 'log';
    xlim(fLim)
    grid minor

    plot(([1 1].*f(b)')',ylim,'-k')
    yLim = ylim;
    plot((f(b)'+[-1 1].*svdStruct.w)',yLim(1).*[1 1]','-g','LineWidth',3)
    xlabel('Hz')

    %% Adjust
%     cb.Location = 'manual';
%     drawnow
    axIm{end}.Position = axIm{end-1}.Position;
%     cb.Position = [sum(axIm{end}.Position([1 3])) axIm{end}.Position(2) 0.02 axIm{end}.Position(4)];
end
linkaxes([axIm{:}])


