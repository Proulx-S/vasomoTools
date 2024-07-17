function F = viewMTsvdKlein3(svdStruct,funPsd,peakFreq,imType,fLim,axIm)

if ~exist("imType",'var') || isempty(imType)
    imType = 'psd'; % svMag, svPhase or psd
end

if exist('funPsd','var') && ~isempty(funPsd(1))
    tmp = strsplit(funPsd(1).fspec,filesep); tmp = strsplit(tmp{end-1},'_');
    sub = replace(tmp{contains(tmp,'sub-')},'sub-','');
    ses = replace(tmp{contains(tmp,'ses-')},'ses-','');
    run = replace(tmp{contains(tmp,'run-')},'run-','');
else
    sub = '???';
    ses = '???';
    run = '???';
end
if ~exist('fLim','var') || isempty(fLim)
    fLim = svdStruct.param.fpass;
end
K = svdStruct.param.tapers(2);
if isfield(svdStruct.param,'augFac')
    augFac = svdStruct.param.augFac;
else
    augFac = 1;
end
f = svdStruct.f;


%% Plot coherence spectrum
F = {figure('WindowStyle','docked')};
hP1 = plot(f,svdStruct.c); hold on
ylim([1/(K*augFac) 1])
ylabel('coherence')
ax = gca; ax.YScale = 'linear';
grid minor
if exist('funPsd','var') && ~isempty(funPsd)
    yyaxis right
    if isfield(funPsd,'roi')
        hP2 = plot(funPsd.roi.f,funPsd.roi.psd); hold on
    else
        tmp = cat(4,funPsd.vec);
        if ~isreal(tmp)
            tmp = conj(tmp).*tmp;
        end
        % tmp = mean(cat(4,funPsd.vec),2);
        hP2 = plot(funPsd(1).psd.f,mean(tmp(:,:),2)); hold on
    end
    axis tight
    ylabel('psd')
    ytickformat('%2.0f')
    ax.YAxis(2).Scale = 'log';
    hP2.Color = 'k';
    ax.YAxis(2).Color = hP2.Color;
end
xlim(fLim)
fOffset1 = fLim(1)+0.1*range(fLim)+svdStruct.w;
yLim = ylim;

hP1.Color = 'r';
ax.YAxis(1).Color = hP1.Color;
plot([-1 1].*svdStruct.w + fOffset1,(yLim(2)-range(yLim)*0.1).*[1 1],'Color',hP1.Color,'LineWidth',3,'LineStyle','-')

yyaxis right
fOffset2 = fLim(1)+0.1*range(fLim)+funPsd(1).psd.w;
plot([-1 1].*funPsd(1).psd.w + fOffset2,(yLim(2)-range(yLim)*0.2).*[1 1],'Color',hP2.Color,'LineWidth',3,'LineStyle','-')

if augFac>1
    title(['coherence spectrum; sub-' sub '; ses-' ses '; ' num2str(augFac) 'runs; K=' num2str(K)])
else
    title(['coherence spectrum; sub-' sub '; ses-' ses '; run-' run '; K=' num2str(K)])
end
% num2str(peakFreq(fInd),'%.3f') 'Hz'

if exist('peakFreq','var') && ~isempty(peakFreq)
    
    %% Plot underlay
    F{end+1} = figure('WindowStyle','docked');
    if exist('funPsd','var') && ~isempty(funPsd(1)) && isfield(funPsd(1),'imMean') && ~isempty(funPsd(1).imMean)
        imagesc(mean(funPsd.imMean,6))
    else
        imagesc(svdStruct.tsMean);
    end
    ax = gca; ax.Colormap = gray; ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off'; ax.DataAspectRatio = [1 1 1];
    title(['timeseries mean; sub-' sub '; ses-' ses '; run-' run]);
    pos = F{2}.Position;
    drawnow

    %% Plot peak frequencies
    for fInd = 1:length(peakFreq)
        [~,b] = min(abs(f-peakFreq(fInd)));
        %% Plot peak line
        figure(F{1})
        if strcmp(F{1}.Children.YAxisLocation,'right'); yyaxis(F{1}.Children,'left'); end
        plot([1 1].*f(b),ylim,'r')
        yLim = ylim;
        % plot(f(b)+[-1 1].*svdStruct.w,yLim(1).*[1 1],'-g','LineWidth',3)
        text(f(b), yLim(2), num2str(f(b),'%.4fHz'),'VerticalAlignment','top')

        %% Plot sv
        F{end+1} = figure('WindowStyle','docked');
        switch imType
            case 'psd'
                error('code that')
                curIm = vec2vol(funPsd); curIm = mean(curIm.vol(:,:,:,b,:,:),6);
                hIm = imagesc(curIm);
                ax = gca;
                ax.Colormap = jet;
                cb = colorbar;
                ylabel(cb,'psd')
                ax.ColorScale = 'log';
                % cLim = curIm(logical(funPsd.roi.mask)); cLim = [1 max(cLim)];
                % axIm{end}.CLim = cLim;
            case 'svMag'
                curIm = nan(size(svdStruct.mask));
                curIm(logical(svdStruct.mask)) = svdStruct.sp(:,b,1);
                hIm = imagesc(abs(curIm));
                ax = gca;
                ax.Colormap = jet;
                ax.ColorScale = 'log';
                % ax.ColorScale = 'linear';
                cb = colorbar;
                ylabel(cb,'spatial sv weigth mag')
                if isfield(funPsd(1),'roi') && ~isempty(funPsd(1).roi)
                    cLim = abs(curIm(logical(funPsd(1).roi.mask))); cLim = [min(cLim) max(cLim)];
                    ax.CLim = cLim;
                else
                    cLim = ax.CLim;
                end
                cbTicks = cb.Ticks; if cbTicks(1) ~= cLim(1); cbTicks = [cLim(1) cbTicks]; end; if cbTicks(end) ~= cLim(2); cbTicks = [cbTicks cLim(2)]; end; cb.Ticks = cbTicks; cb.TickLabels(2:end-1) = {''};
                cb.TickLabels{1} = 'min'; cb.TickLabels{end} = 'max';
            case 'svPhase'
                curIm = nan(size(svdStruct.mask));
                curIm(logical(svdStruct.mask)) = svdStruct.sp(:,b,1);
                hIm = imagesc(wrapToPi(angle(curIm) - angle(mean(curIm(:),'omitnan'))));
                % alpha = log(abs(curIm));
                alpha = abs(curIm);
                alpha = alpha - min(alpha(:));
                alpha = alpha ./ max(alpha(:));
                hIm.AlphaData = alpha;
                ax = gca;
                ax.Color = [1 1 1].*0.5;
                ax.Colormap = hsv;
                cb = colorbar;
                hYl = ylabel(cb,'spatial sv weigth phase');
                cLim = [-pi pi];
                ax.CLim = cLim;
                cb.Ticks = [-pi -pi/2 0 pi/2 pi];
                cb.TickLabels = {'-pi' '-pi/2' '0' 'pi/2' 'pi'};
        end
        ax.DataAspectRatio = [1 1 1];
        ax.XTick = []; ax.YTick = [];
        title(['spatial singular vector; ' imType '; sub-' sub '; ses-' ses '; run-' run ' K=' num2str(K) '; ' num2str(f(b),'%.3fHz')])
        ax.Position = F{2}.Children.Position;
        drawnow
    end

    %% Link axes
    tmpF = [F{2:end}];
    linkaxes(findobj(cat(1,tmpF.Children),'type','axes'))
end


