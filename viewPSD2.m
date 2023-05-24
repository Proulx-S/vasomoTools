function hF = viewPSD2(funPsd,f0,mask,fpass,id)


%% Prepare some stuff
if exist('mask','var') && ~isempty(mask)
    roiFlag = 0;
    maskC = getMaskOutline(mask,5);
elseif isfield(funPsd,'roi')
    roiFlag = 1;
    mask = funPsd.roi.mask;
    maskC = getMaskOutline(mask,5);
else
    error('code that')
end

%% Prepare figure
hF = figure('WindowStyle','docked');
% tl = tiledlayout(3,4);
tl = tiledlayout(5,4);
tl.TileSpacing = "tight"; tl.Padding = "tight";
tl.TileIndexing = 'rowmajor';

%% Brain
nexttile(1,[2 1])
imagesc(funPsd.tMean)
ax = gca;
ax.Colormap = gray;
ax.YTick = []; ax.XTick = [];
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
title('timeseries mean')

%% Mask
nexttile(2,[2 1])
imagesc(funPsd.tMean)
ax = gca;
ax.Colormap = gray;
ax.YTick = []; ax.XTick = [];
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
hold on
hMask1 = plot(maskC);
hMask1.FaceColor = 'r';
hMask1.FaceAlpha = 0.1;
hMask1.EdgeColor = 'none';
title('mask')

%% Normalization
nexttile(3,[2 1])
if isfield(funPsd.psd,'norm')
    normFlag = 1;
    normFact = nan(size(funPsd.vol2vec));
    normFact(funPsd.vol2vec) = funPsd.psd.norm.fact;
    imagesc(normFact)
    hold on
    hMask2 = plot(maskC);
    hMask2.FaceColor = 'none';
    hMask2.EdgeColor = 'k';
    ax = gca;
    ax.YTick = []; ax.XTick = [];
    ax.ColorScale = 'log'; ax.Colormap = jet;
    ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ylabel(colorbar,'noise floor psd')
    title('normalization factor')
else
    normFlag = 0;
    ax = gca;
    ax.Visible = 'off';
end

%% f0
nexttile(4,[2 1])
f = funPsd.psd.f;
[~,f0Ind] = min(abs(f-f0'),[],2);
tmpPsd = vec2vol(funPsd);
imagesc(tmpPsd.vol(:,:,:,f0Ind));
hold on
hMask3 = plot(maskC);
hMask3.FaceColor = 'none';
hMask3.EdgeColor = 'k';
ax = gca;
ax.YTick = []; ax.XTick = [];
ax.ColorScale = 'log'; ax.Colormap = jet;
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
ylabel(colorbar,'raw psd')
title(['f0=' num2str(f0,'%0.3f') 'Hz'])

%% Spectrum
% nexttile(9,[1 4])
nexttile(17,[1 4])
legLabel = {};
h = {};
if roiFlag
    psd = funPsd.roi.psd;
    psdErr = funPsd.roi.psdErr;
    W = funPsd.roi.w;
    h{end+1} = plot(f,psd,'k'); legLabel{end+1} = 'mean spectrum';
    hold on
    h{end+1} = plot(f,psdErr,'r'); legLabel{end+1} = '95%CI';
    h{end} = h{end}(1);
else
    tmpPsd = vec2vol(funPsd);
    tmpPsd = vol2vec(tmpPsd,mask,1);
    psd = mean(tmpPsd.vec,2);
    W = funPsd.psd.w;
    h{end+1} = plot(f,psd,'k'); legLabel{end+1} = 'mean spectrum';
    hold on
end
ax = gca;
ax.YScale = 'log';
ax.YAxisLocation = 'right';
axis tight
ax.XLim = fpass;
yLim = [min(psd) max(psd)];
yLim(1) = mean(psd(f>3.5));
yLim(1) = exp(log(yLim(1)) - range(log(yLim))*0.05);
ylim(yLim);
grid on
ax.XMinorGrid = 'on';
xlabel('Hz')
if normFlag
    ylabel('normalized psd')
else
    ylabel('raw psd')
end
h{end+1} = plot(f0.*[1 1],yLim,':r'); legLabel{end+1} = 'f0';
h{end+1} = plot(f0+W.*[-1 1],yLim(1).*[1 1],'g','LineWidth',5); legLabel{end+1} = 'BW';
legend([h{:}],legLabel,'box','off')
if normFlag
    title('spectrum average within mask, normalized voxel-wise to noise=1')
else
    title('spectrum average within mask')
end

%% Title
titleStr1 = {strjoin(id,'; ')};
if normFlag
    titleStr2 = {'normalized voxel-wise'};
else
    titleStr2 = {'raw'};
end
titleStr2{end+1} = ['W=' num2str(funPsd.psd.w,'%0.4f')];
titleStr2{end+1} = ['K=' num2str(funPsd.psd.param.tapers(2))];
titleStr2 = {strjoin(titleStr2,'; ')};

titleStr = [titleStr1 titleStr2];
tlTitle = title(tl,titleStr,'interpreter','none');
drawnow

