function [hF,tl,axSpec,maskC] = viewPSD2(funPsd,f0,mask,fpass,id,threshFlag)
%% Prepare some stuff
if ~exist("threshFlag",'var') || isempty(threshFlag)
    threshFlag = false;
end
if exist('mask','var') && ~isempty(mask)
    roiFlag = 0;
    maskC = getMaskOutline(mask,5);
elseif isfield(funPsd,'roi')
    roiFlag = 1;
    mask = funPsd.roi.mask;
    maskC = getMaskOutline(mask,5);
else
    roiFlag = 0;
    maskC = [];
end
if ~exist('f0','var')
    f0=0.1;
end
if ~exist('id','var')
    id = {};
end

%% Prepare figure
hF = figure('WindowStyle','docked');
% tl = tiledlayout(3,4);
tl = tiledlayout(5,4);
tl.TileSpacing = "tight"; tl.Padding = "tight";
tl.TileIndexing = 'rowmajor';

%% Brain
nexttile(1,[2 1])
im = funPsd.tMean;
imagesc(im)
ax = gca;
ax.Colormap = gray;
ax.YTick = []; ax.XTick = [];
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
ax.CLim = prctile(im(:),[0 99]);
title('timeseries mean')

%% Mask
nexttile(2,[2 1])
if isempty(maskC)
    ax = gca;
    ax.Visible = 'off';
else
    im = funPsd.tMean;
    imagesc(im)
    ax = gca;
    ax.Colormap = gray;
    ax.YTick = []; ax.XTick = [];
    ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ax.CLim = prctile(im(:),[0 99]);
    hold on
    hMask1 = plot(maskC);
    hMask1.FaceColor = 'r';
    hMask1.FaceAlpha = 0.1;
    hMask1.EdgeColor = 'none';
    title('mask')
end

%% Normalization
% drawnow
nexttile(3,[2 1])
if isfield(funPsd.psd,'norm') && ~isfield(funPsd.psd.norm,'fact2')
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
tmpIm = tmpPsd.vol(:,:,:,f0Ind);
hIm = imagesc(tmpIm);
hold on
if ~isempty(maskC)
    hMask3 = plot(maskC);
    hMask3.FaceColor = 'none';
    hMask3.EdgeColor = 'k';
end
ax = gca;
ax.YTick = []; ax.XTick = [];
ax.ColorScale = 'log'; ax.Colormap = jet;
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
if normFlag
    ylabel(colorbar,'normalized psd')
else
    ylabel(colorbar,'raw psd')
end

if threshFlag && isfield(funPsd.psd,'aboveNoiseInd')
    tmp = zeros(size(mask));
    tmp(:) = funPsd.psd.aboveNoiseInd(:,f0Ind);
    hIm.AlphaData = tmp;
    ax.CLim = prctile(tmpIm(logical(tmp)&logical(mask)),[0 95]);
    ax.Color = 'k';
    hMask3.EdgeColor = 'w';
    title(['f0=' num2str(f0,'%0.3f') 'Hz (thresholded)'])
else
    title(['f0=' num2str(f0,'%0.3f') 'Hz'])
end


%% Spectrum
% nexttile(9,[1 4])
axSpec = nexttile(17,[1 4]);
legLabel = {};
h = {};
if roiFlag
    psd = funPsd.roi.psd;
    psdErr = funPsd.roi.psdErr;
    W = funPsd.roi.w;
    h{end+1} = plot(f,psd,'k'); legLabel{end+1} = 'mean spectrum';
    hold on
%     yyaxis right
%     h{end+1} = plot(f,diff(psdErr,[],2)./psd);
    if ~isempty(psdErr)
        h{end+1} = plot(f,psdErr,'r'); legLabel{end+1} = '95%CI';
        h{end} = h{end}(1);
    end
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
if exist('fpass','var') && ~isempty(fpass)
    ax.XLim = fpass;
end
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
% drawnow

