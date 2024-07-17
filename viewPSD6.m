function [hF,tl,axSpec,maskC] = viewPSD6(funPsd,f0,t0,roi,mask,fpass,id,threshFlag)
%% Prepare some stuff
if ~exist("threshFlag",'var') || isempty(threshFlag)
    threshFlag = false;
end

if ~exist('roi','var') || isempty(roi)
    roi = [];
end
roiC = getMaskOutline(roi,5);

if ~exist('mask','var') || isempty(mask)
    if isfield(funPsd,'vol2vec') && ~isempty(funPsd.vol2vec)
        mask = funPsd.vol2vec;
    else
        mask = [];
    end
end
% maskC = getMaskOutline(mask,5);

% if exist('mask','var') && ~isempty(mask)
%     roiFlag = 0;
%     maskC = getMaskOutline(mask,5);
% elseif isfield(funPsd,'roi')
%     roiFlag = 1;
%     mask = funPsd.roi.mask;
%     maskC = getMaskOutline(mask,5);
% else
%     roiFlag = 0;
%     maskC = [];
% end
if ~exist('f0','var') || isempty(f0)
    f0=0.1;
    f0w = nan;
else
    if iscell(f0)
        f0c = nan(1,length(f0));
        f0w = nan(1,length(f0));
        for i = 1:length(f0)
            f0c(i) = mean(f0{i});
            if length(f0{i})==2
                f0w(i) = diff(f0{i});
            end
        end
        f0 = f0c; clear f0c
    end
    [f0,b] = sort(f0);
    if exist('f0w','var')
        f0w = f0w(b);
    else
        f0w = nan;
    end
end
if ~exist('id','var')
    id = {};
end


%% Prepare figure
hF = figure('WindowStyle','docked');
% tl = tiledlayout(3,4);
tl = tiledlayout(5,5);
tl.TileSpacing = "tight"; tl.Padding = "tight";
tl.TileIndexing = 'rowmajor';
axLink = {};

%% Brain
nexttile(1,[2 1])
if isfield(funPsd,'imMean') && ~all(size(funPsd.imMean)==1)
    im = squeeze(mean(funPsd.imMean,6));
    imagesc(im); hold on
    ax = gca;
    ax.Colormap = gray;
    ax.YTick = []; ax.XTick = [];
    ax.DataAspectRatio = [1 1 1];
    ax.CLim = prctile(im(:),[0 99]);
    title('timeseries mean')
    if ~isempty(roiC)
        plot(roiC,'FaceColor','none','EdgeColor','w');
        plot(roiC,'FaceColor','none','EdgeColor','k','LineStyle','--');
    end
else
    ax = gca;
    ax.Visible = 'off';
end
axLink{end+1} = ax;

% %% Mask and roi
% if ~isempty(roi)
%     roiC = getMaskOutline(roi,5);
% 
%     % if isempty(mask)
%     %     tmp = funPsd.vol2vec(funPsd.vol2vec);
%     %     tmp(any(isnan(funPsd.vec),1) | all(funPsd.vec==0,1)) = false;
%     %     funPsd.vol2vec(funPsd.vol2vec) = tmp;
%     %     funPsd.vec(~tmp) = [];
%     %     % if isfield(funPsd,'normFac') && ~isempty(funPsd.normFac)
%     %     %     mask = ~any(isnan(funPsd.normFac),6);
%     %     % else
%     %     %     tmpPsd = vec2vol(funPsd);
%     %     %     mask = all(~isnan(tmpPsd.vol),4) & ~all(tmpPsd.vol==0,4);
%     %     %     clear tmpPsd
%     %     % end
%     % end
% 
%     nexttile(2,[2 1])
%     if isempty(maskC)
%         ax = gca;
%         ax.Visible = 'off';
%     else
%         im = squeeze(mean(funPsd.imMean(:,:,:,:),4));
%         imagesc(im)
%         ax = gca;
%         ax.Colormap = gray;
%         ax.YTick = []; ax.XTick = [];
%         ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
%         ax.CLim = prctile(im(:),[0 99]);
%         hold on
%         hMask1 = plot(maskC);
%         hMask1.FaceColor = 'r';
%         hMask1.FaceAlpha = 0.1;
%         hMask1.EdgeColor = 'none';
%         title('mask')
%     end
% end

%% Scaling
% drawnow
nexttile(2,[2 1])
if isfield(funPsd,'scale') && ~isempty(funPsd.scale)
% if ( isfield(funPsd.psd,'norm') && ~isfield(funPsd.psd.norm,'fact2') ) || isfield(funPsd,'normFac')
    normFlag = 1;
    normFact = funPsd.scale;
    % if isfield(funPsd.psd,'norm') & ~isfield(funPsd,'normFac')
    %     normFact = nan([size(funPsd.vec,[3 4]) size(funPsd.vol2vec)]);
    %     normFact(:,:,funPsd.vol2vec) = permute(funPsd.psd.norm.fact,[3 4 2 1]);
    %     normFact = permute(normFact,[3 4 5 1 2]);
    % elseif ~isfield(funPsd.psd,'norm') & isfield(funPsd,'normFac')
    %     normFact = funPsd.normFac;
    % else
    %     dbstack; error('X');
    % end

    if ~isreal(normFact)
        normFact = conj(normFact).*normFact;
    end
    imagesc(squeeze(mean(normFact(:,:,:,:,:),5)))
    hold on
    ax = gca;
    ax.YTick = []; ax.XTick = [];
    ax.ColorScale = 'log'; ax.Colormap = jet;
    ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    if ~isempty(roiC)
        plot(roiC,'FaceColor','none','EdgeColor','w');
        plot(roiC,'FaceColor','none','EdgeColor','k','LineStyle','--');
    end
    ylabel(colorbar,'noise floor psd')
    title('normalization factor')
else
    normFlag = 0;
    ax = gca;
    ax.Visible = 'off';
end;
axLink{end+1} = ax;

%% f0
f = funPsd.psd.f;
allF0 = f0;
allF0w = f0w;
for fInd = 1:length(allF0)
    f0 = allF0(fInd);
    f0w = allF0w(fInd);
    if fInd==1
        nexttile(3,[2 1])
    elseif fInd>2
        warning('too many f0 to plot, skipping')
        break
    else
        nexttile(8+fInd-1,[2 1])
    end
    if size(funPsd.vec,2)>1
        if isnan(f0w)
            [~,f0Ind] = min(abs(f-f0'),[],2);
            f0 = f(f0Ind);
        else
            [~,f0l_ind] = min(abs(f-(f0-f0w/2)),[],2);
            [~,f0u_ind] = min(abs(f-(f0+f0w/2)),[],2);
            f0l = f(f0l_ind);
            f0u = f(f0u_ind);
            f0 = mean([f0l f0u]);
            f0Ind = false(size(f));
            f0Ind(f0l_ind:f0u_ind) = true;
        end
        if ~isempty(funPsd.vec)
            tmpIm = nan(size(funPsd.vol2vec));
            tmpIm(funPsd.vol2vec) = mean(mean(conj(funPsd.vec(f0Ind,:,:)).*funPsd.vec(f0Ind,:,:),3),1);
        else
            error('code that')
            tmpIm = mean(conj(funPsd.vol(:,:,:,:)).*funPsd.vol(:,:,:,:),4);
        end
        imagesc(tmpIm);
        hold on
        ax = gca;
        ax.YTick = []; ax.XTick = [];
        ax.ColorScale = 'log';
        ax.DataAspectRatio = [1 1 1];
        if normFlag
            ylabel(colorbar,'normalized psd')
        else
            ylabel(colorbar,'raw psd')
        end
        if ~isempty(roiC)
            plot(roiC,'FaceColor','none','EdgeColor','w');
            plot(roiC,'FaceColor','none','EdgeColor','k','LineStyle','--');
        end

        % if threshFlag && isfield(funPsd.psd,'aboveNoiseInd')
        %     tmp = zeros(size(mask));
        %     tmp(:) = funPsd.psd.aboveNoiseInd(:,f0Ind);
        %     hIm.AlphaData = tmp;
        %     ax.CLim = prctile(tmpIm(logical(tmp)&logical(mask)),[0 95]);
        %     ax.Color = 'k';
        %     hMask3.EdgeColor = 'w';
        %     title(['f0=' num2str(f0,'%0.3f') 'Hz (thresholded)'])
        % else
        title(['Power at f0=' num2str(f0,'%0.3f') 'Hz'])
        % end
        if normFlag
            ax.CLim(1) = 1;
            if exist('mask','var') && ~isempty(mask)
                sz = size(mask,[1 2 3]);
                tmpMask = false(sz); tmpMask(ceil(sz(1)*0.1):round(sz(1)*0.9),ceil(sz(2)*0.1):round(sz(2)*0.9),ceil(sz(3)*0.1):round(sz(3)*0.9)) = true;
                ax.CLim(2) = max(tmpIm(logical(mask)&tmpMask));
            end
            ax.Colormap = hot;
        else
            ax.Colormap = jet;
        end
    else
        ax = gca;
        ax.Visible = 'off';
    end
%     %%% scale every image to the first
%     if fInd==1
%         cLim  = ax.CLim;
%     else
%         ax.CLim = cLim;
%     end
end
f0 = allF0;

axLink{end+1} = ax;


%% Plot first singular vector
if length(f0)==1 && isfield(funPsd,'svd')
    [~,f0Ind] = min(abs(funPsd.svd.f - f0));
    sv = funPsd.svd.u(:,1,f0Ind);
    
    nexttile(4,[2 1])

    im = nan(size(funPsd.vol2vec));
    im(funPsd.vol2vec) = abs(sv);

    % sz = size(tmpMask);
    tmpMask = funPsd.vol2vec;
    tmpMask([1:3 end-2:end],:) = false;
    tmpMask(:,[1:3 end-2:end]) = false;
    cLim = [min(im(tmpMask)) max(im(tmpMask))];
    imagesc(im,cLim);
    
    ax = gca;
    ax.YTick = []; ax.XTick = [];
    % ax.ColorScale = 'log';
    ax.DataAspectRatio = [1 1 1];
    ax.Colormap = hot;
    % if ~isempty(roiC)
    %     plot(roiC,'FaceColor','none','EdgeColor','w');
    %     plot(roiC,'FaceColor','none','EdgeColor','k','LineStyle','--');
    % end
    title(['1st singular vector at f0=' num2str(f0,'%0.3f') 'Hz'])
    ylabel(colorbar,'squared mag')

    axLink{end+1} = ax;


    nexttile(5,[2 1])

    im = nan(size(funPsd.vol2vec));
    im(funPsd.vol2vec) = angle(sv);
    cLim = [-pi pi];
    imagesc(im,cLim);

    ax = gca;
    ax.YTick = []; ax.XTick = [];
    ax.DataAspectRatio = [1 1 1];
    ax.Colormap = hsv;
    % if ~isempty(roiC)
    %     plot(roiC,'FaceColor','none','EdgeColor','w');
    %     plot(roiC,'FaceColor','none','EdgeColor','k','LineStyle','--');
    % end
    title(['1st singular vector at f0=' num2str(f0,'%0.3f') 'Hz'])
    ylabel(colorbar,'phase')

    scale = nan(size(funPsd.vol2vec));
    scale(funPsd.vol2vec) = abs(sv);
    scale = scale.^2;
    % scale = log(scale);
    % scale = sqrt(scale);
    % scale = log(scale);
    scale = scale - min(scale(tmpMask));
    scale = scale ./ max(scale(tmpMask));
    ax.Children.AlphaData = scale;
    ax.Color = [0.5 0.5 0.5];

    axLink{end+1} = ax;

else
    ax = nexttile(4,[2 1]);
    ax.Visible = 'off';
    ax = nexttile(5,[2 1]);
    ax.Visible = 'off';
end

% if length(allF0)==1
%     nexttile(3,[2 1])
% end

linkaxes([axLink{:}])

%% Spectrum
axSpec = nexttile(11,[3 5]);
% axSpec = nexttile(17,[1 4]);
legLabel = {};
h = {};
if ~isempty(funPsd.vec)
    psd = funPsd.vec;
    if ~isempty(roi) && any(funPsd.vol2vec(:)~=roi(:))
        if any(roi(~funPsd.vol2vec)); dbstack; error('roi contains voxels that do not have data'); end
        psd = psd(:,logical(roi(funPsd.vol2vec)),:);
    end
else
    dbstack; error('code that')
end
yyaxis left
psd = mean(conj(psd(:,:)).*psd(:,:),2);
W = funPsd.psd.w;
h{end+1} = plot(f,squeeze(psd));
if ~isempty(roi)
    legLabel{end+1} = 'roi mean spectrum';
else
    legLabel{end+1} = 'image mean spectrum';
end
hold on

% if roiFlag
%     psd = funPsd.roi.psd;
%     psdErr = funPsd.roi.psdErr;
%     W = funPsd.roi.w;
%     h{end+1} = plot(f,psd,'k'); legLabel{end+1} = 'mean spectrum';
%     hold on
% %     yyaxis right
% %     h{end+1} = plot(f,diff(psdErr,[],2)./psd);
%     if ~isempty(psdErr)
%         h{end+1} = plot(f,psdErr,'r'); legLabel{end+1} = '95%CI';
%         h{end} = h{end}(1);
%     end
% else
%     tmpPsd = vec2vol(funPsd);
%     if ~exist('mask','var')
%         mask = true(size(tmpPsd.vol,1:3));
%     end
%     tmpPsd = vol2vec(tmpPsd,mask,1);
%     if ~isreal(tmpPsd.vec)
%         tmpPsd.vec = conj(tmpPsd.vec).*tmpPsd.vec;
%     end
%     psd = mean(tmpPsd.vec(:,:),2);
%     % psd = mean(psd,4);
%     W = funPsd.psd.w;
%     h{end+1} = plot(f,squeeze(psd),'k'); legLabel{end+1} = 'mean spectrum';
%     hold on
% end
ax = gca;
ax.YScale = 'log';
% ax.YAxisLocation = 'right';
axis tight
if exist('fpass','var') && ~isempty(fpass)
    ax.XLim = fpass;
end
yLim = [min(psd(2:end)) max(psd(2:end))];
% if max(f)>3.5
%     yLim(1) = mean(psd(f>3.5));
% end
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
for fInd = 1:length(allF0)
    f0 = allF0(fInd);
    f0w = allF0w(fInd);
    if ~isnan(f0w)
        f0l = f0-f0w/2;
        f0u = f0+f0w/2;
        hPatch = patch([f0l f0u f0u f0l f0l], yLim([1 1 2 2 1]),[1 1 1]*0);
        hPatch.EdgeColor = 'none';
        hPatch.FaceAlpha = 0.12;
        uistack(hPatch,'bottom')
    else
        [~,f0Ind] = min(abs(f-f0'),[],2);
        f0 = f(f0Ind);
        if fInd==1
            h{end+1} = plot(f0.*[1 1],yLim,'-r'); legLabel{end+1} = 'f0';
            % h{end+1} = plot(f0+W.*[-1 1],yLim(1).*[1 1],'-g','LineWidth',5); legLabel{end+1} = 'BW';
        else
            plot(f0.*[1 1],yLim,'-r');
            plot(f0+W.*[-1 1],yLim(1).*[1 1],'-g','LineWidth',5);
        end
    end
end
f0 = allF0;

%%% bandwidth
X = f0+W.*[-1 1]; X = X - mean(X) + 0.1;
Y = [1 1];
h{end+1} = plot(X,Y,'-g','LineWidth',5);
legLabel{end+1} = 'BW';



% if normFlag
%     title('spectrum average within roi, normalized voxel-wise to noise=1')
% else
%     title('spectrum average within mask')
% end

if isfield(funPsd,'svd')
    % yyaxis left
    yyaxis right
    h{end+1} = plot(squeeze(funPsd.svd.f),squeeze(funPsd.svd.coh(:,1,:)));
    legLabel{end+1} = 'coherence spectrum';
    
    K = size(funPsd.svd.s,2);
    ylim([0 1])
    h{end+1} = plot(squeeze(funPsd.svd.f([1 end])),[1 1].*1/K,'--','Color',h{end}.Color);
    legLabel{end+1} = 'baseline';

    ylabel('coherence')

    if isfield(funPsd,'scale') && ~isempty(funPsd.scale)
        yyaxis left
        yLim = ylim;
        m = log(yLim(2));
        z = log(1);
        B = (K*z-m)/(K-1);
        yLim(1) = exp(B);
        ylim(yLim)
    end

    
end


legend([h{:}],legLabel,'box','off')





%% Title
if iscell(id)
    titleStr1 = {strjoin(id,'; ')};
else
    titleStr1 = {id};
end
if isempty(titleStr1{1})
    titleStr1 = strsplit(funPsd.fspec,filesep);
    label = 'sub-';
    tmp = strsplit(titleStr1{find(contains(titleStr1,label),1)},'_'); tmp = strsplit(tmp{find(contains(tmp,label),1)},'-');
    sub = tmp{2};
    label = 'ses-';
    tmp = strsplit(titleStr1{find(contains(titleStr1,label),1)},'_'); tmp = strsplit(tmp{find(contains(tmp,label),1)},'-');
    ses = tmp{2};
    label = 'run-';
    tmp = strsplit(titleStr1{find(contains(titleStr1,label),1)},'_'); tmp = strsplit(tmp{find(contains(tmp,label),1)},'-');
    run = tmp{2};
    titleStr1 = {['sub-' sub]};
    titleStr1{end+1} = ['ses-' ses];
    if funPsd.nruns==1
        titleStr1{end+1} = ['run-' run];   
    end
    titleStr1 = {strjoin(titleStr1,'_')};
end
if normFlag
    titleStr2 = {'normalized voxel-wise'};
else
    titleStr2 = {'raw'};
end
titleStr2{end+1} = ['W=' num2str(funPsd.psd.w,'%0.4f')];
titleStr2{end+1} = ['K=' num2str(funPsd.psd.param.tapers(2))];
if funPsd.nruns>1
    titleStr2{end+1} = ['nRuns=' num2str(funPsd.nruns)];
end

titleStr2 = {strjoin(titleStr2,'; ')};

titleStr = strjoin([titleStr1 titleStr2],'; ');
tlTitle = title(tl,titleStr,'interpreter','none');
% drawnow

