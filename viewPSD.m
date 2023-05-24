function [funTmpName,funPsdNorm,hF] = viewPSD(funPsd,f0,mask,fpass,id,slc,brain)
normMethod = 'psdNoise'; % 'psdNoise', 'smoothPsdAv', 'psdAv' or 'none'
globalNormMethod = 'to1'; % 'toImageAverage' or 'to1';
cLimMethod = 'fullRange';% 'fullRange', '99prctile' or 'highFloor'
fwhm = 5;
if iscell(funPsd)
    psdStruct = funPsd{2};
    funPsd = funPsd{1};
end
if ischar(f0)
    switch f0
        case 'all'
            freeviewFlag = 1;
        otherwise
            error('f0 can only be numeric of ''all''')
    end
    f0 = [];
else
    if isempty(f0)
        freeviewFlag = 0;
    else
        freeviewFlag = 1;
    end
end
if ~exist('id','var') || isempty(id)
    id = '';
elseif isnumeric(id)
    id = num2str(id);
end
if ~exist('mask','var') || isempty(mask)
    mask = funPsd.vol2vec;
    switch funPsd.vol2vecFlag
        case 'allInclusiveMask'
            doMaskOutline = false;
        case 'customMask'
            doMaskOutline = true;
        otherwise
            error('X')
    end
else
    doMaskOutline = true;
end
if ~exist('slc','var') || isempty(slc)
    slc = round(funPsd.depth*0.5);
end
f = funPsd.psd.f;
w = funPsd.psd.w;

% %% Perform spatial normalization
% switch normMethod
%     case 'psdNoise'
%         if isempty(funPsd.vec)
%             funPsd = vol2vec(funPsd);
%         end
%         noiseInd = funPsd.psd.f>4 & funPsd.psd.f<funPsd.psd.f(round(0.95*end));
%         normFact = exp( mean(log(funPsd.vec(noiseInd,:)),1) );
%         normFact0 = exp(mean(log(normFact),2));
%     case 'psdAv'
%         if isempty(funPsd.vec)
%             funPsd = vol2vec(funPsd);
%         end
%         normFact = exp( mean(log(funPsd.vec),1) );
%         normFact0 = exp(mean(log(normFact),2));
%     case 'smoothPsdAv'
%         error('this is probably broken')
%         %%% Compute average power and write to disk
%         funNorm = funPsd;
%         if isempty(funNorm.vol)
%             funNorm.vec = mean(log(funNorm.vec),1);
%             funNorm.nframes = 1;
%             normFact0 = exp(funNorm.vec);
%             funNorm = vec2vol(funNorm);
%         else
%             funNorm.vol = mean(log(funNorm.vol),4);
%             funNorm.nframes = 1;
%             tmp = vol2vec(funNorm);
%             normFact0 = exp(tmp.vec); clear tmp
%         end
%         funNormTmpName = [tempname '.nii.gz'];
%         MRIwrite(funNorm,funNormTmpName);
%         
% 
%         %%% Define mask for smoothing and write to disk
%         funMask = funNorm;
%         funMask.vol = ~isnan(funMask.vol);
%         funMaskTmpName = [tempname '.nii.gz'];
%         MRIwrite(funMask,funMaskTmpName);
% 
%         %%% Apply smoothing (mri_fwhm) to normalization map
%         cmd = {'source /usr/local/freesurfer/fs-dev-env-autoselect'};
%         cmd{end+1} = ['mri_fwhm'...
%             ' --i ' funNormTmpName...
%             ' --o ' funNormTmpName...
%             ' --smooth-only'...
%             ' --mask ' funMaskTmpName...
%             ' --fwhm ' num2str(fwhm)];
%         [status,cmdout] = system(strjoin(cmd,'; '));
% 
%         % funNormSm = MRIread(funNormTmpName);
%         % figure('WindowStyle','docked');
%         % imagesc(exp(funNorm.vol(:,:,46)));
%         % ax = gca; ax.ColorScale = 'log';
%         % cLim = clim;
%         % figure('WindowStyle','docked');
%         % imagesc(exp(funNormSm.vol(:,:,46)));
%         % ax = gca; ax.ColorScale = 'log';
%         % clim(cLim);
% 
%         funNorm = MRIread(funNormTmpName);
%         funNorm = vol2vec(funNorm,funPsd.vol2vec);
%         normFact = exp(funNorm.vec); clear funNorm
%     case 'none'
%         sz = size(funPsd.vec);
%         normFact = ones([1 sz(2)]);
%         normFact0 = ones([1 sz(2)]);
%     otherwise
%         error('')
% end
% 
% %%% Apply normalization
% if isempty(funPsd.vec)
%     funPsd = vol2vec(funPsd);
% end
% funPsdNorm = funPsd;
% switch globalNormMethod
%     case 'toImageAverage'
%         funPsdNorm.vec = funPsd.vec./normFact .* normFact0;
%     case 'to1'
%         funPsdNorm.vec = funPsd.vec./normFact;
% end
% funPsdNorm.normFact = normFact;


%% Compute summary statistics for later
tmpMask = logical(mask(funPsdNorm.vol2vec));
if any(all(funPsdNorm.vec==0,1))
    warning('some voxels are all zeros. fixing it but expect wierd stuff')
end
tmpMask(all(funPsdNorm.vec==0,1)) = false;
psdSpecNormAv = exp(mean(log(funPsdNorm.vec(:,tmpMask)),2));

% psdSpecNormEr = exp(prctile(log(funPsdNorm.vec(:,mask(funPsdNorm.vol2vec))),[2.7 97.5],2));
% psdSpecNormEr(:,1) = psdSpecNormEr(:,1) - psdSpecNormAv; psdSpecNormEr(:,2) = psdSpecNormAv - psdSpecNormEr(:,2);

%% Prepare for ploting
hF = figure('WindowStyle','docked');
% Title
if exist('id','var') && ischar(id) && ~isempty(id)
    titleStr     = {['subj' id]};
elseif iscell(id)
    titleStr = {strjoin(id,'; ')};
else
    titleStr = {''};
end
titleStr{end+1} = ['tw=' num2str(funPsd.psd.tw) ', ' 'w=' num2str(w) 'Hz' ', ' 'k=' num2str(funPsd.psd.param.tapers(2))];
% titleStr{end+1} = ['w=' num2str(w) 'Hz'];
% switch normMethod
%     case 'smoothPsdAv'
%         titleStr{end+1} = ['fwhm=' num2str(fwhm) 'mm'];
% end
[~,b,c] = fileparts(funPsd.fspec);
titleStr{end+1} = [b c];
% titleStr = [strjoin(titleStr(1:end-1),'; ') titleStr(end)];

if doMaskOutline
    maskC = getMaskOutline(mask(:,:,slc),5);
end

if ~isempty(f0)
    funPsd = vec2vol(funPsd);
    funPsdNorm = vec2vol(funPsdNorm);
    [~,f0Ind] = min(abs(f-f0'),[],2);
    if length(f0)==1
        tl = tiledlayout(2,4); tl.TileSpacing = "tight"; tl.Padding = "tight";
        tl.TileIndexing = 'rowmajor';

        %% Plot spatial pattern before and after normalization
        xStr = {['f0=' num2str(f0) 'Hz']};
        xStr{end+1} = ['slc' num2str(slc)];
        % Before spatial normalization
        nexttile
        imagesc(squeeze(funPsd.vol(:,:,slc,f0Ind)));
        ylabel(colorbar,['psd @ f0=' num2str(f0) 'Hz'])
        ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
        ax.Colormap = jet;
        tmp = funPsd.vol(:,:,slc,f0Ind); tmp = tmp(logical(mask));
        switch cLimMethod
            case 'highFloor'
                cLim = exp(prctile(log(tmp),[0.5 100]));
            case '99prctile'
                cLim = exp(prctile(log(tmp),[0.5 99.5]));
            case 'fullRange'
                cLim = [min(tmp) max(tmp)];
            otherwise
                error('X')
        end
        ax.CLim = cLim;
        xlabel({strjoin(xStr,'; ') 'before spatial norm'})
        hold on;
        if doMaskOutline; plot(maskC,'FaceColor','none','EdgeColor','m'); end

        % Average power
        nexttile
        tmp = funPsdNorm; tmp.vec = normFact; tmp.nframes = 1; tmp = vec2vol(tmp);
        imagesc(squeeze(tmp.vol(:,:,slc))); clear tmp;
        ylabel(colorbar,'full spectrum psd')
        ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
        ax.Colormap = jet;
        ax.CLim = cLim;
        if strcmp(normMethod,'psdNoise')
            xlabel('high frequency noise psd')
        else
            xlabel('full spectrum psd')
        end
        hold on;
        if doMaskOutline; plot(maskC,'FaceColor','none','EdgeColor','m'); end

        % Normalization factor
        nexttile
        tmp = funPsdNorm; tmp.vec = normFact; tmp.nframes = 1; tmp = vec2vol(tmp);
        imagesc(squeeze(tmp.vol(:,:,slc))); clear tmp;
        ylabel(colorbar,'full spectrum psd')
        ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
        ax.Colormap = jet;
        ax.CLim = cLim;
        xlabel('normalization factor')
        hold on;
        if doMaskOutline; plot(maskC,'FaceColor','none','EdgeColor','m'); end

        %     % After spatial normalization (same scale)
        %     nexttile
        %     imagesc(squeeze(funPsdNorm.vol(:,:,slc,f0Ind)));
        %     ylabel(colorbar,['psd @ f0=' num2str(f0) 'Hz'])
        %     ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
        %     ax.CLim = cLim;
        %     xlabel({strjoin(xStr,'; ') 'after spatial norm' '(same scale)'})
        %     hold on; plot(maskC,'FaceColor','none','EdgeColor','m')

        % After spatial normalization (auto scale)
        nexttile
        imagesc(squeeze(funPsdNorm.vol(:,:,slc,f0Ind)));
        ylabel(colorbar,['psd @ f0=' num2str(f0) 'Hz'])
        ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
        ax.Colormap = jet;
        xlabel({strjoin(xStr,'; ') 'after spatial norm'})
        hold on;
        if doMaskOutline; plot(maskC,'FaceColor','none','EdgeColor','m'); end
        tmp = funPsdNorm.vol(:,:,slc,f0Ind); tmp = tmp(logical(mask));
        switch normMethod
            case 'psdNoise'
                cLim = [1 exp(prctile(log(tmp),95))];
            otherwise
                switch cLimMethod
                    case 'highFloor'
                        cLim = exp(prctile(log(tmp),[0.5 100]));
                    case '99prctile'
                        cLim = exp(prctile(log(tmp),[0.5 99.5]));
                    case 'fullRange'
                        cLim = [min(tmp) max(tmp)];
                    otherwise
                        error('X')
                end

        end
        ax.CLim = cLim;
    else
        if ~isempty(brain)
            if iscell(brain)
                tl = tiledlayout(2,length(f0)+length(brain)); tl.TileSpacing = "tight"; tl.Padding = "tight";
                for iii = 1:length(brain)
                    nexttile
                    imagesc(brain{iii})
                    ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = [];
                    ax.CLim = prctile(brain{iii}(:),[1 99]);
                    ax.Colormap = gray;
                    hold on;
                    if doMaskOutline; plot(maskC,'FaceColor','none','EdgeColor','r'); end
                end
            else
                tl = tiledlayout(2,length(f0)+1); tl.TileSpacing = "tight"; tl.Padding = "tight";
                nexttile
                imagesc(brain)
                ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = [];
                ax.CLim = prctile(brain(:),[1 99]);
                ax.Colormap = gray;
                hold on;
                if doMaskOutline; plot(maskC,'FaceColor','none','EdgeColor','r'); end
            end
        else
            tl = tiledlayout(2,length(f0)); tl.TileSpacing = "tight"; tl.Padding = "tight";
        end
        tl.TileIndexing = 'rowmajor';
        % After spatial normalization (auto scale)
        for iii = 1:length(f0Ind)
            nexttile
            imagesc(squeeze(funPsdNorm.vol(:,:,slc,f0Ind(iii))));
            ylabel(colorbar,['psd @ f0=' num2str(f0(iii)) 'Hz'])
            ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = [];
%             ax.YDir = 'normal';
            ax.Colormap = jet;
%             xlabel(strjoin(xStr,'; '),'interpreter','none')
            hold on;
            if doMaskOutline; plot(maskC,'FaceColor','none','EdgeColor','m'); end
            cLim = ax.CLim;
        end
    end
else
    %% Plot brain and mask
    xStr = {['slc' num2str(slc)]};
    
    nexttile
    imagesc(brain.vol(:,:,slc))
    ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
    ax.Colormap = gray;
    xlabel({strjoin(xStr,'; ') 'after spatial norm' 'auto scale'})
    hold on;
    if doMaskOutline; plot(maskC,'FaceColor','none','EdgeColor','w','LineWidth',0.2); end
    xlabel({'mean brain' strjoin(xStr,'; ')});

    nexttile
    % eventually put std here
    ax = gca; ax.Visible = 'off';

    nexttile
    % eventually put std normalization factor here
    ax = gca; ax.Visible = 'off';

end
tlTitle = title(tl,titleStr,'interpreter','none');


% %% Plot average spectra on loglog plot
% nexttile([1 tl.GridSize(2)])
% plot(f,psdSpecNormAv,'k'); hold on
% ax = gca; ax.YScale = 'log';
% ax.XScale = 'log';
% ax.YAxisLocation = 'right';
% xlabel('Hz'); ylabel('psd averaged within mask');
% grid on
% if exist('fpass','var') && ~isempty(fpass)
%     xlim(fpass)
% else
%     xlim tight
% end
% ylim tight
% 
% %%% Add w references
% f0w = exp(linspace(log(f(2)),log(f(end)),6)); f0w([1 end]) = [];
% yLim = ylim;
% y = yLim(1).*[1 1];
% for ii = 1:length(f0w)
%     x = f0w(ii)+[-1 1].*w;
%     if x(1)>f(1) && x(2)<f(end)
%         plot(x,y,'g','LineWidth',3)
%     end
% end
% 
% %%% Add f0 reference
% if ~isempty(f0)
%     for iii = 1:length(f0)
%         plot([1 1].*f(f0Ind(iii)),ylim,':r');
%         x = f(f0Ind(iii))+[-1 1].*w;
%         yLim = ylim;
%         y = yLim(2).*[1 1];
%         plot(x,y,'g','LineWidth',3)
%     end
% end

%% Plot average spectra on semilog plot
nexttile([1 tl.GridSize(2)])
if exist('psdStruct','var')
    plot(psdStruct.f,psdStruct.psd,'k'); hold on
    plot(psdStruct.f,squeeze(psdStruct.psdErr),'r');
else
    plot(f,psdSpecNormAv,'k'); hold on
end
ax = gca; ax.YScale = 'log';
ax.YAxisLocation = 'right';
xlabel('Hz'); ylabel('psd averaged within mask');
grid on
grid minor
if exist('fpass','var') && ~isempty(fpass)
    xlim(fpass)
else
    xlim tight
end
noiseBase = mean(psdSpecNormAv(f>2.5));
yLim = ylim; yLim(1) = noiseBase; yLim(1) = yLim(1) - 0.1*range(yLim);
ylim(yLim)
% ylim tight

%%% Add noise floor
if strcmp(normMethod,'psdNoise')
    tmpInd = [find(noiseInd,1,'first') find(noiseInd,1,'last')];
    plot(funPsdNorm.psd.f(tmpInd),exp(mean(log(normFact))).*ones(size(tmpInd)),'r')
end

%%% Add f0 (if provided) and w references
if ~isempty(f0)
    for iii = 1:length(f0)
        plot([1 1].*f(f0Ind(iii)),ylim,':r');
        x = f(f0Ind(iii))+[-1 1].*w;
        yLim = ylim;
        y = yLim(2).*[1 1];
        plot(x,y,'-g','LineWidth',3)
    end
    if strcmp(normMethod,'psdNoise')
        legend({'psd' 'noiseFloor' 'f0' '2*w'},'box','off','Location','southeast')
    else
        legend({'psd' 'f0' '2*w'},'box','off','Location','southeast')
    end
else
%     x = f(round(end/2))+[-1 1].*w;
    x = mean(fpass)+[-1 1].*w;
    
    yLim = ylim;
    y = yLim(1).*[1 1];
    plot(x,y,'-g','LineWidth',3)
    if strcmp(normMethod,'psdNoise')
        error('code that')
    else
        legend({'psd' '2*w'},'box','off','Location','northeast')
    end
end
drawnow

%% Launch 3D viewer
if freeviewFlag
    funPsdNormFS = funPsdNorm;
    if isempty(funPsdNormFS.vol)
        funPsdNormFS = vec2vol(funPsdNormFS);
    end
    if ~isempty(f0)
        f0IndExc = true(size(f));
        f0IndExc(f0Ind) = false;
        funPsdNormFS.vol(:,:,:,f0IndExc) = [];
        funPsdNormFS.nframes = nnz(~f0IndExc);
    end
    funPsdNormFS.vol = log(funPsdNormFS.vol);
    funPsdNormFS.vol(isnan(funPsdNormFS.vol)) = 0;
    funTmpName = [tempname '.nii.gz'];
    display('******************')
    display('Writing log(psd) volumes to disk...')
    MRIwrite(funPsdNormFS,funTmpName);

    %%% Print command to view
    anatTmpName = strsplit(funPsdNormFS.fspec,'/'); anatTmpName = strjoin(anatTmpName(1:11),'/');
    cmd = {'freeview'};
%     cmd{end+1} = fullfile(anatTmpName,'T1w_restore.nii.gz');
%     cmd{end+1} = fullfile(anatTmpName,'T1w_restore.2.nii.gz');
%     cmd{end+1} = fullfile(anatTmpName,'T2w_restore.nii.gz');
%     cmd{end+1} = fullfile(anatTmpName,'T2w_restore.2.nii.gz');
%     cmd{end+1} = fullfile(anatTmpName,'Results/rfMRI_REST__brainMean.nii.gz');
%     cmd{end+1} = fullfile(anatTmpName,'Results/rfMRI_REST__brainStd.nii.gz');
    cmd{end+1} = [funTmpName ':colormap=turbo'];

    if isempty(f0)
        display(['To view 3D volume of log(psd) :'])
    else
        display(['To view 3D volume of log(psd) @' num2str(f(f0Ind),'%0.4f') 'Hz :'])
    end
    display(strjoin(cmd,' '))
%     display(['freeview ' funTmpName ':colormap=turbo'])
    if ~isempty(f0)
        display(['recommended window: ' num2str(log(cLim(1))) ' to ' num2str(log(cLim(2)))])
    end
    display('******************')
end