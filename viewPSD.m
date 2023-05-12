function viewPSD(funPsd,f0,mask,fpass,id,slc,brain)
normMethod = 'psdAv'; % 'smoothPsdAv', 'psdAv' or 'none'
fwhm = 5;
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
end
if ~exist('slc','var') || isempty(slc)
    slc = round(funPsd.depth*0.5);
end
f = funPsd.psd.f;
w = funPsd.psd.w;

%% Perform spatial normalization
switch normMethod
    case 'psdAv'
        if isempty(funPsd.vec)
            funPsd = vol2vec(funPsd);
        end
        normFact0 = exp( mean(log(funPsd.vec),1) );
        normFact = normFact0;
    case 'smoothPsdAv'
        %%% Compute average power and write to disk
        funNorm = funPsd;
        if isempty(funNorm.vol)
            funNorm.vec = mean(log(funNorm.vec),1);
            funNorm.nframes = 1;
            normFact0 = exp(funNorm.vec);
            funNorm = vec2vol(funNorm);
        else
            funNorm.vol = mean(log(funNorm.vol),4);
            funNorm.nframes = 1;
            tmp = vol2vec(funNorm);
            normFact0 = exp(tmp.vec); clear tmp
        end
        funNormTmpName = [tempname '.nii.gz'];
        MRIwrite(funNorm,funNormTmpName);
        

        %%% Define mask for smoothing and write to disk
        funMask = funNorm;
        funMask.vol = ~isnan(funMask.vol);
        funMaskTmpName = [tempname '.nii.gz'];
        MRIwrite(funMask,funMaskTmpName);

        %%% Apply smoothing (mri_fwhm) to normalization map
        cmd = {'source /usr/local/freesurfer/fs-dev-env-autoselect'};
        cmd{end+1} = ['mri_fwhm'...
            ' --i ' funNormTmpName...
            ' --o ' funNormTmpName...
            ' --smooth-only'...
            ' --mask ' funMaskTmpName...
            ' --fwhm ' num2str(fwhm)];
        [status,cmdout] = system(strjoin(cmd,'; '));

        % funNormSm = MRIread(funNormTmpName);
        % figure('WindowStyle','docked');
        % imagesc(exp(funNorm.vol(:,:,46)));
        % ax = gca; ax.ColorScale = 'log';
        % cLim = clim;
        % figure('WindowStyle','docked');
        % imagesc(exp(funNormSm.vol(:,:,46)));
        % ax = gca; ax.ColorScale = 'log';
        % clim(cLim);

        funNorm = MRIread(funNormTmpName);
        funNorm = vol2vec(funNorm,funPsd.vol2vec);
        normFact = exp(funNorm.vec); clear funNorm
    case 'none'
        sz = size(funPsd.vec);
        normFact = ones([1 sz(2)]);
        normFact0 = ones([1 sz(2)]);
    otherwise
        error('')
end

%%% Apply normalization
if isempty(funPsd.vec)
    funPsd = vol2vec(funPsd);
end
funPsdNorm = funPsd;
funPsdNorm.vec = funPsd.vec./normFact .* exp(mean(log(normFact)));

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
figure('WindowStyle','docked');
tl = tiledlayout(3,4); tl.TileSpacing = "tight"; tl.Padding = "tight";
tl.TileIndexing = 'rowmajor';
% Title
if exist('id','var') && ischar(id) && ~isempty(id)
    titleStr     = {['subj' id]};
else
    titleStr = {''};
end
titleStr{end+1} = ['tw=' num2str(funPsd.psd.tw)];
titleStr{end+1} = ['w=' num2str(w) 'Hz'];
switch normMethod
    case 'smoothPsdAv'
        titleStr{end+1} = ['fwhm=' num2str(fwhm) 'mm'];
end
[~,b,c] = fileparts(funPsd.fspec);
titleStr{end+1} = [b c];
titleStr = [strjoin(titleStr(1:end-1),'; ') titleStr(end)];
tlTitle = title(tl,titleStr,'interpreter','none');

if doMaskOutline
    maskC = getMaskOutline(mask(:,:,slc),5);
end

if ~isempty(f0)
    %% Plot spatial pattern before and after normalization
    funPsd = vec2vol(funPsd);
    funPsdNorm = vec2vol(funPsdNorm);
    
    titleStr = {['f0=' num2str(f0) 'Hz']};
    titleStr{end+1} = ['slc' num2str(slc)];
    % Before spatial normalization
    [~,f0Ind] = min(abs(f-f0));
    nexttile
    imagesc(squeeze(funPsd.vol(:,:,slc,f0Ind)));
    ylabel(colorbar,['psd @ f0=' num2str(f0) 'Hz'])
    ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
    ax.Colormap = jet;
    cLim = ax.CLim;
    xlabel({strjoin(titleStr,'; ') 'before spatial norm'})
    hold on;
    if doMaskOutline; plot(maskC,'FaceColor','none'); end

    % Average power
    nexttile
    tmp = funPsdNorm; tmp.vec = normFact0; tmp.nframes = 1; tmp = vec2vol(tmp);
    imagesc(squeeze(tmp.vol(:,:,slc))); clear tmp;
    ylabel(colorbar,'full spectrum psd')
    ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
    ax.Colormap = jet;
    ax.CLim = cLim;
    xlabel('full spectrum psd')
    hold on;
    if doMaskOutline; plot(maskC,'FaceColor','none'); end

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
    if doMaskOutline; plot(maskC,'FaceColor','none'); end

%     % After spatial normalization (same scale)
%     nexttile
%     imagesc(squeeze(funPsdNorm.vol(:,:,slc,f0Ind)));
%     ylabel(colorbar,['psd @ f0=' num2str(f0) 'Hz'])
%     ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
%     ax.CLim = cLim;
%     xlabel({strjoin(titleStr,'; ') 'after spatial norm' '(same scale)'})
%     hold on; plot(maskC,'FaceColor','none')

    % After spatial normalization (auto scale)
    nexttile
    imagesc(squeeze(funPsdNorm.vol(:,:,slc,f0Ind)));
    ylabel(colorbar,['psd @ f0=' num2str(f0) 'Hz'])
    ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
    ax.Colormap = jet;
    xlabel({strjoin(titleStr,'; ') 'after spatial norm'})
    hold on;
    if doMaskOutline; plot(maskC,'FaceColor','none'); end
    cLim = ax.CLim;
else
    %% Plot brain and mask
    titleStr = {['slc' num2str(slc)]};
    
    nexttile
    imagesc(brain.vol(:,:,slc))
    ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
    ax.Colormap = gray;
    xlabel({strjoin(titleStr,'; ') 'after spatial norm' 'auto scale'})
    hold on;
    if doMaskOutline; plot(maskC,'FaceColor','none','EdgeColor','w','LineWidth',0.2); end
    xlabel({'mean brain' strjoin(titleStr,'; ')});

    nexttile
    % eventually put std here
    ax = gca; ax.Visible = 'off';

    nexttile
    % eventually put std normalization factor here
    ax = gca; ax.Visible = 'off';

end

%% Plot average spectra on loglog plot
nexttile(5,[1 4])
plot(f,psdSpecNormAv,'k'); hold on
ax = gca; ax.YScale = 'log';
ax.XScale = 'log';
ax.YAxisLocation = 'right';
xlabel('Hz'); ylabel('psd averaged within mask');
grid on
if exist('fpass','var') && ~isempty(fpass)
    xlim(fpass)
else
    xlim tight
end
ylim tight

%%% Add w references
f0w = exp(linspace(log(f(2)),log(f(end)),6)); f0w([1 end]) = [];
yLim = ylim;
y = yLim(1).*[1 1];
for ii = 1:length(f0w)
    x = f0w(ii)+[-1 1].*w;
    if x(1)>f(1) && x(2)<f(end)
        plot(x,y,'g','LineWidth',3)
    end
end

%%% Add f0 reference
if ~isempty(f0)
    plot([1 1].*f(f0Ind),ylim,':r');
    x = f(f0Ind)+[-1 1].*w;
    yLim = ylim;
    y = yLim(2).*[1 1];
    plot(x,y,'g','LineWidth',3)
end

%% Plot average spectra on semilog plot
nexttile(9,[1 4])
plot(f,psdSpecNormAv,'k'); hold on
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
ylim tight

%%% Add f0 (if provided) and w references
if ~isempty(f0)
    plot([1 1].*f(f0Ind),ylim,':r');
    x = f(f0Ind)+[-1 1].*w;
    yLim = ylim;
    y = yLim(2).*[1 1];
    plot(x,y,'-g','LineWidth',3)
    legend({'psd' 'f0' '2*w'},'box','off','Location','southeast')
else
%     x = f(round(end/2))+[-1 1].*w;
    x = mean(fpass)+[-1 1].*w;
    
    yLim = ylim;
    y = yLim(1).*[1 1];
    plot(x,y,'-g','LineWidth',3)
    legend({'psd' '2*w'},'box','off','Location','northeast')
end
drawnow

%% Launch 3D viewer
if freeviewFlag
    if isempty(funPsdNorm.vol)
        funPsdNorm = vec2vol(funPsdNorm);
    end
    if ~isempty(f0)
        f0IndExc = true(size(f));
        f0IndExc(f0Ind) = false;
        funPsdNorm.vol(:,:,:,f0IndExc) = [];
        funPsdNorm.nframes = nnz(~f0IndExc);
    end
    funPsdNorm.vol = log(funPsdNorm.vol);
    funPsdNorm.vol(isnan(funPsdNorm.vol)) = 0;
    funTmpName = [tempname '.nii.gz'];
    display('******************')
    display('Writing log(psd) volumes to disk...')
    MRIwrite(funPsdNorm,funTmpName);

    %%% Print command to view
    anatTmpName = strsplit(funPsdNorm.fspec,'/'); anatTmpName = strjoin(anatTmpName(1:11),'/');
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