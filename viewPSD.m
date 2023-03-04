function viewPSD(funPsd,f0,mask,id,slc)
if ~exist('f0','var') || isempty(f0)
    f0 = 0.1;
end
if ~exist('id','var') || isempty(id)
    id = '';
elseif isnumeric(id)
    id = num2str(id);
end
if ~exist('mask','var') || isempty(mask)
    mask = funPsd.vol2vec;
end
if ~exist('slc','var') || isempty(slc)
    slc = round(funPsd.depth*0.5);
end
f = funPsd.psd.f;
w = funPsd.psd.w;

%% Apply spatial normalization
if isempty(funPsd.vec)
    funPsd = vol2vec(funPsd);
end
normFact = exp( mean(log(funPsd.vec),1) - mean(log(funPsd.vec(:))));
funPsdNorm = funPsd;
funPsdNorm.vec = funPsd.vec./normFact;

%% Compute summary statistics for later
psdSpecNormAv = exp(mean(log(funPsdNorm.vec),2));
psdSpecNormEr = exp(prctile(log(funPsdNorm.vec),[2.7 97.5],2));
psdSpecNormEr(:,1) = psdSpecNormEr(:,1) - psdSpecNormAv; psdSpecNormEr(:,2) = psdSpecNormAv - psdSpecNormEr(:,2);

%% Plot spatial pattern before and after normalization
funPsd = vec2vol(funPsd);
funPsdNorm = vec2vol(funPsdNorm);
maskC = getMaskOutline(mask(:,:,slc),5);

figure('WindowStyle','docked');
tl = tiledlayout(3,3); tl.TileSpacing = "tight"; tl.Padding = "tight";
tl.TileIndexing = 'rowmajor';
% Title
if exist('id','var') && ischar(id) && ~isempty(id)
    titleStr     = {['subj' id]};
else
    titleStr = {''};
end
titleStr{end+1} = ['tw=' num2str(funPsd.psd.tw)];
titleStr{end+1} = ['w=' num2str(w)];
[~,b,c] = fileparts(funPsd.fspec);
titleStr{end+1} = [b c];
titleStr = [strjoin(titleStr(1:3),'; ') titleStr(4)];
tlTitle = title(tl,titleStr,'interpreter','none');


titleStr = {['f0=' num2str(f0) 'Hz']};
titleStr{end+1} = ['slc' num2str(slc)];
% Before spatial normalization
[~,f0Ind] = min(abs(f-f0));
nexttile
imagesc(squeeze(funPsd.vol(:,:,slc,f0Ind)));
ylabel(colorbar,['psd @ f0=' num2str(f0) 'Hz'])
ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
cLim = ax.CLim;
xlabel({strjoin(titleStr,'; ') 'before spatial norm'})
hold on; plot(maskC,'FaceColor','none')

% After spatial normalization (same scale)
nexttile
imagesc(squeeze(funPsdNorm.vol(:,:,slc,f0Ind)));
ylabel(colorbar,['psd @ f0=' num2str(f0) 'Hz'])
ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
ax.CLim = cLim;
xlabel({strjoin(titleStr,'; ') 'after spatial norm' '(same scale)'})
hold on; plot(maskC,'FaceColor','none')

% After spatial normalization (auto scale)
nexttile
imagesc(squeeze(funPsdNorm.vol(:,:,slc,f0Ind)));
ylabel(colorbar,['psd @ f0=' num2str(f0) 'Hz'])
ax = gca; ax.ColorScale = 'log'; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
xlabel({strjoin(titleStr,'; ') 'after spatial norm' 'auto scale'})
hold on; plot(maskC,'FaceColor','none')


%% Plot average spectra on loglog plot
nexttile(4,[1 3])
plot(f,psdSpecNormAv,'k'); hold on
ax = gca; ax.YScale = 'log';
ax.XScale = 'log';
ax.YAxisLocation = 'right';
xlabel('Hz'); ylabel('psd averaged within mask');
grid on
xlim tight

plot([1 1].*f(f0Ind),ylim,':r');

f0w = exp(linspace(log(f(2)),log(f(end)),6)); f0w([1 end]) = [];
yLim = ylim;
y = yLim(1).*[1 1];
for ii = 1:length(f0w)
    x = f0w(ii)+[-1 1].*w;
    if x(1)>f(1) && x(2)<f(end)
        plot(x,y,'g','LineWidth',3)
    end
end

x = f(f0Ind)+[-1 1].*w;
yLim = ylim;
y = yLim(2).*[1 1];
plot(x,y,'g','LineWidth',3)

%% Plot average spectra on semilog plot
nexttile(7,[1 3])
plot(f,psdSpecNormAv,'k'); hold on
ax = gca; ax.YScale = 'log';
ax.YAxisLocation = 'right';
xlabel('Hz'); ylabel('psd averaged within mask');
grid on
grid minor
xlim tight

plot([1 1].*f(f0Ind),ylim,':r');

x = f(f0Ind)+[-1 1].*w;
yLim = ylim;
y = yLim(2).*[1 1];
plot(x,y,'g','LineWidth',3)
legend({'psd' 'f0' '2*w'},'box','off','Location','southeast')

%% Launch 3D viewer
f0IndExc = true(size(f));
f0IndExc(f0Ind) = false;
funPsdNorm.vol(:,:,:,f0IndExc) = [];
funPsdNorm.vol = log(funPsdNorm.vol);
funPsdNorm.vol(isnan(funPsdNorm.vol)) = 0;
funPsdNorm.nframes = 1;
funTmpName = [tempname '.nii.gz'];
MRIwrite(funPsdNorm,funTmpName);
display('******************')
display(['To view 3D volume of log(psd) @' num2str(f(f0Ind),'%0.3f') 'Hz :'])
display(['freeview ' funTmpName ':colormap=jet'])
display('******************')
