function F = viewMTsvdProulx(svdStruct,modePsd,modeTs,modeOrder,funPsd,funTs)
bidsDir = '/autofs/space/takoyaki_001/users/proulxs/vasomo/vasomoAna/vasomoInflow/bids';
if ~exist('modeOrder','var') || isempty(modeOrder)
    modeOrder = 1;
end
if ~exist('K','var')
    K = 1;
end
fileList = cell(size(funTs))';
for runInd = 1:length(funTs)
    tmp = strsplit(funTs(runInd).fspec,filesep);
    fileList{runInd} = tmp{end-1};
    dir(replace(funTs(runInd).fspec,'.nii.gz',''))
    tmp = strsplit(fileList{runInd},'_');
    system(['jq ''.AcquisitionTime'' ' fullfile(bidsDir,tmp{1},tmp{2},'func',[fileList{runInd} '.json'])],'-echo');
end
fileList

%% Compute psd of the reconstructed timeseries
W = [];
K = 1;
cFlag = 1;
if ( ~exist('modePsd','var') || isempty(modePsd) ) && ( exist('modeTs','var') && ~isempty(modeTs) )
    modePsd = runPSD2(modeTs,W,K,cFlag,modeOrder);
end
%% For comparison, compute psd of original timeseries
if ~exist('funPsd','var') || isempty(funPsd)
    funPsd = vec2vol(runPSD2(funTs,W,K,cFlag));
else
    funPsd = vec2vol(funPsd);
end






F = figure('WindowStyle','docked');
tl = tiledlayout(6,5);
tl.TileSpacing = "tight"; tl.Padding = "tight";
tl.TileIndexing = 'rowmajor';

%%% mean brain
tileInd = 1;
nexttile(tileInd,[2 1])
if isfield(funTs,'imMean')
    im = funTs.imMean;
else
    if length(funTs)>1
        im = mean(cat(4,funTs.vol),4);
    else
        im = mean(funTs.vol,4);
    end
end
imagesc(im)
ax = gca;
ax.Colormap = gray;
ax.YTick = []; ax.XTick = [];
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
ax.CLim = prctile(im(:),[0 99]);
title('timeseries mean')

AX = cell(0,0);
AX{end+1} = ax;

%%% Plot psd at peak
vol = cat(6,funPsd.vol);
if ~isreal(vol); vol = conj(vol).*vol; end
for ii = 1:size(svdStruct.param.fpass,1)
    if ii>3
        break
    end
    tileInd = tileInd + 1;
    ax = nexttile(tileInd,[2 1]);

    fpass = svdStruct.param.fpass(ii,:) + [-1 1].*svdStruct.param.BW(ii);
    fInd = funPsd(1).psd.f>fpass(1) & funPsd(1).psd.f<fpass(2);
    hIm = imagesc(mean(mean(vol(:,:,:,fInd,:),5),4));
    % hIm = imagesc(mean(mean(abs(vol(:,:,:,fInd,:)),4),5));
    % hIm = imagesc(mean(mean(abs(vol(:,:,:,fInd,:,:)),4),6));

    ax.Colormap = jet; ax.ColorScale = 'log';
    ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
    ylabel(colorbar,[num2str(mean(fpass),'psd @ %0.2fHz')])

    axPsd{ii} = ax;
end
AX = [AX axPsd];

%%% Plot first singular vector
%%%% Magnitude map
tileInd = 5;
ax = nexttile(tileInd,[2 1]);
spIm = nan(size(funPsd(1).vol2vec));
spIm(funPsd(1).vol2vec) = svdStruct.sp(:,:,modeOrder);
% if ~isreal(spIm); spIm = conj(spIm).*spIm; end
if ~isreal(spIm)
    hIm = imagesc(conj(spIm).*spIm);
else
    hIm = imagesc(spIm);
end
ax.Colormap = jet; %ax.ColorScale = 'log';
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
ax.ColorScale = 'log';
ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
ylabel(colorbar,'first sv weigths mag (a.u.)')
AX{end+1} = ax;
%%%% Phase map
tileInd = 15;
ax = nexttile(tileInd,[2 1]);
hIm = imagesc(angle(spIm),[-pi pi]);
% scale = conj(spIm).*spIm;
scale = abs(spIm);
% scale = log(scale);
scale = scale - min(scale(:)); scale = scale ./ max(scale(:));
hIm.AlphaData = scale;
ax.Colormap = hsv; ax.ColorScale = 'linear';
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
ax.Color = [0.5 0.5 0.5];
ylabel(colorbar,'first sv weigths phase (rad)')
AX{end+1} = ax;
%%%% Phase distribution
tileInd = 14;
n = round(numel(svdStruct.sp(:,:,modeOrder))/500);
binEdges = linspace(-pi,pi,n+1);
binCent = binEdges(1:end-1)+diff(binEdges);
[~,~,bin] = histcounts(angle(svdStruct.sp(:,:,modeOrder)),binEdges);
spBin = nan(1,n);
for i = 1:n
    spBin(i) = mean(svdStruct.sp(bin==i,:,modeOrder));
    % spBin(i) = mean(sp(bin==i));
end
ax = nexttile(tileInd,[2 1]);
[X,Y] = pol2cart(binCent,abs(spBin));
spBin = complex(X,Y);
% polarplot(angle(spBin([1:end 1])),abs(spBin([1:end 1])));
polarplot(angle(spBin([1:end 1])),conj(spBin([1:end 1])).*spBin([1:end 1]));

ax = gca;
ax.ThetaAxisUnits = 'radians';
thetaTick = ax.ThetaTick(1:end-1);
minusPiInd = thetaTick<pi & thetaTick~=0;
plusPiInd = thetaTick>pi & thetaTick~=0;
ax.ThetaTickLabel(plusPiInd) = flip(strcat({'-'},ax.ThetaTickLabel(minusPiInd)));
ax.ThetaAxis.Label.String = 'sv weigth phase';
ax.RAxis.Label.String = 'sv weigth binned mag';

%%% Spectra
%%%% Simple spatial average
tileInd = 21;
ax = nexttile(tileInd,[1 4]);
tmpPsd = catRuns(funPsd);
tmpPsd = vol2vec(tmpPsd);
% tmpPsd = vol2vec(funPsd(1));
f = tmpPsd.psd.f;
y = tmpPsd.vec;
if ~isreal(y)
    y = conj(y).*y;
end
y = mean(y(:,:),2);
% y = mean(y(:,:),2,"omitmissing");
% y = mean(mean(y,4,"omitmissing"),2,"omitmissing");
hP = plot(f,y,'k');
ax.YScale = 'log';
yLim = [min(y(2:end)) max(y(2:end))];
xLim = [min(f) max(f)];
xlim(xLim)
grid on; ax.XMinorGrid = 'on'; ax.YMinorGrid = 'on';
hold on

%%%% weighted avereage of power (no phasing)
scale = svdStruct.sp(:,:,modeOrder);
scale = abs(scale);
scale = scale./sum(scale);
tmp = vol2vec(funPsd);
vec = tmp.vec;
vec = abs(vec);
vec = mean(vec(:,:,:),3);
spec = vec*scale;
hP(end+1) = plot(f,spec,'-');
yLim(1) = min([yLim(1); spec(2:end)]);
yLim(2) = max([yLim(2); spec(2:end)]);


%%%% modePsd
% yyaxis right
hold on
tmpPsd = vol2vec(catRuns(modePsd));
f = tmpPsd.psd.f;
% y = abs(tmpPsd.vec(:,modeOrder,:,:));
y = tmpPsd.vec(:,modeOrder,:,:);
if ~isreal(y)
    % y = abs(y);
    y = conj(y).*y;
end
y = mean(y(:,:,:),3);
% y = abs(tmpPsd.vec(:,modeOrder,:,:));
% y = mean(mean(y,4,"omitmissing"),2,"omitmissing");
hP(end+1) = plot(f,y,'-');
yLim(1) = min([yLim(1); y(2:end)]);
yLim(2) = max([yLim(2); y(2:end)]);

ylim(yLim)
y = [yLim(1) yLim(1) yLim(2) yLim(2) yLim(1)];
if any(diff(svdStruct.param.fpass,[],2)); dbstack; error('X'); end
fRange = svdStruct.param.fpass(:,1) + svdStruct.param.BW.*[-1 1];
for i = 1:size(fRange,1)
    x = [fRange(i,1) fRange(i,2) fRange(i,2) fRange(i,1) fRange(i,1)];
    hPatch = patch(x,y,[1 1 1].*0.7); hPatch.EdgeColor = 'none'; uistack(hPatch,'bottom'); hPatch.FaceAlpha = 0.4;
end


% ax.YScale = 'log';
% yLim = [min(y(2:end)) max(y(2:end))];
% xLim = [min(f) max(f)];
% ylim(yLim)
% xlim(xLim)
% grid on; ax.XMinorGrid = 'on'; ax.YMinorGrid = 'on';

legend(hP,{'unweighted mean' ['mode-' num2str(modeOrder) ' weigthed mean'] ['mode-' num2str(modeOrder) ' filtered']})

% %%%% psd of modeTs
% yyaxis right
% hold on
% tmpPsd = vol2vec(catRuns(runPSD3(modeTs,W,K,modeOrder)));
% f = tmpPsd.psd.f;
% y = tmpPsd.vec;
% if ~isreal(y); y = conj(y).*y; end
% % if ~isreal(y); y = abs(y); end
% hP(2) = plot(f,mean(y(:,:),2),'-');
% ax.YScale = 'linear';
% yLim = [min(y(2:end)) max(y(2:end))];
% xLim = [min(f) max(f)];
% ylim(yLim)
% xlim(xLim)
% grid on; ax.XMinorGrid = 'on'; ax.YMinorGrid = 'on';
% legend(hP,{'unweighted mean' ['mode-' num2str(modeOrder) ' modeTs PSD']})

%%% Eigen spectrum
tileInd = 25;
ax = nexttile(tileInd,[1 1]);
plot(squeeze(svdStruct.sv/sum(svdStruct.sv).*100),'-ok'); hold on
% plot(squeeze(svdStruct.sv),'-ok'); hold on
hPtmp = plot(modeOrder,squeeze(svdStruct.sv(modeOrder)/sum(svdStruct.sv).*100),'-or');
hPtmp.MarkerFaceColor = hPtmp.Color;
% plot(svdStruct.sv,'-ok')
ax.XScale = 'log';
grid on
ylabel('% variance explained')
xlabel('mode order')
% yLim = ylim; yLim(1) = 0; ylim(yLim)


%%% Timeseries
tileInd = 26;
ax = nexttile(tileInd,[1 5]);
t = 0:modeTs(1).tr/1000:modeTs(1).tr/1000*(modeTs(1).nframes*length(modeTs)-1);
% t = 0:svdTs.tr/1000:svdTs.tr/1000.*(svdTs.nframes-1);
y = cat(4,modeTs.vol);
y = y(modeOrder,:);
% y = squeeze(svdTs.vol(modeOrder,:,:,:));
c = colororder;
plot(t,imag(y),'Color',c(1,:)); hold on
plot(t,real(y),'Color',hP(2).Color);
xLim = [min([real(t) imag(t)]) max([real(t) imag(t)])];
if length(modeTs)>1
    x = [1 1] * -modeTs(1).tr/1000/2;
    y = ylim;
    plot(x,y,'k')
    for runInd = 1:length(modeTs)
        x = [1 1] * modeTs(runInd).tr/1000 * (modeTs(runInd).nframes*runInd-0.5);
        plot(x,y,'k')
    end

end
xlim(xLim)

%%% Spectrogram
Sx = cell(1,length(modeTs));
tx = cell(1,length(modeTs));
fx = cell(1,length(modeTs));
TW = (K+1)/2;
for runInd = 1:length(modeTs)
    % y = cat(4,svdTs.vol);
    % y = y(modeOrder,:);
    T = funTs(runInd).tr/1000.*funTs(runInd).nframes;
    W = TW/T;
    param.tapers = [TW K];
    param.Fs = 1/(funTs(runInd).tr/1000);
    winSz = 10;
    stepSz = winSz/4;
    y = squeeze(modeTs(runInd).vol(modeOrder,:,:,:));
    [Sx{runInd},tx{runInd},fx{runInd}]=mtspecgramc(real(y),[winSz stepSz],param);
end
% TW = (K+1)/2;
% T = funTs(1).tr/1000.*funTs(1).nframes;
% W = TW/T;
% param.tapers = [TW K];
% param.Fs = 1/(funTs(1).tr/1000);
% winSz = 10;
% stepSz = winSz/4;
% [Sx,tx,fx]=mtspecgramc(real(y),[winSz stepSz],param);

tileInd = 11;
ax = nexttile(tileInd,[2 3]);
tx2 = tx;
for runInd = 2:length(tx)
    tx2{runInd} = tx2{runInd} + tx2{runInd-1}(end) + modeTs(runInd-1).tr/1000;
end
hIm = imagesc(cat(2,tx2{:}),fx{1},cat(1,Sx{:})');
% runInd = 4;
% hIm = imagesc(tx{runInd},fx{runInd},Sx{runInd}');
ax.YDir = 'normal';
ax.ColorScale = 'log';
ax.Colormap = jet;
xlabel('time (s)')
ylabel('frequency (Hz)')
ylabel(colorbar,'PSD')

hold on
if length(modeTs)>1
    x = [1 1] * -modeTs(1).tr/1000/2;
    y = ylim;
    plot(x,y,'w')
    for runInd = 1:length(modeTs)
        x = [1 1] * modeTs(runInd).tr/1000 * (modeTs(runInd).nframes*runInd-0.5);
        plot(x,y,'w')
    end
end


linkaxes([AX{:}])
