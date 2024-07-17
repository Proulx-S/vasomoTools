function [svdTs,svdStruct] = hranPhysio3(funTs,fRange,modeOrder,K)
% addpath(genpath('/autofs/space/takoyaki_001/users/proulxs/tools/HRAN'))
if ~exist('K','var'); K = 1; end
if ~exist('modeOrder','var'); modeOrder = 1; end

% W=[];
% cFlag = 0;
% funPsd = runPSD(funTs,W,K,cFlag);
% 
% peakFreq = fRange(:,1) + diff(fRange,[],2)/2;
% hF = viewPSD2(funPsd,[peakFreq(:)' fRange(:)'],[],[],funTs.id);

%% Perform MTsvd
anaType = 'svdProulx';
svdStruct = runMTsvd2(anaType,funTs,fRange);

%% Reconstruct timeseries of the first component
sp = svdStruct.sp(:,:,modeOrder); % spatial sv (space x 1 x mode)
fm = svdStruct.fm(:,:,modeOrder); % taper sv (taper x 1 x mode)
tp = svdStruct.MTS.proj; % tapers (taper x time)
f = svdStruct.MTS.bandFreq; % frequencies (1 x freq)
fInd = svdStruct.MTS.bandInd; % frequency indices (taper x 1)
% whos sp fm tp f fInd
tsRec = fm'*tp; % WARNING! taking the conjugate here with the ' operator

%% Inspect original and reconstructed timeseries
W=[];
cFlag = 1;
funPsd = runPSD(funTs,W,K,cFlag); funPsd = vec2vol(funPsd);
% svdTs = funTs; svdTs.vol = permute(conj(tsRec).*tsRec,[1 3 4 2]);
svdTs = funTs; svdTs.vol = permute(tsRec,[1 3 4 2]);
svdTs.volsize([1 2]) = 1; svdTs.height = 1; svdTs.width = 1; svdTs.nvoxels = 1;
svdPsd = runPSD(svdTs,W,K,cFlag);


%% Plot
hF = figure('WindowStyle','docked');
tl = tiledlayout(6,5);
tl.TileSpacing = "tight"; tl.Padding = "tight";
tl.TileIndexing = 'rowmajor';

%%% mean brain
tileInd = 1;
nexttile(tileInd,[2 1])
if isfield(funPsd,'tMean') && ~all(size(funPsd.tMean)==1)
    im = squeeze(funPsd.tMean);
    imagesc(im)
    ax = gca;
    ax.Colormap = gray;
    ax.YTick = []; ax.XTick = [];
    ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ax.CLim = prctile(im(:),[0 99]);
    title('timeseries mean')
else
    ax = gca;
    ax.Visible = 'off';
end

%%% Plot psd at peak
for ii = 1:size(svdStruct.param.fpass,1)
    if ii>3
        break
    end
    tileInd = tileInd + 1;
    ax = nexttile(tileInd,[2 1]);

    fpass = svdStruct.param.fpass(ii,:) + [-1 1].*svdStruct.param.BW(ii);
    fInd = funPsd.psd.f>fpass(1) & funPsd.psd.f<fpass(2);
    
    if cFlag
        hIm = imagesc(mean(abs(funPsd.vol(:,:,:,fInd)),4));
    else
        hIm = imagesc(mean(funPsd.vol(:,:,:,fInd),4));
    end
    ax.Colormap = jet; ax.ColorScale = 'log';
    ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
    ylabel(colorbar,[num2str(mean(fpass),'psd @ %0.2fHz')])

    axPsd{ii} = ax;
end

%%% Plot first singular vector
%%%% Magnitude map
tileInd = 5;
ax = nexttile(tileInd,[2 1]);
spIm = nan(size(funPsd.vol2vec));
spIm(funPsd.vol2vec) = sp;
hIm = imagesc(abs(spIm));
ax.Colormap = jet; ax.ColorScale = 'log';
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
ylabel(colorbar,'first sv weigths mag (a.u.)')
%%%% Phase map
tileInd = 15;
ax = nexttile(tileInd,[2 1]);
hIm = imagesc(angle(spIm),[-pi pi]);
scale = abs(spIm); scale = scale - min(scale(:)); scale = scale ./ max(scale(:));
hIm.AlphaData = scale;
ax.Colormap = hsv; ax.ColorScale = 'linear';
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
ax.Color = [0.5 0.5 0.5];
ylabel(colorbar,'first sv weigths phase (rad)')
%%%% Phase distribution
tileInd = 14;
n = round(numel(sp)/500);
binEdges = linspace(-pi,pi,n+1);
binCent = binEdges(1:end-1)+diff(binEdges);
[~,~,bin] = histcounts(angle(sp),binEdges);
spBin = nan(1,n);
for i = 1:n
    spBin(i) = mean(sp(bin==i));
end
ax = nexttile(tileInd,[2 1]);
[X,Y] = pol2cart(binCent,abs(spBin)); spBin = complex(X,Y);
polarplot(angle(spBin([1:end 1])),abs(spBin([1:end 1])));
ax = gca;
ax.ThetaAxisUnits = 'radians';
thetaTick = ax.ThetaTick(1:end-1);
minusPiInd = thetaTick<pi & thetaTick~=0;
plusPiInd = thetaTick>pi & thetaTick~=0;
ax.ThetaTickLabel(plusPiInd) = flip(strcat({'-'},ax.ThetaTickLabel(minusPiInd)));
ax.ThetaAxis.Label.String = 'sv weigth phase';
ax.RAxis.Label.String = 'sv weigth binned mag';

%%% Spectra
%%%% Spatial average
tileInd = 21;
ax = nexttile(tileInd,[1 4]);
tmpPsd = vol2vec(funPsd);
f = tmpPsd.psd.f;
y = mean(tmpPsd.vec,2,"omitmissing");
hP = plot(f,y,'k');
ax.YScale = 'log';
yLim = [min(y(2:end)) max(y(2:end))];
xLim = [min(f) max(f)];
ylim(yLim)
xlim(xLim)
grid on; ax.XMinorGrid = 'on'; ax.YMinorGrid = 'on';

y = [yLim(1) yLim(1) yLim(2) yLim(2) yLim(1)];
for i = 1:size(fRange,1)
    x = [fRange(i,1) fRange(i,2) fRange(i,2) fRange(i,1) fRange(i,1)];
    hPatch = patch(x,y,[1 1 1].*0.7); hPatch.EdgeColor = 'none'; uistack(hPatch,'bottom'); hPatch.FaceAlpha = 0.4;
end

% %%%% weighted avereage
% scale = abs(sp); scale = scale./norm(scale);
% spec = funPsd.psd.spec.*scale;
% spec = conj(spec).*spec;
% spec = mean(spec,3);
% spec = sum(spec,1);
% yyaxis right
% plot(f,spec)
%%%% Of reconstructed timeseries
yyaxis right
hold on
tmpPsd = vol2vec(svdPsd);
f = tmpPsd.psd.f;
y = mean(tmpPsd.vec,2,"omitmissing");
hP(2) = plot(f,y,'-');
ax.YScale = 'linear';
yLim = [min(y(2:end)) max(y(2:end))];
xLim = [min(f) max(f)];
ylim(yLim)
xlim(xLim)
grid on; ax.XMinorGrid = 'on'; ax.YMinorGrid = 'on';

legend(hP,{'unweighted mean' 'svd recon'})

%%% Eigen spectrum
tileInd = 25;
ax = nexttile(tileInd,[1 1]);
plot(svdStruct.sv/sum(svdStruct.sv).*100,'-ok'); hold on
hPtmp = plot(modeOrder,svdStruct.sv(modeOrder)/sum(svdStruct.sv).*100,'-or');
hPtmp.MarkerFaceColor = hPtmp.Color;
% plot(svdStruct.sv,'-ok')
ax.XScale = 'log';
grid on
ylabel('% variance explained')
xlabel('mode order')


%%% Timeseries
tileInd = 26;
ax = nexttile(tileInd,[1 5]);
t = 0:svdTs.tr/1000:svdTs.tr/1000.*(svdTs.nframes-1);
y = squeeze(svdTs.vol);
plot(t,real(y)); hold on
plot(t,imag(y))
xLim = [min([real(t) imag(t)]) max([real(t) imag(t)])];
xlim(xLim)

%%% Spectrogram
% K = 1;
TW = (K+1)/2;
T = funTs.tr/1000.*funTs.nframes;
W = TW/T;
param.tapers = [TW K];
param.Fs = 1/(funTs.tr/1000);
winSz = 10;
stepSz = winSz/4;
[Sx,tx,fx]=mtspecgramc(real(y),[winSz stepSz],param);

tileInd = 11;
ax = nexttile(tileInd,[2 3]);
hIm = imagesc(tx,fx,Sx');
ax.YDir = 'normal';
ax.ColorScale = 'log';
ax.Colormap = jet;
xlabel('time (s)')
ylabel('frequency (Hz)')
ylabel(colorbar,'PSD')



