function physioTs = hranPhysio(funTs)
addpath(genpath('/autofs/space/takoyaki_001/users/proulxs/tools/HRAN'))

K=8;
W=[];
cFlag = 1;
funPsd = runPSD(funTs,W,K,cFlag);

peakFreq = 0.338;
peakFreqRange = peakFreq+[-1 1].*0.065;
hF = viewPSD2(funPsd,[peakFreq peakFreqRange],[],[],funTs.id);

%% Extract physio spatial filter
anaType = 'svdMitra';
fpass = peakFreqRange;
svdStruct = runMTsvd(anaType,funTs,fpass);

sp = svdStruct.sp(:,:,1); % spatial sv (space x 1 x mode)
fm = svdStruct.fm(:,:,1); % taper sv (taper x 1 x mode)
tp = svdStruct.tp; % tapers (taper x time)
tv = svdStruct.tv; % time (1 x time)
f  = svdStruct.f; % freq (1 x freq)
i = 1;
tf = tp.*exp(-f(i)*tv);
% whos sp fm tp tv f tf
spaceTaper = sp*permute(fm,[2 1]);
data2 = spaceTaper*tf;

spIm = zeros(size(funPsd.vol2vec));
spIm(funPsd.vol2vec) = sp;

figure(hF)
ax = nexttile(2,[2 1]); ax.Visible = 'on';
imagesc(abs(spIm));
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
ax.Colormap = jet; ax.ColorScale = 'log';
ylabel(colorbar,'sv weights');

ax = nexttile(3,[2 1]); ax.Visible = 'on';
hIm = imagesc(angle(spIm));
ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
ax.Colormap = hsv; ax.ColorScale = 'linear';
trans = abs(spIm);
trans = trans - min(trans(:));
trans = trans./max(trans(:));
hIm.AlphaData = trans;
ax.CLim = [-pi pi];
ylabel(colorbar,'sv weights phase (rad)');



n = 100;
binEdges = linspace(-pi,pi,n+1);
binCent = binEdges(1:end-1)+diff(binEdges);
[~,~,bin] = histcounts(angle(sp),binEdges);
spBin = nan(1,n);
for i = 1:n
    spBin(i) = mean(sp(bin==i));
end
ax = nexttile(11,[2 1]);
[X,Y] = pol2cart(binCent,abs(spBin)); spBin = complex(X,Y);
polarplot(angle(spBin),abs(spBin));
hold on
spMean = mean(sp);
ax = gca;
polarplot([0 angle(spMean)],ax.RLim,'LineWidth',3,'Color','r')
polarplot(angle(spMean)+[-pi/2 pi/2],ax.RLim(2).*[1 1],'LineWidth',1,'Color','r')
ax.ThetaAxisUnits = 'radians';

% ts = dtrnd4psd(vol2vec(funTs));
ts = vol2vec(funTs);
ts = ts.vec;


tmpPh = wrapToPi(angle(sp) - angle(spMean));
ind = tmpPh>-pi/2 & tmpPh<pi/2;
w = log(abs(sp));
w = w - min(w(ind));
w(~ind) = 0;
w = w ./ norm(w(ind)); %./ max(w(ind));
tsW1 = ts*w;
funTs1 = funTs;
funTs1.vol = permute(tsW1,[2 3 4 1]);
funPsd1 = runPSD(funTs1,W,K,cFlag);


tmpPh = wrapToPi(angle(sp) - angle(spMean) + pi);
ind = tmpPh>-pi/2 & tmpPh<pi/2;
w = log(abs(sp));
w = w - min(w(ind));
w(~ind) = 0;
w = w ./ norm(w(ind)); %./ max(w(ind));
tsW2 = ts*w;
funTs2 = funTs;
funTs2.vol = permute(tsW2,[2 3 4 1]);
funPsd2 = runPSD(funTs2,W,K,cFlag);

ax = nexttile(12,[2 1]);
plot(tsW1); hold on
plot(tsW2)

ax = findobj(hF.Children.Children,'type','axes'); tmp = [ax.Title]; ax = ax(ismember({tmp.String},'spectrum average within mask'));
axes(ax);

plot(funPsd1.psd.f,funPsd1.vec)
plot(funPsd2.psd.f,funPsd2.vec)

physioTs = cat(1,funTs1.vol,funTs2.vol);


return

% 
% 
% angle(mean(sp(ind)))
% 
% ind = angle(sp) > angle(spMean)-pi/2;
% physTrace1 = mean(ts(:,ind),2);
% ind = angle(sp) <= angle(spMean)-pi/2;
% physTrace1 = mean(ts(:,ind),2);
% 
% 
% %% Visualize binned timeseries
% magThresh = exp(-4);
% figure('WindowStyle','docked');
% hT = tiledlayout(2,2); hT.TileSpacing = 'tight'; hT.Padding = 'tight';
% nexttile
% histogram(log(abs(spIm(:)))); hold on
% plot([1 1].*log(magThresh),ylim,'r')
% nexttile
% histogram(angle(spIm(abs(spIm)>magThresh))); hold on
% nexttile
% hIm = imagesc(angle(spIm));
% ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
% ax.Colormap = hsv; ax.ColorScale = 'linear';
% hIm.AlphaData = abs(spIm)>magThresh;
% colorbar
% ax.CLim = [-pi pi];
% 
% mask = abs(spIm)>magThresh;
% ts = vec2vol(dtrnd4psd(vol2vec(funTs)));
% ts = permute(ts.vol,[4 1 2 3]);
% ts = ts(:,mask);
% cmpl = spIm(mask)';
% mag = abs(spIm(mask))';
% ph = angle(spIm(mask))';
% ts = ts./mag;
% [~,b] = sort(ph);
% nexttile
% hIm = imagesc(ts(:,b));
% ax = gca; %ax.CLim = [0.5 1.5];
% % ax.Colormap = jet;
% colorbar
% tmp = mag - min(mag);
% tmp = tmp./max(tmp);
% hIm.AlphaData = repmat(tmp(b),[size(ts,1) 1]);
% 
% plot(mean(ts(:,ph<0),2)); hold on
% plot(mean(ts(:,ph>0),2)); hold on
% 
% n = 100;
% binEdges = linspace(-pi,pi,n+1);
% binCent = binEdges(1:end-1)+diff(binEdges);
% [~,~,bin] = histcounts(ph,binEdges);
% physGraph = nan(size(ts,1),n);
% physGraphAlpha = nan(size(ts,1),n);
% cmplBin = nan(1,n);
% for i = 1:n
%     physGraph(:,i) = mean(ts(:,bin==i),2);
%     physGraphAlpha(:,i) = abs(mean(cmpl(bin==i),2));
%     cmplBin(i) = mean(cmpl(bin==i),2);
% end
% 
% figure('WindowStyle','docked');
% hIm = imagesc(physGraph');
% physGraphAlpha = physGraphAlpha - min(physGraphAlpha(:));
% physGraphAlpha = physGraphAlpha ./ max(physGraphAlpha(:));
% hIm.AlphaData = physGraphAlpha';
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Extract physio spec
% spec = funPsd.psd.spec;
% funPsd = vec2vol(funPsd);
% spec = permute(funPsd.vol,[4 1 2 3]);
% spec = spec(:,:);
% physSpec = spec*physMap(:);
% 
% figure(hF);
% ax = hF.Children.Children(end-5);
% axes(ax);
% yyaxis left
% plot(funPsd.psd.f,conj(physSpec).*physSpec)
% % ax
% % 
% % size(tmp(1,:)*physMap(:));
% % size(funPsd.vol(:,:,1,1)*physMap);
% 
% 
% % nexttile(10,[2 1])
% % imagesc(trans);
% % ax = gca;
% % ax.DataAspectRatio = [1 1 1];
% % ax.PlotBoxAspectRatio = [1 1 1];
% % ax.Colormap = jet;
% % ax.ColorScale = 'log';
% 
% 
% 
% 
% 


% HRAN_demo_nifti
% prpFileList %','prpAvFile','oriFileList','oriAvFile');
% for runInd = 1:length(prpFileList)
%     funTs = MRIread(prpFileList{runInd});
% 
% 
% 
% end

physiologicalData = squeeze(resp(2,:,:,:)); % time x 1: a single fMRI timeseries that should contain physiological signals

inputParams = struct;
% 2) Set TR, moving window length, percent overlap
TR = funTs.tr/1000; %0.347s in the demo data
inputParams.TR = TR; %TR of fMRI data (s)
inputParams.windowLength = 30; % length of moving window (s) (e.g. 24-30s)
inputParams.percentOverlap = .75; % percent overlap of windows (e.g. .5 or .75)

% 3) Set neural regressors of design matrix (time x # regressors)
% In general, one can convolve the neural stimulus with HRF and input into
% design matrix (or use whatever design matrix they wish)

% Here, we use a cos and sin at the particular neural frequency (which
% approximate the convolution at steady state)
time = [0:inputParams.TR:funTs.nframes*inputParams.TR-inputParams.TR];
% time = [0:inputParams.TR:size(V,4)*inputParams.TR-inputParams.TR];
neuralFrequency = .13; % Hz
neuralSignal = 60*cos(2*pi*neuralFrequency.*time);
neuralZ = [cos(2*pi*neuralFrequency.*time)';sin(2*pi*neuralFrequency.*time)'];
inputParams.neuralZ = neuralZ;

% 4) Set physiological frequency range
inputParams.cardiacFreqRange = [50:80]; % cardiac fundamental freq range (bpm)
inputParams.respFreqRange = [6:27]; %res    p fundamental freq range (bpm)

% 5) Set model orders used to estimate the cardiac and respiratory
% frequencies from the physiological time series
inputParams.P_freq = 2; % AR order used to estimate physio frequencies (e.g. 1-4)
inputParams.R_freq = 1; % Resp order used to estimate physio frequencies (e.g. 1-4)
inputParams.C_freq = 1; % Cardiac order used to estimate physio frequencies (e.g. 1-3)
inputParams.N_freq = 0; % Number of neural regressors used to estimate physio freqs (e.g. 0 bc assuming no neural signal in this region)
inputParams.X_freq = 0; % Interaction order used to estimate physio frequencies (e.g. 0-1)

% 5) Set model orders used to create design matrix to regress physio noise
% from each voxel in the brain
inputParams.P_Z = 2; % AR order used to regress physio noise (e.g. 1-4)
inputParams.R_Z = 2; % Resp order used to regress physio noise (e.g. 1-4)
inputParams.C_Z = 3; % Cardiac order used to regress physio noise (e.g. 1-3)
inputParams.N_Z = size(neuralZ,2); % Number of neural regressors used to regress physio noise (e.g. number of columns in neuralZ)
inputParams.X_Z = 0; % Interaction order used to regress physio noise (e.g. 0-1)

% Toggle for waitbar on or off
waitbarBoolean = 1; % 1 - on, 0 - off
HRAN_estFreqs_outputParams = HRAN_estimatePhysFreqs(physiologicalData,inputParams,waitbarBoolean);


fig = figure('Position',[1 1 800 600]);

% Spectrogram of physiological data
subplot(1,2,1)
movingwin = [24 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[spec_original, stime, sfreq] = mtspecgramc(detrend(physiologicalData), movingwin, params);
hold on
imagesc(stime, sfreq, 10*log10(spec_original)')
scatter(HRAN_estFreqs_outputParams.t,HRAN_estFreqs_outputParams.w_hr_hat,'w','filled')
scatter(HRAN_estFreqs_outputParams.t,HRAN_estFreqs_outputParams.w_rr_hat,'w','filled')
hold off
set(gca,'YDir','normal')
colormap('jet')
colorbar
xlim([stime(1) stime(end)])
ylim([sfreq(1) sfreq(end)])
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('Physiological Noise ROI')
