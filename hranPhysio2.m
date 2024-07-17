function [physioTs,reconPhysioTs,reconPhysioTs2] = hranPhysio2(funTs,fRange,nBin)
addpath(genpath('/autofs/space/takoyaki_001/users/proulxs/tools/HRAN'))
if ~exist('nBin','var')
    nBin = 1;
end

K=8;
W=[];
cFlag = 1;
funPsd = runPSD(funTs,W,K,cFlag);

peakFreq = fRange(1) + diff(fRange)/2;
hF = viewPSD2(funPsd,[peakFreq fRange],[],[],funTs.id);

%% Perform mt-svd
anaType = 'svdMitra';
fpass = fRange;
svdStruct = runMTsvd(anaType,funTs,fpass);

sp = svdStruct.sp(:,:,1); % spatial sv (space x 1 x mode)
fm = svdStruct.fm(:,:,1); % taper sv (taper x 1 x mode)
tp = svdStruct.tp; % tapers (taper x time)
tv = svdStruct.tv; % time (1 x time)
f  = svdStruct.f; % freq (1 x freq)

%% Reconstruct component timecourse
i = 1;
tf = tp.*exp(-f(i)*tv);
reconPhysioTs = permute(tf,[2 1])*fm;
funTsX = funTs;
funTsX.vol = permute(real(reconPhysioTs),[2 3 4 1]);
physioReconPsd = runPSD(funTsX,W,K,cFlag);
% whos sp fm tp tv f tf
spaceTaper = sp*permute(fm,[2 1]);
reconPhysioTs2 = spaceTaper*tf;
funTsX = funTs;
funTsX.vol = permute(mean(real(reconPhysioTs2),1),[1 3 4 2]);
physioRecon2Psd = runPSD(funTsX,W,K,cFlag);

% figure('WindowStyle','docked');
% yyaxis left
% plot(real(reconPhysioTs))
% yyaxis right
% plot(mean(real(reconPhysioTs2),1))

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
ax.Color = [0.5 0.5 0.5];


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
polarplot([0 angle(spMean)],ax.RLim,'LineWidth',3,'Color','g')
% polarplot(angle(spMean)+[-pi/2 pi/2],ax.RLim(2).*[1 1],'LineWidth',1,'Color','r')
ax.ThetaAxisUnits = 'radians';

% ts = dtrnd4psd(vol2vec(funTs));
ts = vol2vec(funTs);
ts = ts.vec;


if mod(nBin,2)
    phOffset = angle(spMean);
else
    phOffset = angle(spMean) + 2*pi/nBin/2;
end
ph = wrapToPi(angle(sp) - phOffset);
binEdges = linspace(-pi,pi,nBin+1);
binCent = binEdges(1:end-1)+diff(binEdges);
[~,~,bin] = histcounts(ph,binEdges);
physioTs = nan(size(ts,1),nBin);
physioPsd = cell(1,nBin);
for i = 1:nBin
    curTs = ts(:,bin==i)*abs(sp(bin==i));
    funTsX = funTs;
    funTsX.vol = permute(curTs,[2 3 4 1]);
    physioPsd{i} = runPSD(funTsX,W,K,cFlag);    
    physioTs(:,i) = curTs;
end
binEdges = wrapToPi(binEdges + phOffset);
binCent = wrapToPi(binCent + phOffset);


if nBin>1
    ax = gca;
    theta = cat(1,zeros(1,length(binEdges)-1),binEdges(1:end-1));
    rho = repmat(ax.RLim',[1 length(binEdges)-1]);
    polarplot(theta,rho,'LineWidth',1.5,'Color','r');
end
% polarplot([0 angle(spMean)],ax.RLim,'LineWidth',3,'Color','g')


% ind = tmpPh>-pi/2 & tmpPh<pi/2;
% w = log(abs(sp));
% w = w - min(w(ind));
% w(~ind) = 0;
% w = w ./ norm(w(ind)); %./ max(w(ind));
% tsW1 = ts*w;
% funTs1 = funTs;
% funTs1.vol = permute(tsW1,[2 3 4 1]);
% funPsd1 = runPSD(funTs1,W,K,cFlag);


% tmpPh = wrapToPi(angle(sp) - angle(spMean) + pi);
% ind = tmpPh>-pi/2 & tmpPh<pi/2;
% w = log(abs(sp));
% w = w - min(w(ind));
% w(~ind) = 0;
% w = w ./ norm(w(ind)); %./ max(w(ind));
% tsW2 = ts*w;
% funTs2 = funTs;
% funTs2.vol = permute(tsW2,[2 3 4 1]);
% funPsd2 = runPSD(funTs2,W,K,cFlag);

ax = nexttile(12,[2 1]);
t = 0:funTs.tr/1000:funTs.tr/1000*(funTs.nframes-1);
hPlot = plot(t,physioTs);
ylabel('binned signal')
yyaxis right
hPlot2 = plot(t,real(reconPhysioTs));
ylabel('recon signal')
grid on
xlabel('time (s)')
% plot(tsW1); hold on
% plot(tsW2)

ax = findobj(hF.Children.Children,'type','axes'); tmp = [ax.Title]; ax = ax(ismember({tmp.String},'spectrum average within mask'));
axes(ax);

yyaxis right
for i = 1:length(physioPsd)
    plot(physioPsd{i}.psd.f,physioPsd{i}.vec,'-','Color',hPlot(i).Color); hold on
end
hPlotX = plot(physioReconPsd.psd.f,physioReconPsd.vec./max(physioReconPsd.vec).*max([physioPsd{:}.vec]),'-','Color',hPlot2.Color); hold on
% hPlotX = plot(physioReconPsd.psd.f,physioReconPsd.vec.*30000000,'-','Color',hPlot2.Color); hold on

% figure('WindowStyle','docked');
% yyaxis left
% plot(physioReconPsd.psd.f,physioReconPsd.vec,'-'); hold on
% yyaxis right
% plot(physioRecon2Psd.psd.f,physioRecon2Psd.vec,'-'); hold on

% ylim auto

% plot(funPsd1.psd.f,funPsd1.vec)
% plot(funPsd2.psd.f,funPsd2.vec)
% 
% physioTs = cat(1,funTs1.vol,funTs2.vol);


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

