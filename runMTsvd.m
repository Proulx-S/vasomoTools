function svdStruct = runMTsvd(anaType,funTs,fpass,W,mask,normFact)
% Similar to Mitra 1997. A single svd is run on data tapered for
% sensitivity over user-defined frequency band (fpass).
tsMean = mean(funTs.vol,4);
if ~exist('mask','var') || isempty(mask)
    funTs = vol2vec(funTs);
    mask = funTs.vol2vec;
else
    funTs = vol2vec(funTs,mask,1);
end

%% Apply timeseries normalization
% normFact based on psd so we need to use its square root here
if exist('normFact','var')
    funTs.vec = funTs.vec ./ sqrt(normFact(logical(mask))');
else
    error('code that')
end

%% Detrend time series (detrend up to order-2 polynomial, since this is the highest order not fitting a sinwave)
funTs = dtrnd4psd(funTs);
% funTs.vec = funTs.vec - mean(funTs.vec,2);

% %% Zscore time series
% funTs.vec = zscore(funTs.vec,[],1);

%% Set multitaper parameters
tr = funTs.tr/1000;
nFrame = funTs.nframes;
param.tapers = [];
param.Fs = 1/tr;
switch anaType
    case 'svdProulx'
%         error('double-check all that')
        param.anaType = 'proulx';
        T = tr.*funTs.nframes;
%         TW = funTs.nframes/32-1;
%         K = round(TW*2-1);
        K = 10;
        TW = (K+1)/2;
        param.tapers = [TW K];
        if ~isempty(fpass)
            param.fpass = fpass;
        end
        [~,f] = mtspectrumc(funTs.vec(:,1), param);
        %%% Display actual half-widht used
        Wreal = TW/T;
%         display(['w  (halfwidth) requested  : ' num2str(W,'%0.5f ')])
        display(['w  (halfwidth) used       : ' num2str(Wreal,'%0.5f ')])
        display(['tw (time-halfwidth) used  : ' num2str(TW)])
        display(['k  (number of tapers) used: ' num2str(K)])
        mdkp = [];
%         tic
%         [sv,sp,fm] = spsvd(funTs.vec,param,mdkp);
%         % sv:     1 X mode
%         % sp: space X 1    X mode
%         % fm: taper X freq X mode
%         toc
% 
%         figure('WindowStyle','docked');
%         plot(sv)
%         
%         figure('WindowStyle','docked');
%         offsetFac = 1;
%         for kfInd = 1:size(fm,3)
%             y = squeeze(mean(abs(fm(:,:,kfInd)),1));
% %             y = squeeze(exp(mean(log(abs(fm(:,:,kfInd))),1))) + offsetFac*kfInd;
%             plot(f,log(y) + offsetFac*kfInd);
%             text(0.5,mean(y(end-10:end)),num2str(kfInd))
%             hold on
%         end
%         ax = gca;
%         ax.XScale = 'log';
%         %     ax.YScale = 'log';
%         grid on
%         xlim([0.01 0.5])
% 
% 
%         figure('WindowStyle','docked');
%         kfInd = 2;
%         y = abs(fm(:,:,kfInd));
%         plot(f,y,':');
%         hold on
%         plot(f,mean(y,1),'k');
%         ax = gca;
%         ax.YScale = 'log';
%         ax.XScale = 'log';
%         grid on
%         
%         
% 
%         ax = gca;
%         ax.XScale = 'log';
%         %     ax.YScale = 'log';
%         grid on
%         xlim([0.01 0.5])
% 
% 
%         
%         funTs.vec = permute(abs(sp),[3 1 2]);
%         funTs.nframes = size(funTs.vec,1);
%         funTs = vec2vol(funTs);
%         funTmpName = [tempname '.nii.gz'];
%         MRIwrite(funTs,funTmpName);
        


    case 'svdKlein'
        error('double-check all that')
        % if exist('W','var') && ~isempty(W)
        %     anaType = 'svdKlein';
        T = tr.*nFrame;
        TW = T*W;
        K = round(TW*2-1);
        TW = (K+1)/2;
        param.tapers = [TW K];
        if ~isempty(fpass)
            param.fpass = fpass;
        end
        mdkp = [];
        [~,f] = mtspectrumc(funTs.vec(:,1), param);
        %%% Display actual half-widht used
        Wreal = TW/T;
        display(['w  (halfwidth) requested  : ' num2str(W,'%0.5f ')])
        display(['w  (halfwidth) used       : ' num2str(Wreal,'%0.5f ')])
        display(['tw (time-halfwidth) used  : ' num2str(TW)])
        display(['k  (number of tapers) used: ' num2str(K)])
    case 'svdMitra'
        % else
        %     anaType = 'svdMitra';
        %%% Set parameters for the user-defined frequency band
        W = diff(fpass)/2;
        T = tr.*nFrame;
        TW = T*W;
        K = round(TW*2-1);
        TW = (K+1)/2;
        param.tapers = [TW K];
        [~,f] = mtspectrumc(funTs.vec(:,1), param);
        f0 = fpass(1)+W; [~,b] = min(abs(f - f0)); f0 = f(b);
        param.fpass = [f0 f0];
        mdkp = [];
        %%% Display actual frequency band used
        Wreal = W;
        fpassReal = f0+[-1 1].*(TW/T);
        display(['frequency band requested: fpass=[' num2str(fpass,'%0.5f ') ']'])
        display(['frequency band used     : fpass=[' num2str(fpassReal,'%0.5f ') ']'])
    otherwise
        error('Invalid anaType. Choose one of ''svdMitra'' or ''svdKlein''')
end

%% Run the decomposition
tic
[sv,sp,fm,u,s,v,a,proj] = spsvd(funTs.vec,param,mdkp);
toc

% %% Reconstruct reduced psd
% maxModeInd = 1:2;
% A = u(:,maxModeInd)*s(maxModeInd,maxModeInd)*v(:,maxModeInd)';
% sz = size(A);
% sz(3) = sz(2)/param.tapers(2);
% sz(2) = param.tapers(2);
% tmp = reshape(A,sz); % vox X taper X freq
% psdRed = squeeze(mean(conj(tmp).*tmp,2)); % vox x freq
% 
% %% Reconstruct psd at each mode
% nTaper = param.tapers(2);
% nFreq = size(v,1)./nTaper;
% nVox = size(u,1);
% nMode = size(s,2);
% psdRec = nan(nVox,nFreq,nMode);
% for modeInd = 1:size(s,2)
%     A = u(:,modeInd)*s(modeInd,modeInd)*v(:,modeInd)';
%     sz = size(A);
%     sz(3) = sz(2)/nTaper;
%     sz(2) = param.tapers(2);
%     tmp = reshape(A,sz); % vox X taper X freq
%     psdRec(:,:,modeInd) = permute(mean(conj(tmp).*tmp,2),[1 3 2]); % vox x freq x mode
% end

%% Reconstruct psd at each mode (simply from v)
nTaper = param.tapers(2);
tmp = permute(v,[2 1]);
sz = size(tmp);
sz(3) = sz(2)/nTaper;
sz(2) = param.tapers(2);
tmp = reshape(tmp,sz);
psdRecSimple = permute(mean(conj(tmp).*tmp,2),[3 1 2]);

%% Plot reconstructions
[~,pad,~,~,~,~,~]=getparams(param);
N=size(funTs.vec,1);
nfft=max(2^(nextpow2(N)+pad),N);
[f,~]=getfgrid(param.Fs,nfft,param.fpass);

close all
figure('WindowStyle','docked');
plot(diag(s))
ax = gca;
ax.YScale = 'log';
for modInd = 1:30
    figure('WindowStyle','docked');
    tl = tiledlayout(3,2);
    tl.TileSpacing = "tight"; tl.Padding = "tight";
    tl.TileIndexing = 'rowmajor';
    title(tl,['mode ' num2str(modInd)])
    
    tmp = nan(size(mask));
    tmp(logical(mask)) = u(:,modInd);
    
    nexttile(1,[2 1])
    hIm = imagesc(abs(tmp));
    hIm.AlphaData = ~isnan(tmp);
    ax = gca;
    ax.YTick = []; ax.XTick = [];
    ax.ColorScale = 'linear'; ax.Colormap = jet;
    ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ylabel(colorbar,'psd')

    nexttile(2,[2 1])
    hIm = imagesc(angle(tmp));
    hIm.AlphaData = ~isnan(tmp);
    ax = gca;
    ax.YTick = []; ax.XTick = [];
    ax.ColorScale = 'linear'; ax.Colormap = hsv;
    ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ax.CLim = [-pi pi];
    ylabel(colorbar,'phase')

    nexttile(5,[1 2])
%     plot(f,mean(psdRec(:,:,modInd),1))
    plot(f,psdRecSimple(:,modInd))
    ylabel('psd')
    xlabel('Hz')
    drawnow
end


% 
% % [U,S,V] = svd(A)
% % A = u*s*v';
% % A = u(:,1:2)*s(1:2,1:2)*v(:,1:2)';
% % A = s*v';
% A = s(1:2,1:2)*v(:,1:2)';
% sz = size(A);
% sz(3) = sz(2)/param.tapers(2);
% sz(2) = param.tapers(2);
% A = reshape(A,sz); % mode X taper X freq
% 
% % [tapers,pad,Fs,fpass,err,trialave,params]=getparams(param);
% N=size(funTs.vec,1);
% nfft=max(2^(nextpow2(N)+pad),N);
% [f,~]=getfgrid(param.Fs,nfft,param.fpass); 
% 
% S=permute(mean(conj(A).*A,2),[1 3 2]);
% 
% 



%% Output
svdStruct.mask = mask;
svdStruct.normFact = normFact;
svdStruct.tsMean = tsMean;
svdStruct.dim = strjoin({'space/taper' 'freq' 'modes'},' X ');
svdStruct.sv = permute(sv,[3 1 2]);
svdStruct.sp = sp;
svdStruct.fm = fm;
svdStruct.c = sv(:,1)'.^2./sum(sv.^2,2)';
svdStruct.f = f;
svdStruct.w = Wreal;
svdStruct.param = param;


return

%% Plot

% cLim = [min(abs(svdStruct.sp(:))) max(abs(svdStruct.sp(:)))];
cLim = 'auto';

figure('WindowStyle','docked');
tl = tiledlayout(2,2);
tl.Padding = 'tight'; tl.TileSpacing = 'tight';
for ind = 1:4
    nexttile
    tmp = svdStruct.tsMean;
    imagesc(tmp)
    ax = gca;
    ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ax.XTick = []; ax.YTick = [];
    ax.Colormap = gray;
    ylabel(colorbar,'a.u.')
    tmp_cLim = clim; tmp_cLim(2) = tmp_cLim(2)./3;
    clim(tmp_cLim)
end

for ind = 1:size(svdStruct.sp,3)
    figure('WindowStyle','docked');
    tl = tiledlayout(2,2);
    tl.Padding = 'tight'; tl.TileSpacing = 'tight';

    nexttile
    tmp = svdStruct.tsMean;
    imagesc(tmp)
    ax = gca;
    ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ax.XTick = []; ax.YTick = [];
    ax.Colormap = gray;
    ylabel(colorbar,'a.u.')
    tmp_cLim = clim; tmp_cLim(2) = tmp_cLim(2)./3;
    clim(tmp_cLim)

    nexttile
    tmp = nan(size(mask));
    tmp(logical(mask)) = svdStruct.sp(:,:,ind);
    imagesc(abs(tmp))
    ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ax.XTick = []; ax.YTick = [];
    ax.Colormap = jet;
    ylabel(colorbar,'mag')
    clim(cLim)

    nexttile
    imagesc(angle(tmp))
    ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ax.XTick = []; ax.YTick = [];
    ax.Colormap = hsv;
    ax.CLim = [-pi pi];
    ylabel(colorbar,'abs')

    nexttile
    tmp = nan(size(mask));
    tmp(logical(mask)) = svdStruct.sp(:,:,ind);
    imagesc(abs(tmp))
    ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
    ax.XTick = []; ax.YTick = [];
    ax.Colormap = jet;
    ylabel(colorbar,'mag')
    ax.ColorScale = 'log';
    clim(cLim)
end


