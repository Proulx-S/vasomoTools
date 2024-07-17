function svdStruct = runMTsvd(anaType,funTs,fpass,W,K,mask,vecNorm)
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
% normFact (vector) based on psd so we need to use its square root here
if exist('normFact','var') && exist('vecNorm','var') && ~isempty(vecNorm)
    funTs.vec = funTs.vec ./ sqrt(vecNorm(logical(mask(funTs.vol2vec))));
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
%     case 'svdProulx'
% %         error('double-check all that')
%         param.anaType = 'proulx';
%         T = tr.*funTs.nframes;
% %         TW = funTs.nframes/32-1;
% %         K = round(TW*2-1);
%         K = 10;
%         TW = (K+1)/2;
%         param.tapers = [TW K];
%         if ~isempty(fpass)
%             param.fpass = fpass;
%         end
%         [~,f] = mtspectrumc(funTs.vec(:,1), param);
%         %%% Display actual half-widht used
%         Wreal = TW/T;
% %         display(['w  (halfwidth) requested  : ' num2str(W,'%0.5f ')])
%         display(['w  (halfwidth) used       : ' num2str(Wreal,'%0.5f ')])
%         display(['tw (time-halfwidth) used  : ' num2str(TW)])
%         display(['k  (number of tapers) used: ' num2str(K)])
%         mdkp = [];
% %         tic
% %         [sv,sp,fm] = spsvd(funTs.vec,param,mdkp);
% %         % sv:     1 X mode
% %         % sp: space X 1    X mode
% %         % fm: taper X freq X mode
% %         toc
% % 
% %         figure('WindowStyle','docked');
% %         plot(sv)
% %         
% %         figure('WindowStyle','docked');
% %         offsetFac = 1;
% %         for kfInd = 1:size(fm,3)
% %             y = squeeze(mean(abs(fm(:,:,kfInd)),1));
% % %             y = squeeze(exp(mean(log(abs(fm(:,:,kfInd))),1))) + offsetFac*kfInd;
% %             plot(f,log(y) + offsetFac*kfInd);
% %             text(0.5,mean(y(end-10:end)),num2str(kfInd))
% %             hold on
% %         end
% %         ax = gca;
% %         ax.XScale = 'log';
% %         %     ax.YScale = 'log';
% %         grid on
% %         xlim([0.01 0.5])
% % 
% % 
% %         figure('WindowStyle','docked');
% %         kfInd = 2;
% %         y = abs(fm(:,:,kfInd));
% %         plot(f,y,':');
% %         hold on
% %         plot(f,mean(y,1),'k');
% %         ax = gca;
% %         ax.YScale = 'log';
% %         ax.XScale = 'log';
% %         grid on
% %         
% %         
% % 
% %         ax = gca;
% %         ax.XScale = 'log';
% %         %     ax.YScale = 'log';
% %         grid on
% %         xlim([0.01 0.5])
% % 
% % 
% %         
% %         funTs.vec = permute(abs(sp),[3 1 2]);
% %         funTs.nframes = size(funTs.vec,1);
% %         funTs = vec2vol(funTs);
% %         funTmpName = [tempname '.nii.gz'];
% %         MRIwrite(funTs,funTmpName);
        


    case 'svdKlein'
% %         error('double-check all that')
%         % if exist('W','var') && ~isempty(W)
%         %     anaType = 'svdKlein';
%         T = tr.*nFrame;
%         TW = T*W;
%         K = round(TW*2-1);
%         TW = (K+1)/2;
%         param.tapers = [TW K];
%         if ~isempty(fpass)
%             param.fpass = fpass;
%         end
%         mdkp = [];
%         [~,f] = mtspectrumc(funTs.vec(:,1), param);
%         %%% Display actual half-widht used
%         Wreal = TW/T;
%         display(['w  (halfwidth) requested  : ' num2str(W,'%0.5f ')])
%         display(['w  (halfwidth) used       : ' num2str(Wreal,'%0.5f ')])
%         display(['tw (time-halfwidth) used  : ' num2str(TW)])
%         display(['k  (number of tapers) used: ' num2str(K)])




        Wflag = exist('W','var') && ~isempty(W);
        Kflag = exist('K','var') && ~isempty(K);
        if Wflag && Kflag
            error('Cannot specify both W and K');
        elseif Wflag
            T = tr.*funTs.nframes;
            TW = T*W;
            K = round(TW*2-1);
            TW = (K+1)/2;
            param.tapers = [TW K];
            Wreal = TW/T;
            display(['w  (halfwidth) requested  : ' num2str(W,'%0.5f ')])
            display(['w  (halfwidth) used       : ' num2str(Wreal,'%0.5f ')])
            display(['tw (time-halfwidth) used  : ' num2str(TW)])
            display(['k  (number of tapers) used: ' num2str(K)])
        elseif Kflag
            TW = (K+1)/2;
            T = tr.*funTs.nframes;
            W = TW/T; Wreal = W;
            param.tapers = [TW K];
            display(['k  (number of tapers) requested  : ' num2str(K)])
            display(['w  (halfwidth) used       : ' num2str(W,'%0.5f ')])
            display(['tw (time-halfwidth) used  : ' num2str(TW)])
        end
        if ~isempty(fpass)
            param.fpass = fpass;
        end
        mdkp = [];
        [~,f] = mtspectrumc(funTs.vec(:,1), param);

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
    case 'svdProulx'
        % warning('code that')
        % keyboard
        % else
        %     anaType = 'svdMitra';
        %%% Set parameters for each user-defined frequency band
        fpassOrig = fpass;
        paramOrig = param;
        for bandInd = 1:size(fpassOrig,1)
            fpass_targ = fpassOrig(bandInd,:);
            W_targ = diff(fpass_targ)/2;
            T = tr.*nFrame;
            TW_targ = T*W_targ;
            K_actual = round(TW_targ*2-1);
            TW_actual = (K_actual+1)/2;
            W_actual = TW_actual/T;
            
            paramCur = paramOrig;
            paramCur.tapers = [TW_actual K_actual];
            [~,f] = mtspectrumc(funTs.vec(:,1), paramCur);
            
            f0_targ = mean(fpass_targ); [~,b] = min(abs(f - f0_targ));
            f0_actual = f(b);
            fpass_actual = f0_actual+[-1 1].*W_actual;

            paramCur.fpass = [f0_actual f0_actual];
            mdkp = [];
            
            % display(['frequency band requested: fpass=[' num2str(fpass,'%0.5f ') ']'])
            % display(['frequency band used     : fpass=[' num2str(fpassReal,'%0.5f ') ']'])

            param.tapers(bandInd,:) = paramCur.tapers;
            param.fpass(bandInd,:) = paramCur.fpass;
            param.BW(bandInd,1) = W_actual;
        end
    otherwise
        error('Invalid anaType. Choose one of ''svdMitra'', ''svdKlein'' or ''svdProulx''')
end

%% Run the decomposition
tic
% param.fpass = fpass;
% [u,s,v,f,bandV] = spsvd2(funTs,param); sp = []; sv = []; fm = [];
% [sv,sp,fm,u,s,v,a,proj] = spsvd(funTs.vec,param,mdkp);
% [sv,sp,fm] = spsvd(funTs.vec,param,mdkp);
switch anaType
    case 'svdProulx'
        [sv,sp,fm,tapers,tvec,f,data] = spsvd3(funTs.vec,param,mdkp);
    otherwise
        [sv,sp,fm,tapers,tvec,f,data] = spsvd(funTs.vec,param,mdkp);
end
toc

% i=1;
% timeTaper=tapers.*exp(-f(i)*tvec);
% spaceTaper = sp(:,:,1)*permute(fm(:,:,1),[2 1]);
% time = spaceTaper(1,:)*permute(timeTaper,[2 1]);
% 
% t = 0:funTs.tr/1000:funTs.tr/1000*(size(timeTaper,1)-1);
% plot(t,real(time)); hold on
% 
% data2 = (data'*timeTaper)*permute(timeTaper,[2 1]);


if 0
%% Get component timecourse
[~,fInd] = min(abs(f-0.09343));
cInd = 1;
sX = diag(sv(fInd,cInd));
vX = conj(permute(fm(:,fInd,cInd),[1 3 2]));
tm = proj*(vX*sX);

%% Reconstructed data (best guess)
[~,fInd] = min(abs(f-0.09343));
cInd = 1;
uX = conj(permute(sp(:,fInd,cInd),[1 3 2]));
sX = diag(sv(fInd,cInd));
vX = conj(permute(fm(:,fInd,cInd),[1 3 2]));
aX = uX*sX*vX';
fvec = exp(-f(fInd)*tvec);
proj=tapers.*fvec;
dataRec = proj*aX';
whos aX proj dataRec data


%% Reconstruct reduced psd
[~,fInd] = min(abs(f-0.09343));
cInd = 1;
whos sv sp fm tapers tvec f data
sp(:,fInd,cInd);
sv(fInd,cInd);
fm(:,fInd,cInd);

uX = conj(permute(sp(:,fInd,cInd),[1 3 2]));
sX = diag(sv(fInd,cInd));
vX = conj(permute(fm(:,fInd,cInd),[1 3 2]));
whos uX sX vX

tmpX = proj*(vX*sX);
whos proj vX tmpX
cIndX = cInd;
% plot(real(tmpX(:,cIndX))); hold on
% plot(imag(tmpX(:,cIndX))); hold on
plot(abs(tmpX(:,cIndX)))
plot(abs(tmpX(:,1)),'k','LineWidth',8)
sum(abs(tmpX),1)


fvec = exp(-f(fInd)*tvec);
proj=tapers.*fvec;
a=data'*proj; % projected data
[u,s,v]= svd(a,0); % svd
dataX = (a/proj)';
imagesc(data); colorbar
imagesc(real(dataX)); colorbar

voxInd = 1;
plot(data(:,voxInd))
tmpX = a(voxInd,:).*proj;
tmpX = sum(tmpX,2);
whos tmpX a proj
plot(abs(tmpX))



aX = uX*sX*vX';
whos uX sX vX aX

dataX = aX/proj;
whos data dataX


for mk=1:mdkp,
    sp(:,j,mk)=u(:,mk)';
    fm(:,j,mk)=v(:,mk)';
end
sv(j,:)=diag(s);
end


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

% %% Reconstruct psd at each mode (simply from v)
% nTaper = param.tapers(2);
% tmp = permute(v,[2 1]);
% sz = size(tmp);
% sz(3) = sz(2)/nTaper;
% sz(2) = param.tapers(2);
% tmp = reshape(tmp,sz);
% psdRecSimple = permute(mean(conj(tmp).*tmp,2),[3 1 2]);
% 
% %% Plot reconstructions
% [~,pad,~,~,~,~,~]=getparams(param);
% N=size(funTs.vec,1);
% nfft=max(2^(nextpow2(N)+pad),N);
% [f,~]=getfgrid(param.Fs,nfft,param.fpass);
% 
% close all
% figure('WindowStyle','docked');
% plot(diag(s))
% ax = gca;
% ax.YScale = 'log';
% for modInd = 1:30
%     figure('WindowStyle','docked');
%     tl = tiledlayout(3,2);
%     tl.TileSpacing = "tight"; tl.Padding = "tight";
%     tl.TileIndexing = 'rowmajor';
%     title(tl,['mode ' num2str(modInd)])
%     
%     tmp = nan(size(mask));
%     tmp(logical(mask)) = u(:,modInd);
%     
%     nexttile(1,[2 1])
%     hIm = imagesc(abs(tmp));
%     hIm.AlphaData = ~isnan(tmp);
%     ax = gca;
%     ax.YTick = []; ax.XTick = [];
%     ax.ColorScale = 'linear'; ax.Colormap = jet;
%     ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
%     ylabel(colorbar,'psd')
% 
%     nexttile(2,[2 1])
%     hIm = imagesc(angle(tmp));
%     hIm.AlphaData = ~isnan(tmp);
%     ax = gca;
%     ax.YTick = []; ax.XTick = [];
%     ax.ColorScale = 'linear'; ax.Colormap = hsv;
%     ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
%     ax.CLim = [-pi pi];
%     ylabel(colorbar,'phase')
% 
%     nexttile(5,[1 2])
% %     plot(f,mean(psdRec(:,:,modInd),1))
%     plot(f,psdRecSimple(:,modInd))
%     ylabel('psd')
%     xlabel('Hz')
%     drawnow
% end


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
if exist('vecNorm','var')
    svdStruct.normFact = vecNorm;
else
    svdStruct.normFact = [];
end
svdStruct.tsMean = tsMean;
svdStruct.dim = strjoin({'space/taper' 'freq/time' 'modes'},' X ');
svdStruct.sv = permute(sv,[3 1 2]);
svdStruct.sp = sp; %spatial singular vectors
svdStruct.fm = fm; %taper singular vectors
svdStruct.tp = permute(tapers,[2 1]); %tapers
svdStruct.tv = permute(tvec(:,1),[2 1]); %time vector
svdStruct.f = f; %freq vector
svdStruct.tf = 'tp.*exp(-f(i)*tv)'; %taper frequency time vector
svdStruct.c = sv(:,1)'.^2./sum(sv.^2,2)';
svdStruct.w = Wreal;
svdStruct.param = param;


if 0 
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

nF = min(size(svdStruct.sp,3),30)
for ind = 1:nF
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
    ylabel(colorbar,'rad')

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
end

