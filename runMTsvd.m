function svdStruct = runMTsvd(anaType,funTs,fpass,W,mask)
% Similar to Mitra 1997. A single svd is run on data tapered for
% sensitivity over user-defined frequency band (fpass).
funTs = vol2vec(funTs,mask);
tr = funTs.tr/1000;

%% Detrend time series (detrend up to order-2 polynomial, since this is the highest order not fitting a sinwave)
funTs = dtrnd4psd(funTs);

%% Zscore time series
funTs.vec = zscore(funTs.vec,[],1);

%% Set multitaper parameters
param.tapers = [];
param.Fs = 1/tr;
switch anaType
    case 'svdProulx'
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
        mdkp = 25;
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
        


    case 'svdKlein' | 'svdProulx'
        % if exist('W','var') && ~isempty(W)
        %     anaType = 'svdKlein';
        T = tr.*funTs.nframes;
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
        T = tr.*funTs.nframes;
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
[sv,sp,fm] = spsvd(funTs.vec,param,mdkp);
toc

%% Output
svdStruct.mask = mask;
svdStruct.dim = strjoin({'space/taper' 'freq' 'modes'},' X ');
svdStruct.sv = permute(sv,[3 1 2]);
svdStruct.sp = sp;
svdStruct.fm = fm;
svdStruct.c = sv(:,1)'.^2./sum(sv.^2,2)';
svdStruct.f = f;
svdStruct.w = Wreal;
svdStruct.param = param;
