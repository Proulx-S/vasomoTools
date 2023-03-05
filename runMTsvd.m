function funTs = runMTsvd(funTs,fpass,W)
% Similar to Mitra 1997. A single svd is run on data tapered for
% sensitivity over user-defined frequency band (fpass).

if isempty(funTs.vec)
    funTs = vol2vec(funTs);
end
tr = funTs.tr/1000;

%% Detrend time series (detrend up to order-2 polynomial, since this is the highest order not fitting a sinwave)
funTs = dtrnd4psd(funTs);

%% Zscore time series
funTs.vec = zscore(funTs.vec,[],1);

%% Set multitaper parameters
param.tapers = [];
param.Fs = 1/tr;
if exist('W','var') && ~isempty(W)
    flag = 'svdKlein';
    T = tr.*funTs.nframes;
    TW = T*W;
    K = round(TW*2-1);
    TW = (K+1)/2;
    param.tapers = [TW K];
    if ~isempty(fpass)
        param.fpass = fpass;
    end
    %%% Display actual half-widht used
    Wreal = TW/T;
    display(['w  (halfwidth) requested  : ' num2str(W,'%0.5f ')])
    display(['w  (halfwidth) used       : ' num2str(Wreal,'%0.5f ')])
    display(['tw (time-halfwidth) used  : ' num2str(TW)])
    display(['k  (number of tapers) used: ' num2str(K)])
else
    flag = 'svdMitra';
    %%% Dummy mtspec just to get f
    param.tapers = 1.5; param.tapers(2) = param.tapers(1)*2-1;
    [~,f] = mtspectrumc(funTs.vec(:,1), param);
    %%% Set parameters for the user-defined frequency band
    W = diff(fpass)/2;
    f0 = fpass(1)+W; [~,b] = min(abs(f - f0)); f0 = f(b);
    T = tr.*funTs.nframes;
    TW = T*W;
    K = round(TW*2-1);
    TW = (K+1)/2;
    param.tapers = [TW K];
    param.fpass = [f0 f0];
    %%% Display actual frequency band used
    fpassReal = f0+[-1 1].*(TW/T);
    display(['frequency band requested: fpass=[' num2str(fpass,'%0.5f ') ']'])
    display(['frequency band used     : fpass=[' num2str(fpassReal,'%0.5f ') ']'])
end

%% Run the decomposition
[sv,sp,fm] = spsvd(funTs.vec,param);

%% Output
funTs.(flag).sv = sv;
funTs.(flag).sp = sp;
funTs.(flag).fm = fm;
funTs.(flag).param = param;

% 
% %% Visualize
% figure('WindowStyle','docked');
% h1 = plot(sv,'k'); hold on
% grid on
% ylabel('singular value')
% xlabel('k');
% title(['subj' subjId_list{subjInd} 'run' num2str(rInd) '; ' preprocStr(2:end) '; tw=' num2str(tw)],'Interpreter','none')
% K=19; % 283543-> 24
% 
% %         nPerm = 1000;
% %         sz = size(sv); sz(3) = nPerm;
% %         svPerm = nan(sz);
% %         nframes = funTs.nframes;
% %         tic
% %         parfor indPerm = 1:nPerm
% %             rPerm = randperm(nframes);
% %             [svPerm(:,:,indPerm),~,~] = spsvd(funTs.vec(rPerm,:),param);
% %         end
% %         svPerm(:,:,any(isnan(svPerm(1,:,:)),2)) = [];
% %         nPerm = size(svPerm,3);
% %         tH = toc;
% %         svPermPrct95 = prctile(svPerm,95,3);
% %         svPermAv = mean(svPerm,3);
% %
% %         h2 = plot(svPermAv,':r');
% %         h3 = plot(svPermPrct95,'r');
% %         yyaxis right
% %         percAbovePerm = (sv-svPermAv)./svPermAv.*100;
% %         h4 = plot(percAbovePerm);
% %         h5 = plot(xlim,[1 1].*thresh,'--','Color',h4.Color);
% %         ylabel('% above mean of permuted sv')
% %         thresh = 10;
% %         K = find(percAbovePerm<thresh,1)-1;
% h6 = plot([1 1].*K+0.5,ylim,'--k');
% %         legend([h1 h2 h3 h5 h6],{'real' [num2str(nPerm) 'permAv'] [num2str(nPerm) 'perm95'] 'thresh' 'cutoff'},'box','off')
% 
% sv(:,K+1:end) = [];
% sp(:,:,K+1:end) = [];
% fm(:,:,K+1:end) = [];
% 
% % Extract singular-vector-weighted component psd
% cmpntPsd = nan(length(f),K);
% for k = 1:K
%     spCur = abs(sp(:,:,k));
%     spCur = spCur./sum(spCur);
%     cmpntPsd(:,k) = funPsd.vec*spCur;
% end
% 
% figure('WindowStyle','docked');
% plot(f,mean(funPsd.vec,2),'k'); hold on
% plot(f,cmpntPsd(:,1))
% plot(f,cmpntPsd(:,2))
% plot(f,cmpntPsd(:,3))
% plot(f,cmpntPsd(:,4))
% plot(f,cmpntPsd(:,end))
% ax = gca; ax.YScale = 'log'; ax.XScale = 'log';
% grid on
% xlabel('Hz'); ylabel('PSD');
% legend({'unweighted spatial averaged' 'cmpnt 1 weighted' 'cmpnt 2 weighted' 'cmpnt 3 weighted' 'cmpnt 4 weighted' 'last cmpnt weighted'},'box','off')
% xlim(f([1 end]))
% title(['subj: ' num2str(subjId_list{subjInd}) '; ' preprocStr(2:end)],'Interpreter','none')
% 
% % Decompose these component psd to remove the mean spectrum
% fInd = fpassReal(1)<=f & fpassReal(2)>=f;
% [U,S,V] = svd(log(cmpntPsd(fInd,:)),0);
% figure('WindowStyle','docked');
% plot(diag(S))
% k = 1;
% plot(f,-0.23*S(k,k).*U(:,k))
% hold on
% plot(f,log(mean(cmpntPsd,2)))
% 
% psdRecon = exp(U(:,2:end)*S(2:end,2:end)*V(:,2:end)');
% figure('WindowStyle','docked');
% plot(f,psdRecon(:,1)); hold on
% plot(f,psdRecon(:,2))
% plot(f,psdRecon(:,3))
% plot(f,psdRecon(:,4))
% plot(f,psdRecon(:,end))
% ax = gca; ax.YScale = 'log'; ax.XScale = 'log';
% grid on
% xlabel('Hz'); ylabel('PSD');
% legend({'cmpnt 1 weighted' 'cmpnt 2 weighted' 'cmpnt 3 weighted' 'cmpnt 4 weighted' 'last cmpnt weighted'},'box','off')
% xlim(f([1 end]))
% title(['subj: ' num2str(subjId_list{subjInd}) '; ' preprocStr(2:end)],'Interpreter','none')


