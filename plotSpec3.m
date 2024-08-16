function [ax,F] = plotSpec3(H,volPsd,metricLabel,dsgn,mask,volResp,respQthresh,tWin)
if ~exist('H','var'); H = [];                             end
if isempty(H);        H = figure('WindowStyle','docked'); end

% if isfield(volPsd,'dsgn') && isfield(volPsd.dsgn,'f0'); f0 = volPsd.dsgn.f0; else; f0 = []; end

if ~exist('tWin','var');                                     tWin = []; end
if ~exist('metricLabel','var');                       metricLabel = []; end
if ~exist('volTs','var');                                   volTs = []; end
if ~exist('respQthresh','var');                       respQthresh = []; end
if ~exist('dsgn','var');                                     dsgn = []; end
if ~exist('mask','var');                                     mask = []; end
if isempty(metricLabel);                              metricLabel = 'psd'; end % 'psd' 'coh'
if isempty(respQthresh) && strcmp(metricLabel,'psd'); respQthresh = 1; end % 'psd' 'coh'
% threshAvFlag = ismember(metricLabel,{'psd' 'psdEPC'}) && ~isempty(volTs) && isfield(volTs,'resp') && ~isempty(respQthresh) && respQthresh~=inf && respQthresh~=1;
if isempty(dsgn)
    if isfield(volTs,'dsgn');       dsgn = volTs.dsgn;
    elseif isfield(volResp,'dsgn'); dsgn = volResp.dsgn;
    elseif isfield(volPsd,'dsgn');  dsgn = volPsd.dsgn;
    end
end

switch class(H)
    case 'matlab.graphics.layout.TiledChartLayout'
        F = H.Parent;
    case 'matlab.ui.Figure'
        F = H;
    otherwise
end

figure(F);
ax = {};
ax{end+1} = nexttile;


%% Select approrpiate data
switch metricLabel
    case 'psd'
        mt    = volPsd.psd;
        spc   = mt.PSD;
        label = 'spatially averaged spectrum';
    case 'coh'
        mt    = volPsd.svd;
        spc   = mt.COH;
        label = 'coherence spectrum';
    otherwise
        dbstack; error('code trhat')
end




switch metricLabel
    case 'psd'
        %% average across space
        %%% simple crop
        if isempty(mask)
            mask = getCropMask(volPsd);
        else
            if isMRI(mask)
                mask = logical(mask.vol);
            else
                mask = MRIread(mask);
                mask = logical(mask.vol);
            end
            mask = mask & getCropMask(volPsd);
        end
        %%% stat thresh
        if respQthresh==0; respQthresh = 0.05; end % respQthresh = 0 does not mean anything, reverts to default 0.05
        if respQthresh~=inf && respQthresh~=1
            mask = mask & volResp.Fq.vol<respQthresh;
        end
        %%% apply
        if any(size(spc,[1:4 7 8])~=1); dbstack; error('something unexpected here'); end
        spc = permute(spc(:,:,:,:,:,mask(volPsd.vol2vec),:,:),[5 6 1 2 3 4 7 8]);
        spc = mean(spc,2);
        f = permute(mt.f,[5 1 2 3 4 6 7 8]);
    case 'coh'
        %% coherence is already summarizing space
        m = 1;
        if any(size(spc,[1:4 6:7])~=1); dbstack; error('something unexpected here'); end
        spc = permute(spc(:,:,:,:,:,:,:,m),[5 1 2 3 4 6 7 8]);
        f = permute(mt.f,[5 1 2 3 4 6 7 8]);
    otherwise
        dbstack; error('code trhat')
end





%% Plot
plot(f,spc,'k')
grid on
axis tight
xlabel('f (Hz)')
switch metricLabel
    case 'psd'
        ylabel('psd')
        ax{end}.YScale = 'log';
    case 'coh'
        ylabel('coherence')
    otherwise
        dbstack; error('code trhat')
end


K = mt.K;
T = mt.T;
[TW,W] = K2W(T,K,0);

paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec; TW=' num2str(TW) ')'];
switch metricLabel
    case 'psd'
        if respQthresh~=inf && respQthresh~=1
            title(['spectrum ' paramStr ' cross-vox (Q<=' num2str(respQthresh,'%0.2f') ') mean'])
        else
            title(['spectrum ' paramStr ' cross-vox mean'])
        end
    case 'coh'
        title(['coherence spectrum ' paramStr])
end

yLim = spc(f>0.01);
yLim = [min(yLim) max(yLim)];
ylim(yLim)
xlim([0 f(end)])

addW([],volPsd.psd)



% if isfield(volPsd,'psdTrialGram') && ~isempty(volPsd.psdTrialGram)
%     addFreq([],volPsd.psdTrialGram.onsetList,volPsd.psdTrialGram.durList,2)
% end
if ~isempty(dsgn)
    addFreq([],dsgn.onsetList,dsgn.ondurList)
end

% if ~isempty(f0)
%     xline(f0,'--b')
% end

if isfield(volPsd,'psdTrialGramMD') && ~isempty(volPsd.psdTrialGramMD) && ~isempty(tWin)
    dbstack; error('double-check that');
    t   = volPsd.psdTrialGramMD.t(1,1,1,1,1,1,:,1);
    f   = volPsd.psdTrialGramMD.f(1,1,1,1,:,1,1,1);
    if tWin==inf
        psd = mean(mean(volPsd.psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,:,:),6),7);
    else
        [~,wInd] = min(abs(t-tWin));
        psd = mean(volPsd.psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,wInd,:),6);
    end
    plot(squeeze(f),squeeze(psd),'--k')
end


ax = [ax{:}];




