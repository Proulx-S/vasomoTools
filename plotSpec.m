function [ax,F] = plotSpec(volPsd,metricLabel,H,tWin,volTs,respQthresh)
if ~exist('H','var'); H = [];                             end
if isempty(H);        H = figure('WindowStyle','docked'); end

if isfield(volPsd,'dsgn') && isfield(volPsd.dsgn,'f0'); f0 = volPsd.dsgn.f0; else; f0 = []; end

if ~exist('tWin','var');                                     tWin = []; end
if ~exist('metricLabel','var');                       metricLabel = []; end
if ~exist('volTs','var');                                   volTs = []; end
if ~exist('respQthresh','var');                       respQthresh = []; end
if isempty(metricLabel);                              metricLabel = 'psd'; end % 'psd' 'coh'
if isempty(respQthresh) && strcmp(metricLabel,'psd'); respQthresh = 0.05; end % 'psd' 'coh'
threshAvFlag = ismember(metricLabel,{'psd' 'psdEPC'}) && ~isempty(volTs) && isfield(volTs,'resp') && ~isempty(respQthresh) && respQthresh~=inf && respQthresh~=1;

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
        vec   = mt.PSD;
        label = 'spatially averaged spectrum';
    case 'coh'
        mt    = volPsd.svd;
        vec   = mt.COH;
        label = 'coherence spectrum';
    otherwise
        dbstack; error('code trhat')
end




%%% average across space
if threshAvFlag
    vec = squeeze(mean(vec(:,:,:,:,:,volTs.resp.Fq.vol(volPsd.vol2vec)<=respQthresh,:,1),6));
else
    vec = squeeze(mean(vec(:,:,:,:,:,:,:,1),6));
end
f = squeeze(mt.f);
% f   = squeeze(volPsd.f);
% psd = squeeze(mean(volPsd.vec,2));
plot(f,vec,'k')
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
if threshAvFlag
    title(['spectrum ' paramStr ' averaged across voxels with Q<=' num2str(respQthresh,'%0.2f')])
else
    switch metricLabel
        case 'psd'
            title(['spectrum ' paramStr ' averaged across voxels'])
        case 'coh'
            title(['coherence spectrum ' paramStr])
    end
end

yLim = vec(f>0.01);
yLim = [min(yLim) max(yLim)];
ylim(yLim)
xlim([0 f(end)])

addW([],volPsd.psd)



if isfield(volPsd,'psdTrialGram') && ~isempty(volPsd.psdTrialGram)
    addFreq([],volPsd.psdTrialGram.onsetList,volPsd.psdTrialGram.durList,2)
end

if ~isempty(f0)
    xline(f0,'--b')
end

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




