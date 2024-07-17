function volTs = volTsRecon(modeTs,modeOrder)
% modeTs = volTs; modeTs.vol = []; if isfield(modeTs,'vec'); modeTs.vec = []; end
if ~exist('modeOrder','var') || isempty(modeOrder); modeOrder = 1:length(modeTs.modeInd); end
volTs = modeTs;
volTs.vec = []; volTs.vol = [];

% u = permute(modeTs.svd.u(:,modeOrder,:,:,:,:),[3 1 6 4 2 5]);
% u = permute(u,[3 1 6 4 2 5]);
% modeInfo = strsplit(modeTs.vecInfo,' x ');
% uInfo = strsplit(modeTs.svd.info,' x ');
% uInfo = uInfo([3 1 6 4 2 5]);
% modeInfo
% uInfo

% volTs.vec = modeTs.vec(:,:,:,:,modeOrder) .* permute(modeTs.svd.u(:,modeOrder,:,:,:,:),[3 1 6 4 2 5]);
volTs.vec = conj(modeTs.vec(:,:,:,:,modeOrder)) .* permute(modeTs.svd.u(:,modeOrder,:,:,:,:),[3 1 6 4 2 5]);

%% Housekeeping
[volTs.height,volTs.width,volTs.depth] = size(modeTs.svd.mask,[1 2 3]);
volTs.nvoxels = prod(size(modeTs.svd.mask,[1 2 3]));
volTs.ntapers = size(volTs.vec,3);
volTs.vol2vec = modeTs.svd.mask;
volTs.modeInd = modeOrder;


return

%% Get indexes for the time window
[runInd,winCentTimeInd] = find(ismember(squeeze(volPsd.svdGram.tWin),winCentTime));
[~,winFreqInd] = min(abs(volPsd.svdGram.f-winFreq));

%% Get frequency modulated tapers
t = permute(time2winTime(1/volPsd.svdGram.param.Fs,winCentTime,volPsd.svdGram.lWin),[1 3 2]);
f = volPsd.svdGram.f(winFreqInd);
tp = volPsd.svdGram.tp;
% tp = permute(volPsd.svdGram.tp,[6 3 1 2 4 5]);
ftp = tp .* exp(-f*(t*2*pi*1i));

%% Get svd coefficients
s = volPsd.svdGram.s(:,modeOrder,winFreqInd,runInd,winCentTimeInd,:); % singular values
v = volPsd.svdGram.v(:,modeOrder,winFreqInd,runInd,winCentTimeInd,:); % freqTime singular vectors
% s = permute(volPsd.svdGram.s(:,modeOrder,winFreqInd,runInd,winCentTimeInd,:),[6 2 3 1 4 5]); % freqTime singular vector
% v = permute(volPsd.svdGram.v(:,modeOrder,winFreqInd,runInd,winCentTimeInd,:),[6 2 3 1 4 5]); % freqTime singular vector

%% Recon mode Ts
ts = v.*ftp.*s;


%% Housekeeping
tmp = strsplit(volPsd.svdGram.info,' x '); tmp = tmp([3 1 6 4 2 5]);
modeTs.vecInfo = strsplit(modeTs.vecInfo,' x '); modeTs.vecInfo = [modeTs.vecInfo tmp(length(modeTs.vecInfo)+1:end)]; modeTs.vecInfo = strjoin(modeTs.vecInfo,' x ');
modeTs.vec = permute(ts,[3 1 6 4 2 5]); clear ts
modeTs.t = permute(t,[3 1 6 4 2 5]);
modeTs.height = 1; modeTs.width = 1; modeTs.depth = 1; modeTs.nvoxels = 1; modeTs.nframes = length(t); modeTs.nruns = 1;
modeTs.modeInd = modeOrder;
modeTs.K = volPsd.svdGram.K;
% modeTs.svdGram = rmfield(volPsd.svdGram,'tp','u','s','v','coh','T','tWin','lWin',)
if isfield(modeTs,'vol2vec'); modeTs = rmfield(modeTs,'vol2vec'); end



return
sv = svdStruct.sv(:,:,1:modeOrder);
sp = svdStruct.sp(:,:,1:modeOrder); % spatial sv (space x 1 x mode)
fm = svdStruct.fm(:,:,1:modeOrder); % taper sv (taper x 1 x mode)
tp = svdStruct.MTS.proj; % tapers (taper x time)
f = svdStruct.MTS.bandFreq; % frequencies (1 x freq)
fInd = svdStruct.MTS.bandInd; % frequency indices (taper x 1)

rec = funTs;
rec.vol = squeeze(fm)'*tp; % WARNING! taking the conjugate here with the ' operator
rec.vol = rec.vol .* squeeze(sv); % scale according to singular value
rec.vol = permute(rec.vol,[1 3 4 2]); rec.volMean = mean(rec.vol,4);
rec.volsize = size(rec.vol,1:3); rec.height = rec.volsize(1); rec.width = rec.volsize(2); rec.depth = rec.volsize(3); rec.nvoxels = prod(rec.volsize);
