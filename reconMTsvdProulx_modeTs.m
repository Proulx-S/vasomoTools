function svdTs = reconMTsvdProulx_modeTs(svdStruct,funTs,modeOrder)
if ~exist('modeOrder','var') || isempty(modeOrder)
    modeOrder = length(svdStruct.sv);
end

sv = svdStruct.sv(:,:,1:modeOrder); % spatial sv (space x 1 x mode)
sp = svdStruct.sp(:,:,1:modeOrder); % spatial sv (space x 1 x mode)
fm = svdStruct.fm(:,:,1:modeOrder); % taper sv (taper x 1 x mode)
tp = svdStruct.MTS.proj; % tapers (taper x time)
f = svdStruct.MTS.bandFreq; % frequencies (1 x freq)
fInd = svdStruct.MTS.bandInd; % frequency indices (taper x 1)
svdTs = rmfield(funTs,'vol'); clear funTs
if ~isfield(svdStruct.MTS,'NRUN') || isempty(svdStruct.MTS.NRUN); svdStruct.MTS.NRUN = 1; end
for run = 1:svdStruct.MTS.NRUN
    runInd = (run-1)*svdStruct.MTS.K+1:(run)*svdStruct.MTS.K;
    tsRec = permute(squeeze(fm(runInd,:,:)),[2 1])*tp;
    % tsRec = squeeze(fm(runInd,:,:))'*tp; % WARNING! taking the conjugate here with the ' operator
    tsRec = tsRec .* squeeze(sv); % scale according to singular value
    % 
    svdTs(run).vol = permute(tsRec,[1 3 4 2]);
    svdTs(run).volsize = size(svdTs(run).vol,1:3); svdTs(run).height = svdTs(run).volsize(1); svdTs(run).width = svdTs(run).volsize(2); svdTs(run).depth = svdTs(run).volsize(3); svdTs(run).nvoxels = prod(svdTs(run).volsize);
    svdTs(run).nruns = 1;
    svdTs(run).svd.param = svdStruct.param;
    svdTs(run).svd.MTS = svdStruct.MTS; svdTs(run).svd.MTS.proj = [];
end
% % tsRec = squeeze(fm)'*tp; % WARNING! taking the conjugate here with the ' operator
% % tsRec = tsRec .* squeeze(sv); % scale according to singular value
% svdTs = funTs(1);
% svdTs.vol = permute(cat(6,tsRec{:}),[1 3 4 2 5 6]);
% % svdTs.vol = permute(tsRec,[1 3 4 2]);
% % svdTs.volsize = size(svdTs.vol,1:3); svdTs.height = svdTs.volsize(1); svdTs.width = svdTs.volsize(2); svdTs.depth = svdTs.volsize(3); svdTs.nvoxels = prod(svdTs.volsize);
% svdTs.volsize = size(svdTs.vol,1:3); svdTs.height = svdTs.volsize(1); svdTs.width = svdTs.volsize(2); svdTs.depth = svdTs.volsize(3); svdTs.nvoxels = prod(svdTs.volsize);
% svdTs.nruns = svdStruct.MTS.NRUN;
% svdTs.svd.param = svdStruct.param;
% svdTs.svd.MTS = svdStruct.MTS; svdTs.svd.MTS.proj = [];
