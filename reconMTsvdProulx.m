function svdTs = reconMTsvdProulx(svdStruct,funTs,modeOrder)
if ~exist('modeOrder','var') || isempty(modeOrder)
    modeOrder = length(svdStruct.sv);
end

sv = svdStruct.sv(:,:,1:modeOrder); % spatial sv (space x 1 x mode)
sp = svdStruct.sp(:,:,1:modeOrder); % spatial sv (space x 1 x mode)
fm = svdStruct.fm(:,:,1:modeOrder); % taper sv (taper x 1 x mode)
tp = svdStruct.MTS.proj; % tapers (taper x time)
if isfield(svdStruct.MTS,'NRUN') && ~isempty(svdStruct.MTS.NRUN)
    tp = repmat(tp,[svdStruct.MTS.NRUN 1]);
end
f = svdStruct.MTS.bandFreq; % frequencies (1 x freq)
fInd = svdStruct.MTS.bandInd; % frequency indices (taper x 1)
tsRec = squeeze(fm)'*tp; % WARNING! taking the conjugate here with the ' operator
tsRec = tsRec .* squeeze(sv); % scale according to singular value
svdTs = funTs(1);
svdTs.vol = permute(tsRec,[1 3 4 2]);
svdTs.volsize = size(svdTs.vol,1:3); svdTs.height = svdTs.volsize(1); svdTs.width = svdTs.volsize(2); svdTs.depth = svdTs.volsize(3); svdTs.nvoxels = prod(svdTs.volsize);
svdTs.param.NRUN = svdStruct.MTS.NRUN;
svdTs.param.K = svdStruct.MTS.K;
