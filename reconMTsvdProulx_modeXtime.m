function rec = reconMTsvdProulx_modeXtime(svdStruct,funTs,modeOrder)
if ~exist('modeOrder','var') || isempty(modeOrder)
    modeOrder = length(svdStruct.sv);
end

sv = svdStruct.sv(:,:,1:modeOrder); % spatial sv (space x 1 x mode)
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
