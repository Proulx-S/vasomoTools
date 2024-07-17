function svdPsd = reconMTsvdProulx_modePsd(svdStruct,funPsd,modeOrder)
if ~exist('modeOrder','var') || isempty(modeOrder)
    modeOrder = length(svdStruct.sv);
end

sv = svdStruct.sv(:,:,1:modeOrder); % spatial sv (space x 1 x mode)
sp = svdStruct.sp(:,:,1:modeOrder); % spatial sv (space x 1 x mode)
fm = svdStruct.fm(:,:,1:modeOrder); % taper sv (taper x 1 x mode)
tp = svdStruct.MTS.proj; % tapers (taper x time)
f = svdStruct.MTS.bandFreq; % frequencies (1 x freq)
fInd = svdStruct.MTS.bandInd; % frequency indices (taper x 1)
svdPsd = rmfield(vol2vec(catRuns(funPsd)),'vec');
sp = squeeze(sp);
sp = sp./sum(abs(sp),1);
if ~isfield(svdStruct.MTS,'NRUN') || isempty(svdStruct.MTS.NRUN); svdStruct.MTS.NRUN = 1; end
for run = 1:svdStruct.MTS.NRUN
    for k = 1:funPsd(run).ntapers
        vec = squeeze(funPsd(run).vec(:,:,k));
        svdPsd.vec(:,:,k,run) = sqrt(vec*sp);
    end
end
svdPsd.volsize = [modeOrder 1 1]; svdPsd.height = svdPsd.volsize(1); svdPsd.width = svdPsd.volsize(2); svdPsd.depth = svdPsd.volsize(3); svdPsd.nvoxels = prod(svdPsd.volsize);
svdPsd.ntapers = funPsd(1).ntapers;
svdPsd.svd.param = svdStruct.param;
svdPsd.svd.MTS = svdStruct.MTS; svdPsd.svd.MTS.proj = [];