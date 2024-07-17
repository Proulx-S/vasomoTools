function mri = applyRoi(mri,volRoi,roiInd,combineFlag)
if ~exist('combineFlag','var'); combineFlag = false; end
if ~exist('roiInd','var');           roiInd = []; end
if isempty(roiInd);                  roiInd = 1:length(volRoi.label); end

if combineFlag
    mri = vol2vec(vec2vol(mri),any(volRoi.mri.vol,4),1);
    mri.vol2vecFlag = strjoin(volRoi.label,'+');
else
    mri2 = vec2vol(mri); clear mri
    for i = 1:length(roiInd)
        mri(i) = vol2vec(mri2,volRoi.mri.vol(:,:,:,roiInd(i)),1);
    end
    [mri(:).vol2vecFlag] = deal(volRoi.label{roiInd});    
end