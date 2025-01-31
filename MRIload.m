function mri = MRIload(f,volInfo,r)
if issctruct(f)
    if isMRI(f)
        f
    else
        f
    end
    return
end
if iscell(f)
    for r = 1:length(f)
        mri(r,1) = MRIload(f{r},volInfo,r);
    end
    return
end

%% Read header
mri = MRIread(f,1);

%% Copy relevant fields from volInfo
if exist('volInfo','var') && ~isempty(volInfo)
    copyFieldList = {'sub' 'ses' 'label' 'labelAcq' 'dsgn' 'wd' 'bidsDir' 'bidsDerivDir' 'ppLabelList' 'dataType' 'fOrigList' 'fPreprocList' 'fTransList' 'fTransCatList' 'bidsList' 'acqTime' 'nDummy'  'volAnat' 'volAnatSes' 'volAnatSub' 'nFrame'};

    for i = 1:length(copyFieldList)
        if isfield(volInfo,copyFieldList{i})
            if size(volInfo.(copyFieldList{i}),1)==1
                mri.(copyFieldList{i}) = volInfo.(copyFieldList{i});
            else
                mri.(copyFieldList{i}) = volInfo.(copyFieldList{i})(r);
            end
        end
    end
    nFrame = volInfo.nFrame(r);
else
    nFrame = mri.nframes;
end
nFrameOrig = MRIread(volInfo.fOrigList{r},1); nFrameOrig = nFrameOrig.nframes;
mri.nDummyRemoved = nFrameOrig-nFrame;
