function mri = applyMask(mri,mask)

if ischar(mask)
    mask = MRIread(mask);
    mask = mask.vol;
elseif isMRI(mask)
    mask = mask.vol;
end

mask = logical(mask);


for sInd = 1:length(mri)
    if ~isMRI(mri(sInd)) && isfield(mri(sInd),'mri')
        mri(sInd).mri = applyMask(mri(sInd).mri,mask);
        continue
    end

    if ~isempty(mri(sInd).vol)   ||   ( isfield(mri(sInd),'vec') && ~isempty(mri(sInd).vec) )
        if ~isempty(mri(sInd).vol)
            mri(sInd).vol = permute(mri(sInd).vol,[4 5 6 1 2 3]);
            mri(sInd).vol(:,:,:,~mask) = nan;
            mri(sInd).vol = permute(mri(sInd).vol,[4 5 6 1 2 3]);
            mri(sInd).vec = [];
        else
            mri(sInd).vec = mri(sInd).vec(:,mask(mri(sInd).vol2vec),:,:);
        end
    end

    mri(sInd).vol2vec = mask;
end