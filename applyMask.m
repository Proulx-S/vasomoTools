function mri = applyMask(mri,mask)

mask = logical(mask);

for sInd = 1:length(mri)
    if ~isempty(mri(sInd).vol)
        mri(sInd).vol = permute(mri(sInd).vol,[4 5 6 1 2 3]);
        mri(sInd).vol(:,:,:,~mask) = nan;
        mri(sInd).vol = permute(mri(sInd).vol,[4 5 6 1 2 3]);
    else
        mri(sInd).vec = mri(sInd).vec(:,mask(mri(sInd).vol2vec),:,:);
    end
    mri(sInd).vol2vec = mask;
end