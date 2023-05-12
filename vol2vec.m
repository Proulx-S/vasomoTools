function mri = vol2vec(mri,mask)
if exist('mask','var') && ~isempty(mask) && ~isempty(mri.vol)
    if isfield(mri,'vol2vec')
        error('mask already exists')
    end
    mri.vol2vec = mask;
    mri.vec = [];
    curNames = fieldnames(mri);
    curInd = find(ismember(curNames,{'vol' 'vol2vec' 'vec'}));
    tmpNames1 = curNames(1:(curInd(1)-1));
    tmpNames2 = curNames((curInd(1)):end); tmpNames2(ismember(tmpNames2,{'vol' 'vol2vec' 'vec'})) = [];
    newNames = [tmpNames1; {'vol' 'vol2vec' 'vec'}'; tmpNames2];
    [~,Locb] = ismember(newNames,curNames);
    mri = orderfields(mri,Locb);
end

mri.vol = permute(mri.vol,[4 1 2 3]);
mri.vec = mri.vol(:,mri.vol2vec);
mri.vol = [];

