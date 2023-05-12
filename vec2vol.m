function mri = vec2vol(mri)

mri.vol = nan(mri.nframes,mri.height,mri.width,mri.depth);
mri.vol(:,mri.vol2vec) = mri.vec;
mri.vec = [];
mri.vol = permute(mri.vol,[2 3 4 1]);