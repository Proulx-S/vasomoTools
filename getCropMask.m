function mask = getCropMask(mri,crop)
if ~exist('crop','var'); crop = []; end
if isempty(crop);        crop = 5; end

if ismember(fields(mri),'vec'); mri = rmfield(mri,'vec'); end

trim = ceil(crop./mri.volres); %voxel
trimActual = trim.*mri.volres; %mm
mask = true(mri.volsize);
mask([1:trim(1) end-trim(1)+1:end],:                            ,:) = false;
mask(:                            ,[1:trim(2) end-trim(2)+1:end],:) = false;

trimZ = trim(3);
while mri.volsize(3) - trimZ*2 < 5
    trimZ = trimZ-1;
end
mask(:,:,[1:trimZ end-trimZ+1:end]) = false;




