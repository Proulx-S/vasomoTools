function mri = setNiceFieldOrder(mri,fieldList)
% fieldList = {'vol' 'vol2vec' 'vec'};
%% Add fields in fieldList that don't exist
for fieldInd = 1:length(fieldList)
    if ~isfield(mri,fieldList{fieldInd}); mri.(fieldList{fieldInd}) = []; end
end
%% Order fieldList fields
curNames = fieldnames(mri);
[tmp1,~] = ismember(curNames,fieldList);
mri = orderfields(mri,[curNames(~tmp1); fieldList']);
% [tmp1,tmp2] = ismember(curNames,fieldList);
% tmp3 = curNames(tmp1);
% tmp3 = tmp3(tmp2(tmp2~=0));
% mri = orderfields(mri,[curNames(~tmp1); tmp3]);