function mri = ses2ser(mri)
if ~isfield(mri,'f')
    if isfield(mri,'vol') && ~isempty(mri.vol)
        perm1 = [1 2 3 5 4 6]; perm2(perm1) = 1:length(perm1);
        fieldList = {'vol'};
        for fieldCur = fieldList
            fieldCur = char(fieldCur);
            sz = size(mri.(fieldCur),min(perm1):max(perm1)); sz = sz(perm1); sz(end-1) = prod(sz(end-1:end)); sz(end) = [];
            mri.(fieldCur) = permute(reshape(permute(mri.(fieldCur),perm1),sz),perm2);
        end
        % mri.vol = permute(mri.vol,[1 2 3 5 4 6]);
        % mri.vol = permute(mri.vol(:,:,:,:,:),[1 2 3 5 4]);
    end
    if isfield(mri,'vec') && ~isempty(mri.vec)
        perm1 = [2 3 1 4]; perm2(perm1) = 1:length(perm1);
        fieldList = {'vec'};
        for fieldCur = fieldList
            fieldCur = char(fieldCur);
            sz = size(mri.(fieldCur),min(perm1):max(perm1)); sz = sz(perm1); sz(end-1) = prod(sz(end-1:end)); sz(end) = [];
            mri.(fieldCur) = permute(reshape(permute(mri.(fieldCur),perm1),sz),perm2);
        end
        % mri.vec = permute(mri.vec,[2 3 1 4]);
        % mri.vec = permute(mri.vec(:,:,:),[3 1 2]);
    end
end

perm1 = [1 2 3 5 4 6]; perm2(perm1) = 1:length(perm1);
fieldList = {'t'};
for fieldCur = fieldList
    fieldCur = char(fieldCur);
    sz = size(mri.(fieldCur),min(perm1):max(perm1)); sz = sz(perm1); sz(end-1) = prod(sz(end-1:end)); sz(end) = [];
    mri.(fieldCur) = permute(reshape(permute(mri.(fieldCur),perm1),sz),perm2);
end


if isfield(mri,'psdGram')
    perm1 = [1 2 3 6 5 4]; perm2(perm1) = 1:length(perm1);
    fieldList = {'T' 'W' 't' 'tWin' 'tp' 'vec'};    
    for fieldCur = fieldList
        fieldCur = char(fieldCur);
        if isfield(mri.psdGram,fieldCur)
            sz = size(mri.psdGram.(fieldCur),min(perm1):max(perm1)); sz = sz(perm1); sz(end-1) = prod(sz(end-1:end)); sz(end) = [];
            mri.psdGram.(fieldCur) = permute(reshape(permute(mri.psdGram.(fieldCur),perm1),sz),perm2);
        end
    end
    % perm1 = [1 2 3 6 4 5];
    % perm2 = perm1; perm2(end-1:end) = perm2([end end-1]);
    % fieldList = {'T' 'W' 't' 'tWin' 'tp' 'vec'};    
    % for fieldCur = fieldList
    %     fieldCur = char(fieldCur);
    %     sz = size(mri.psdGram.(fieldCur),min(perm1):max(perm1)); sz = sz(perm1); sz(end-1) = prod(sz(end-1:end)); sz(end) = [];
    %     mri.psdGram.(fieldCur) = permute(reshape(permute(mri.psdGram.(fieldCur),perm1),sz),perm2);
    % end
    tmp = strsplit(mri.psdGram.info,' x ');
    tmp = tmp(perm1); tmp{end-1} = strjoin(tmp(end-1:end),'&'); tmp{end} = '1';
    mri.psdGram.info = strjoin(tmp(perm2),' x ');
end
if isfield(mri,'svdGram')
    perm1 = [1 2 3 6 5 4]; perm2(perm1) = 1:length(perm1);
    fieldList = {'T' 'W' 't' 'tWin' 'tp' 'u' 's' 'v' 'coh'};
    for fieldCur = fieldList
        fieldCur = char(fieldCur);
        if isfield(mri.svdGram,fieldCur)
            sz = size(mri.svdGram.(fieldCur),min(perm1):max(perm1)); sz = sz(perm1); sz(end-1) = prod(sz(end-1:end)); sz(end) = [];
            mri.svdGram.(fieldCur) = permute(reshape(permute(mri.svdGram.(fieldCur),perm1),sz),perm2);
        end
    end
    % perm1 = [1 2 3 6 4 5];
    % perm2 = perm1; perm2(end-1:end) = perm2([end end-1]);
    % fieldList = {'T' 'W' 't' 'tWin' 'tp' 'u' 's' 'v' 'coh'};
    % for fieldCur = fieldList
    %     fieldCur = char(fieldCur);
    %     sz = size(mri.svdGram.(fieldCur),min(perm1):max(perm1)); sz = sz(perm1); sz(end-1) = prod(sz(end-1:end)); sz(end) = [];
    %     mri.svdGram.(fieldCur) = permute(reshape(permute(mri.svdGram.(fieldCur),perm1),sz),perm2);
    % end
    tmp = strsplit(mri.svdGram.info,' x ');
    tmp = tmp(perm1); tmp{end-1} = strjoin(tmp(end-1:end),'&'); tmp{end} = '1';
    mri.svdGram.info = strjoin(tmp(perm2),' x ');
end
