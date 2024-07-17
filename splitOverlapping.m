function mri = splitOverlapping(mri)

%% Sort
if ~issorted([mri.acqTime])
    warning('sorting is altered')
end
[~,b] = sort([mri.acqTime]);
mri = mri(b);

%% Identify overlapping runs (e.g. multi-echo)
acqTime = [mri.acqTime]';
dur = [mri.T]';
indOverl = false(length(acqTime));
for I = 1:length(acqTime)
    Irange = [acqTime(I) acqTime(I)+seconds(dur(I))];
    for i = I:length(acqTime)
        irange = [acqTime(i) acqTime(i)+seconds(dur(i))];
        if ~issorted([Irange irange])
            indOverl(I,i) = true;
        end
    end
end

%% Identify blocks of overlapping runs
indOverl2 = indOverl;
indOverl3 = logical([]);
indOverl4 = logical([]);
while any(diag(indOverl2,1))
    indOverl3 = cat(1,indOverl3,indOverl2(find(diag(indOverl2,1),1),:));
    indOverl2(indOverl3(end,:),:) = false;
    indOverl4 = cat(1,indOverl4,indOverl3(end,:));
    indOverl3(end,find(indOverl3(end,:),1,'first')) = false;
    indOverl4(end,indOverl3(end,:)) = false;
end

%% Rearange the matrix with overlaping runs in dim1 and consecutive runs in dim2
mri = num2cell(mri);
for i = 1:size(indOverl3,1)
    pad = nnz(indOverl3(i,:))+1 - size(mri,1);
    if pad>0
        mri = cat(1,mri,cell(pad,size(mri,2)));
    end
    tmp = mri(1,indOverl3(i,:))'; tmp = [tmp{:}]; [~,b] = sort([tmp.acqTime]); tmp = num2cell(tmp(b))';
    mri(2:nnz(indOverl3(i,:))+1,indOverl4(i,:)) = tmp; mri(1,indOverl3(i,:)) = cell(size(mri(1,indOverl3(i,:)))); clear tmp
end
mri(:,all(cellfun('isempty',mri),1)) = [];
