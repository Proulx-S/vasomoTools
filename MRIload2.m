function mri = MRIload2(f,volInfo,r)

if isstruct(f) || isa(f,'runCond')
    %this is not a list of files
    if isfield(f,'name') && isfield(f,'folder') && isfield(f,'date') && isfield(f,'bytes') && isfield(f,'isdir') && isfield(f,'datenum')
        %this is the output of dir
        mri = MRIload2(fullfile({f.folder}',{f.name}'),volInfo);
        return
    elseif isMRI(f)
        %this is a header, read the data
        if exist('volInfo','var')
            if ischar(volInfo)
                mask = MRIread(volInfo);
                mask = mask.vol;
            elseif isMRI(volInfo)
                mask = volInfo;
                mask = mask.vol;
            elseif isnumeric(volInfo)
                mask = volInfo;
            end
        else
            mask = [];
        end
        for ii = 1:numel(f)
            if ~isempty(mask)
                mri = vol2vec(MRIread(f(ii).fspec),mask);
            else
                if isfield(f(ii),'vol2vec') && ~isempty(f(ii).vol2vec)
                    mri = vol2vec(MRIread(f(ii).fspec),f(ii).vol2vec);
                else
                    mri = MRIread(f(ii).fspec);
                    mri.imMean = mean(mri.vol,4);
                    % if mask not specified, use all voxels except those with nans or all zeroes
                    zeroMask = all(mri.vol(:,:,:,:)==0,4);
                    nanMask = any(isnan(mri.vol(:,:,:,:)),4);
                    mri.vol2vec = ~zeroMask & ~nanMask;
                    if any(~mri.vol2vec(:))
                        mri.vol2vecFlag = 'validVoxMask';
                    else
                        mri.vol2vecFlag = 'allVoxMask';
                    end
                end
            end
            fieldList = fields(mri);
            fieldList(ismember(fieldList,'nDummy')) = [];
            for i = 1:length(fieldList)
                f(ii).(fieldList{i}) = mri.(fieldList{i});
            end
            if ~isfield(f(ii),'nDummy')
                if isfield(mri,'nDummy')
                    f(ii).nDummy = mri.nDummy;
                else
                    f(ii).nDummy = nan;
                end
            end
        end
        mri = f;
        return
    else
        %this is not a header, read header
        fieldList = fields(f);
        if nnz(ismember(fieldList,'fPreprocList'))==1
            %this contains a list of files, read their headers
            if any(ismember(f.dataType,'volTs'))
                %this is a volTs, read its headers
                if ~isfield(f,'sub')  || isempty(f.sub);  f.sub  = '?';                              end
                if ~isfield(f,'ses')  || isempty(f.ses);  f.ses  = repmat('?',size(f.fPreprocList)); end
                if ~isfield(f,'acq')  || isempty(f.acq);  f.acq  = '?';                              end
                if ~isfield(f,'task') || isempty(f.task); f.task = '?';                              end
                mri = MRIload2(f.fPreprocList,f);
                for r = 1:length(f.fPreprocList)
                    f.volTs(r,1).mri = mri(r);
                end
            else
                %this is not a volTs
                dbstack; error('unrecognized datatype')
            end
            mri = f;
            return
        elseif ismember('mri',fieldList) && isMRI(cat(1,f.mri))
            %this is not a list of files, but contains an mri subfield so read that
            if exist('volInfo','var') && ~isempty(volInfo)
                for i = 1:numel(f)
                    f(i).mri = MRIload2(f(i).mri,volInfo);
                end
            else
                disp('MRIread')
                for i = 1:numel(f)
                    disp([num2str(i) '/' num2str(numel(f))])
                    f(i).mri = MRIload2(f(i).mri);
                end
            end
            mri = f;
            return
        elseif isPHYS(f)
            mri = f;
            return
        elseif isfield(f,'label') && strcmp(f.label,'fs')
            % this is for freesurfer surfaces, skip this
            mri = f;
            return
        elseif isfield(f,'label') && strcmp(f.label,'avMap')
            if ~isfield(f,'fList') || isempty(f.fList)
                mri = f; return
            end
            % this is unprocessed avMap
            if ~isfield(f,'sub')  || isempty(f.sub);  f.sub  = '?';                              end
            if ~isfield(f,'ses')  || isempty(f.ses);  f.ses  = repmat('?',size(f.fPreprocList)); end
            if ~isfield(f,'acq')  || isempty(f.acq);  f.acq  = 'avMap';                          end
            % if ~isfield(f,'task') || isempty(f.task); f.task = '?';                              end
            if size(f.ses,1)==1 && size(f.fList,1)~=1
                f.ses = repmat(f.ses,size(f.fList));
            end
            mri = MRIload2(f.fList,f);
            return
        else
            %this is not a list of files and does not contain an mri subfield, so loop over subfields
            for i = 1:length(fieldList)
                % disp(fieldList{i})
                if isempty(f.(fieldList{i})); continue; end
                f.(fieldList{i}) = MRIload2(f.(fieldList{i}));
            end
            mri = f;
            return
        end
    end
elseif iscell(f)
    for r = 1:length(f)
        if ischar(f{r})
            %this is a list of files, read their header and copy info
            if isa(volInfo,'runCond')
                label1 = strjoin({...
                    ['sub-'  volInfo.sub]
                    ['ses-'  volInfo.ses(r,:)]
                    ['acq-'  volInfo.acq]
                    ['task-' volInfo.task]...
                    },'_');
            elseif isfield(volInfo,'labelAcq')
                label1 = strjoin({...
                    ['sub-'  volInfo.sub]
                    ['ses-'  volInfo.ses(r,:)]
                    ['acq-'  volInfo.labelAcq]
                    ['task-' volInfo.label]...
                    },'_');
            else
                label1 = strjoin({...
                    ['sub-'  volInfo.sub]
                    ['ses-'  volInfo.ses(r,:)]
                    ['acq-'  volInfo.label]
                    ['task-' 'none']...
                    },'_');
            end
            if isfield(volInfo,'fPreprocList') && ~isempty(volInfo.fPreprocList)
                label2 = [num2str(r) '/' num2str(length(volInfo.fPreprocList))];
            else
                label2 = [num2str(r) '/' num2str(length(volInfo.fList))];
            end
            disp(['   ' label1 '; ' label2])
            % disp([' ' num2str(r) '/' num2str(length(f))]);
            mri(r,1) = MRIload2(f{r},volInfo,r);
        else
            %this not a list of files, loop over cells
            disp(['---dataset ' num2str(r) '/' num2str(length(f))]);
            mri{r} = MRIload2(f{r});
        end
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
    if isfield(volInfo,'nFrame') && ~isempty(volInfo.nFrame)
        nFrame = volInfo.nFrame(r);
    else
        mri.nFrame = nan;
    end
    if isfield(volInfo,'fOrigList')
        nFrameOrig = MRIread(volInfo.fOrigList{r},1); nFrameOrig = nFrameOrig.nframes;
        mri.nDummyRemoved = nFrameOrig-nFrame;
    else
        mri.nDummyRemoved = nan;
    end
else
    if ~isfield(mri,'nDummyRemoved')
        mri.nDummyRemoved = [];
    end
end
