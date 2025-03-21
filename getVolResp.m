function volResp = getVolResp(info,volTs,dsgn,volAnat,force,verbose)

%% Mask
if isstruct(volAnat) && isfield(volAnat,'f')
    fMask = volAnat.f;
    % mask = MRIread(fMask);
    % mask.vol([1:5 end-4:end],[1:5 end-4:end]) = 0;
else
    dbstack; error('code that');
    if isstruct(volAnat) && isfield(volAnat,'f')
        mask = MRIload3(volAnat.f,[],[],0);
        mask = logical(mask.vol);
        mask([1:5 end-4:end],[1:5 end-4:end]) = false;
    else
        dbstack; error('double-check that')
        if ~isempty(volAnat)
            if length(volAnat)==1
                if isMRI(volAnat)
                    mask = volAnat;
                    % if ~isempty(volAnat.vol)
                    %     mask = volAnat.vol;
                    % end
                else
                    if length(volAnat.mask)>1
                        volAnat.mask = volAnat.mask{1};
                    end
                    if isfield(volAnat,'mask') && isfield(volAnat.mask,'crop') && isfield(volAnat.mask.crop,'mri') && isfield(volAnat.mask.crop.mri,'vol') && ~isempty(volAnat.mask.crop.mri.vol)
                        %%%crop
                        mask = volAnat.mask.crop.mri.vol;
                        %%%head
                        mask = mask & any(volAnat.mask.head.mri.vol,4);
                        % %%%brain
                        % mask = mask & any(volAnat.mask.brain.mri.vol,4);
                        %%%apply

                        volTs = applyMask(volTs,mask);
                    else
                        dbstack; error('double-check that')
                    end
                end
            elseif ischar(volAnat)
                mask = MRIload2(MRIload2(volAnat));
            else
                dbstack; error('double-check that')
                volTs = vol2vec(volTs);
                for I = 1:length(volAnat)
                    if isfield(volAnat(I),'fun') && isfield(volAnat(I).fun,'mask') && isfield(volAnat(I).fun.mask,'crop') && ~isempty(volAnat(I).fun.mask.crop.vol)
                        % volTs(I) = applyMask(volTs(I),volAnat(I).fun.mask.crop.vol);

                        %%%crop
                        mask = volAnat(I).fun.mask.crop.vol;
                        if isfield(volAnat(I).fun.mask,'head') && ~isempty(volAnat(I).fun.mask.head)
                            %%%head
                            mask = mask & any(volAnat(I).fun.mask.head.mri.vol,4);
                        elseif isfield(volAnat(I).fun.mask,'brain') && ~isempty(volAnat(I).fun.mask.brain)
                            %%%brain
                            mask = mask & any(volAnat(I).fun.mask.brain.mri.vol,4);
                        end
                        %%%apply
                        volTs(I) = applyMask(volTs(I),mask);
                    end
                end
                % if isfield(volAnat,'fun') && isfield(volAnat.fun,'mask') && isfield(volAnat.fun.mask,'crop') && ~isempty(volAnat.fun.mask.crop.vol)
                %     volTs = applyMask(volTs,volAnat.fun.mask.crop.vol);
                % end
            end
        else
            mask = volTs.vol2vec;
        end
    end
end


if any(diff([volTs.nFrame])) || any(diff([volTs.nFrameOrig]))
    error('not all runs have the same number of frames')
end
param.nDummyRemoved = [volTs.nFrameOrig]' - [volTs.nFrame]';


if ~all(diff([volTs.tr])<0.01); dbstack; error('runs have different tr'); end
tr = [volTs.tr]./1000;


if mean(tr) == dsgn.dt
    param.trDecon = mean(tr);
else
    param.trDecon = dsgn.dt;
    warning(strjoin({''...
        ['volume TR   =  ' sprintf('%7.6f ',mean(tr)) 'sec']...
        ['stim dt     =  ' sprintf('%7.6f ',dsgn.dt) 'sec']...
        ['stim onsets = [' sprintf('%7.3f ',dsgn.onsetList) ']sec']...
        ['            = [' sprintf('%7.3f ',(dsgn.onsetList / mean(tr))) ']vol']...
        ['Defaulting to stim dt (not TR) for deconvolution = ' num2str(param.trDecon,'%7.6f') 'sec']},newline))
end




%% Adjust dsgn
if isfield(dsgn,'nullTrial') && ~isempty(dsgn.nullTrial)
    dsgn.cond = dsgn.nullTrial + 1;
    dsgn.condLabel = {'stim' 'catch'}';
elseif isempty(dsgn.cond)
    dsgn.cond = ones(size(dsgn.onsetList));
    dsgn.condLabel = {'stim'}';
end


%% Compute response and activation
forceThis = force;
verboseThis = verbose;
if isfield(info,'dryRun')
    param.dryRun = info.dryRun;
end
param.nFrame = [volTs.nFrame]';
[fRespCat,fRespRun,fActCat,fActRun] = getRespAndAct(volTs,dsgn,fMask,param,forceThis,verboseThis);
volResp.respCat = fRespCat;
volResp.respRun = fRespRun;
volResp.actCat = fActCat;
volResp.actRun = fActRun;


