function volResp = getVolResp2(volTs,volAnat,dsgn,info,force,verbose)
if ~exist('volAnat','var'); volAnat = []; end
if ~exist('dsgn','var');       dsgn = []; end
if ~exist('info','var');       info = []; end
if ~exist('force','var');     force = []; end
if ~exist('verbose','var'); verbose = []; end

if isempty(force);     force = 0; end
if isempty(verbose); verbose = 1; end
        


%% Assert inputs
if isa(volTs,'runCond')
    disp('volTs is runCond format')
end
mList = volTs.fPreprocMaskList(:,1);
for R = 1:length(mList)
    if contains(mList{R},'Inv.nii.gz')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        disp('inverse brain mask provided, trying to replace with brain mask')
        if exist(replace(mList{R},'Inv.nii.gz','.nii.gz'),'file')
            mList{R} = replace(mList{R},'Inv.nii.gz','.nii.gz');
            disp('great success')
        else
            dbstack; error('failed to replace inverse brain mask with brain mask')
        end
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    end
end

fList = volTs.fPreprocList;
if isempty(dsgn)
    dsgn = volTs.dsgn;
end


%% Assert parameters
param.nFrameOrig    = volTs.nFrameOrig;
param.nFrame        = volTs.nFrame;
param.tr            = volTs.tr;
if (abs(dsgn.dt-mean(param.tr))./mean(param.tr))<0.1
    param.trDecon       = dsgn.dt;
else
    param.trDecon       = mean(param.tr);
end
% param.trDecon       = dsgn.dt;



% if mean(tr) == dsgn.dt
%     param.trDecon = mean(tr);
% else
%     param.trDecon = dsgn.dt;
%     warning(strjoin({''...
%         ['volume TR   =  ' sprintf('%7.6f ',mean(tr)) 'sec']...
%         ['stim dt     =  ' sprintf('%7.6f ',dsgn.dt) 'sec']...
%         ['stim onsets = [' sprintf('%7.3f ',dsgn.onsetList) ']sec']...
%         ['            = [' sprintf('%7.3f ',(dsgn.onsetList / mean(tr))) ']vol']...
%         ['Defaulting to stim dt (not TR) for deconvolution = ' num2str(param.trDecon,'%7.6f') 'sec']},newline))
% end




% %% Adjust dsgn
% if isfield(dsgn,'nullTrial') && ~isempty(dsgn.nullTrial)
%     dsgn.cond = dsgn.nullTrial + 1;
%     dsgn.condLabel = {'stim' 'catch'}';
% elseif isempty(dsgn.cond)
%     dsgn.cond = ones(size(dsgn.onsetList));
%     dsgn.condLabel = {'stim'}';
% end





%% Compute response and activation --- magnitude-only data
[fRespCat,fRespRun,fActCat,fActRun] = getRespAndAct2(fList(:,1),dsgn,mList,param,force,verbose);

volResp.respCat = fRespCat;
volResp.respRun = fRespRun;
volResp.actCat = fActCat;
volResp.actRun = fActRun;





%% Compute response --- phase-contrast + magnitude data in complex domain
% Detect phase contrast data
realInd = contains(fList(1,:),'part-real');
imagInd = contains(fList(1,:),'part-imag');
if any(realInd) && any(imagInd)
    disp('phase contrast data detected')
    realInd = contains(fList,'part-real');
    imagInd = contains(fList,'part-imag');
    realInd = find(all(realInd,1));
    imagInd = find(all(imagInd,1));
    if length(realInd)==1 && length(imagInd)==1
        disp('and well-defined')
        param.PCflag = true;
    else
        disp('but not well-defined... skipping')
        param.PCflag = false;
    end
else
    param.PCflag = false;
    clear realInd imagInd
end

% Fit timeseries in complex domain
if param.PCflag
    [fRespCat,fRespRun,fActCat,fActRun] = getRespAndAct2(permute(fList(:,[realInd imagInd]),[1 3 2]),dsgn,mList,param,force,verbose);

    volRespCmplx.respCat = fRespCat;
    volRespCmplx.respRun = fRespRun;
    volRespCmplx.actCat = [];
    volRespCmplx.actRun = [];
end


%% Compute response --- phase-contrast-only data (mag=1) in complex domain
if param.PCflag
    % Remove magnitude data from complex-domain data
    disp('Converting to phase-only complex data (setting magnitude to 1)...');
    fListMag1 = replace(replace(fList(:,[realInd imagInd]),'part-real','part-realMag1'),'part-imag','part-imagMag1');
    for r = 1:size(fListMag1,1)
        disp(['run ' num2str(r) ' of ' num2str(size(fListMag1,1))])
        if ~exist(fileparts(fListMag1{r,1}),'dir'); mkdir(fileparts(fListMag1{r,1})); end
        if ~exist(fileparts(fListMag1{r,2}),'dir'); mkdir(fileparts(fListMag1{r,2})); end
            
        if force || ~exist(fListMag1{r,1},'file') || ~exist(fListMag1{r,2},'file')
            % Read real and imaginary parts
            mriR = MRIread(fList{r,1});
            mriI = MRIread(fList{r,2});
            
            % Convert to polar
            mriPhase = rmfield(mriR,'vol');
            [mriPhase.vol,~] = cart2pol(mriR.vol,mriI.vol);
            mriMag = rmfield(mriR,'vol');
            mriMag.vol = ones(size(mriR.vol));
            
            % Convert back to cartesian
            [mriR.vol,mriI.vol] = pol2cart(mriPhase.vol,mriMag.vol);

            % Write
            MRIwrite(mriR,fListMag1{r,1});
            MRIwrite(mriI,fListMag1{r,2});
        else
            disp('already done, skipping')
        end
    end
    
    % Fit timeseries in complex domain
    [fRespCat,fRespRun,fActCat,fActRun] = getRespAndAct2(permute(fListMag1,[1 3 2]),dsgn,mList,param,force,verbose);

    volRespPC.respCat = fRespCat;
    volRespPC.respRun = fRespRun;
    volRespPC.actCat = [];
    volRespPC.actRun = [];
end











return
    
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


% if any(diff([volTs.nFrame])) || any(diff([volTs.nFrameOrig]))
%     error('not all runs have the same number of frames')
% end
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


