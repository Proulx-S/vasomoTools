function [volAnat,rCond] = volAnatPreproc6(rCond,force,verbose)
    global src
    if ~exist('force','var');     force = []; end
    if ~exist('verbose','var'); verbose = []; end
    if isempty(force);            force = 0; end
    if isempty(verbose);        verbose = 0; end

        

    %% Summarize across runs

    %%% Combine tasks
    [fList,fMaskList,nDummy,taskList,acqTime] = combineRunsAcrossTasks(rCond);
    taskList = unique(taskList);
    % taskList = fields(rCond); taskList(~contains(taskList,'task_')) = [];
    % if length(taskList)>1; dbstack; error('need to combine data from multiple files'); end
    
    sub     = rCond.(taskList{1}).sub;
    acq     = rCond.(taskList{1}).acq;
    prsc    = rCond.(taskList{1}).prsc;
    vencAcq = rCond.(taskList{1}).vencAcq;
    disp('--------------------------------')
    disp(['volAnat: ' rCond.(taskList{1}).sub ' acq-' acq ' prsc-' prsc ' venc-' vencAcq])
    disp('--------------------------------')
    acqLabel = ['acq-' acq '_prsc-' prsc '_venc-' vencAcq];


    % volAnat in all subfields should now be exactly the same
    volAnat = rCond.(taskList{1}).volAnat;


    % %%% Combine runs
    % fList = cell(size(taskList));
    % for T = 1:length(taskList)
    %     fList{T} = rCond.(char(taskList{T})).fPreprocList;
    % end
    % fList = cat(1,fList{:});
    
    fAvList = cell(size(fList));
    for R = 1:size(fList,1)
        fAvList{R} = strsplit(fList{R,1,1},filesep); fAvList{R}{end} = ['av_' fAvList{R}{end}]; fAvList{R} = strjoin(fAvList{R},filesep);
    end

    fName = strsplit(fAvList{1},filesep); fName = fName{end};
    fCatAv = unique(fileparts(fileparts(fAvList)));
    if length(fCatAv)>1
        fCatAv = strsplit(fCatAv{1},filesep);
        fCatAv{contains(fCatAv,'ses-')} = 'ses-cat';
        fCatAv = strjoin(fCatAv,filesep);
    else
        fCatAv = char(fCatAv);
    end
    fCatAv = fullfile(fCatAv,['cat_' fName]);
    fAvCatAv = strsplit(fCatAv,filesep); fAvCatAv{end} = ['av_' fAvCatAv{end}]; fAvCatAv = strjoin(fAvCatAv,filesep);
    
    if ~exist(fileparts(fAvCatAv),'dir'); mkdir(fileparts(fAvCatAv)); end

    % fAvList = cell(size(rCond.(char(taskList)).fPreprocList,1),1);
    % for R = 1:size(rCond.(char(taskList)).fPreprocList,1)
    %     fAvList{R} = strsplit(rCond.(char(taskList)).fPreprocList{R,1,1},filesep); fAvList{R}{end} = ['av_' fAvList{R}{end}]; fAvList{R} = strjoin(fAvList{R},filesep);
    % end
    % fCatAv = unique(fileparts(fileparts(fAvList)));
    % if length(fCatAv)>1; dbstack; error('not sure where to store this'); end
    % fCatAv = char(fCatAv);
    % [~,b,~] = fileparts(replace(fAvList,'.nii.gz','')); b = unique(b);
    % if length(b)>1; dbstack; error('not sure who to name this'); end
    % fCatAv = fullfile(fCatAv,['cat_' char(b) '.nii.gz']);
    % fAvCatAv = strsplit(fCatAv,filesep); fAvCatAv{end} = ['av_' fAvCatAv{end}]; fAvCatAv = strjoin(fAvCatAv,filesep);

    if force || ~exist(fCatAv,'file') || ~exist(fAvCatAv,'file')
        cmd = {src.afni};
        cmd{end+1} = ['3dTcat -overwrite \'];
        cmd{end+1} = ['-prefix ' fCatAv ' \'];
        cmd{end+1} = strjoin(fAvList,' ');
        cmd{end+1} = ['3dTstat -overwrite -mean \'];
        cmd{end+1} = ['-prefix ' fAvCatAv ' \'];
        cmd{end+1} = fCatAv;
        if verbose
            [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
        else
            [status,cmdout] = system(strjoin(cmd,newline)        ); if status; dbstack; error(cmdout); error('x'); end
        end
    end



    %%% Combine masks
    fMaskBrainInv = unique(fMaskList); fMaskBrainInv(cellfun('isempty',fMaskBrainInv)) = [];
    maskDir = strsplit(rCond.(taskList{1}).dirs.bidsDeriv,filesep);
    ind = contains(maskDir,'sub-');
    if nnz(ind)>1; dbstack; error('more than one mask found'); end
    maskDir = strjoin([maskDir(1:find(ind)) {acqLabel 'anat'}],filesep);
    if ~exist(maskDir,'dir'); mkdir(maskDir); end
    
    mriOblq = MRIread(fAvCatAv,1); mriOblq.vol = false(mriOblq.volsize);
    for i = 1:length(fMaskBrainInv)
        mri = MRIread(fMaskBrainInv{i});
        mriOblq.vol = mriOblq.vol | mri.vol;
    end

    fMaskBrainInv = fullfile(maskDir,'brainMaskInv.nii.gz');
    fMaskBrain    = fullfile(maskDir,'brainMask.nii.gz');
    MRIwrite(mriOblq,fMaskBrainInv);
    mriOblq.vol = mriOblq.vol==0;
    MRIwrite(mriOblq,fMaskBrain);

    


    % if length(fMaskBrainInv)>1; dbstack; error('more than one mask found'); end; fMaskBrainInv = char(fMaskBrainInv);
    % %%%% conform
    % fMaskBrainInv2 = replace(fAvCatAv,'_volTs.nii.gz','_brainMaskInv.nii.gz'); copyfile(fMaskBrainInv,fMaskBrainInv2);
    % fMaskBrainInv = fMaskBrainInv2;
    % MRIconform(fMaskBrainInv,fAvCatAv);
    % fMaskBrain = replace(fMaskBrainInv,'_brainMaskInv.nii.gz','_brainMask.nii.gz');
    % if force || ~exist(fMaskBrain,'file')
    %     cmd = {src.afni};
    %     cmd{end+1} = ['3dcalc -overwrite -a ' fMaskBrainInv ' -expr ''-(a-1)'' -prefix ' fMaskBrain];
    %     if verbose
    %         [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    %     else
    %         [status,cmdout] = system(strjoin(cmd,newline)        ); if status; dbstack; error(cmdout); error('x'); end
    %     end
    % end
    % MRIconform(fMaskBrain,fAvCatAv);



    %%% Output file index
    volAnat.mask.brain.f     = fMaskBrain;
    volAnat.mask.brain.fInv  = fMaskBrainInv;
    volAnat.mask.brain.fBase = fAvCatAv;
    

    % task = char(taskList);
    % acq = rCond.(task).acq;
    % prsc = rCond.(task).prsc;
    

    %% Preprocess anat for vesselness map (for a starting point to vessel drawing)
    forceThis = force;
    %%%% correct bias field
    fVolCorr   = cell(size(fAvList));
    for r = 1:length(fAvList)
        fVolCorr{r} = correctBiasField(fAvList{r}, volAnat.mask.brain.f, fAvList{r}, [], forceThis, verbose);
    end
    
    forceThis = force;
    %%%% compute vesselness
    [fComp,fNonComp,fSegFig] = computeVesselness(fVolCorr,fMaskBrain,forceThis,verbose);
    

    %%%% summarize vesselness across runs
    derivDir = fullfile(rCond.(taskList{1}).dirsOrig.bidsDeriv,acqLabel);
    if ~exist(derivDir,'dir'); mkdir(derivDir); end
    
    fVesselness = fullfile(derivDir,'vesselness.nii.gz');
    mri = MRIread(fComp{r,1},1); mri.vol = zeros(mri.volsize);
    for r = 1:size(fComp,1)
        mri.vol = mri.vol + MRIread(fComp{r,1}).vol;
    end
    mri.vol = mri.vol/length(fComp);
    MRIwrite(mri,fVesselness);
    


    

    %% Process other anatomical images
    


    % % %%% find volAnat in rCond
    % % volAnat  = cell(size(taskList));
    % % nonEmpty = false(size(taskList));
    % % for T = 1:length(taskList)
    % %     volAnat{T} = rCond.(taskList{T}).volAnat;
    % %     acqAnatList = fields(volAnat{T});
    % %     tmp = false;
    % %     for A = 1:length(acqAnatList)
    % %         tmp = tmp | numel(volAnat{T}.(acqAnatList{A}))>0;
    % %     end
    % %     nonEmpty(T) = tmp;
    % % end
    % % if nnz(nonEmpty)>1; dbstack; error('more than one volAnat found'); end
    % % volAnat = rCond.(taskList{nonEmpty}).volAnat;
    % taskList2 = fields(rCond); taskList2 = taskList2(contains(taskList2,'task_'));
    % volAnatTmp = [];
    % for T = 1:length(taskList2)
    %     % acqAnatList = fields(rCond.(taskList2{T}).volAnat);
    %     volAnatTmpTmp = rCond.(taskList2{T}).volAnat;
    %     if isfield(volAnatTmpTmp,'label'); volAnatTmpTmp = rmfield(volAnatTmpTmp,'label'); end
    %     if isfield(volAnatTmpTmp,'roi');   volAnatTmpTmp = rmfield(volAnatTmpTmp,'roi');   end
    %     volAnatTmp = cat(1,volAnatTmp,volAnatTmpTmp);
    % end
    % % combine volAnat
    % imList = fields(volAnatTmp);
    % volAnat = [];
    % for i = 1:length(imList)
    %     volAnat.(imList{i}) = cat(1,volAnatTmp.(imList{i}));
    %     [~,b,~] = unique(fullfile({volAnat.(imList{i}).folder},{volAnat.(imList{i}).name}));
    %     volAnat.(imList{i}) = volAnat.(imList{i})(b);
    % end
    % % redistribute volAnat
    % for T = 1:length(taskList2)
    %     rCond.(taskList2{T}).volAnat = volAnat;
    % end




    %%% find and sort avMap
    fAvMapEchoCat     = fullfile(derivDir,'avMapEchoCat.nii.gz');
    fAvMapEchoCat_reg = fullfile(derivDir,'avMapEchoCat.lta');
    fTmp = fullfile({volAnat.avMap.folder},{volAnat.avMap.name});
    [a,b,c] = unique(regexprep(fTmp, '_echo-\d+_', '_'));
    if length(a)>1
        % pick one
        switch sub
            case {'vsmDiamCenSurP1' 'vsmDiamCenSurP7' 'vsmDiamCenSurP10'}
                f = fTmp(c==2);
            case {'vsmDiamCenSurP2' 'vsmDiamCenSurP4' 'vsmDiamCenSurP5'}
                f = fTmp(c==1);
            otherwise
                fTmp
                dbstack; keyboard;
        end
    else
        f = fTmp;
    end
    mriOut = MRIread(f{1},1);
    for E = 1:length(f)
        mri = MRIread(f{E});
        mriOut.vol = cat(4,mriOut.vol,mri.vol);
    end
    MRIwrite(mriOut,fAvMapEchoCat);
    % Update header information of fAvMapEchoCat to match fAvCatAv without interpolation
    fAvMapEchoCatCnfrm = fullfile(derivDir,'avMapEchoCat_conform.nii.gz');
    cmd = {src.fs};
    cmd{end+1} = ['mri_copy_params ' fAvMapEchoCat ' ' fAvCatAv ' ' fAvMapEchoCatCnfrm];
    disp('Updating header information of echo cat volume to match average...');
    [status, result] = system(strjoin(cmd, newline), '-echo');
    if status ~= 0
        warning('Failed to update header information: %s', result);
    end
    
    %%% find tof
    fTof     = fullfile(derivDir,'tof.nii.gz');
    fTof_reg = fullfile(derivDir,'tof.lta');
    f = fullfile({volAnat.tof.folder},{volAnat.tof.name});
    % just pick the first one if more than one
    copyfile(f{1},fTof);
    
    %%% find memprage
    fMemprage     = fullfile(derivDir,'memprage.nii.gz');
    fMemprage_reg = fullfile(derivDir,'memprage.lta');
    f = fullfile({volAnat.memprage.folder},{volAnat.memprage.name});
    f = f(contains(f,'proc-rms','IgnoreCase',true));
    % just pick the first one if more than one
    if isempty(f)
        fMemprage = []; fMemprage_reg = [];
    else
        copyfile(f{1},fMemprage);
    end


    %%% memprage N4 bias correction
    if ~isempty(fMemprage)
        fMemprageBiasCorr = fullfile(derivDir, 'memprage_n4.nii.gz');
        if force || ~exist(fMemprageBiasCorr, 'file')
            cmd = {src.ants};
            cmd{end+1} = ['N4BiasFieldCorrection -d 3 -i ' fMemprage ' -o ' fMemprageBiasCorr];
            disp('Running ANTs N4 bias field correction...');
            [status, result] = system(strjoin(cmd,newline),'-echo');
            if status ~= 0
                warning('N4 bias correction failed with message: %s', result);
                % Fall back to original image if bias correction fails
                fMemprageBiasCorr = fMemprage;
                disp('Using original image for brain extraction due to N4 failure.');
            else
                disp(['Bias corrected file created: ' fMemprageBiasCorr]);
            end
        else
            disp(['Bias corrected file already exists: ' fMemprageBiasCorr]);
        end
        
        %%% memprage brain extraction using FreeSurfer's machine learning based tools
        fMemprageBrain = fullfile(derivDir, 'memprage_brain.nii.gz');
        fMemprageBrainMask = fullfile(derivDir, 'memprage_brain_mask.nii.gz');
        if force || ~exist(fMemprageBrain, 'file')
            cmd = {src.fs};
            cmd{end+1} = ['mri_synthstrip -i ' fMemprageBiasCorr ' -o ' fMemprageBrain ' -m ' fMemprageBrainMask];
            disp('Running FreeSurfer SynthStrip brain extraction (machine learning based)...');
            [status, result] = system(strjoin(cmd,newline),'-echo');
            if status ~= 0
                warning('FreeSurfer brain extraction failed with message: %s', result);
            else
                disp(['Brain extracted file created: ' fMemprageBrain]);
                disp(['Brain mask file created: ' fMemprageBrainMask]);
            end
        else
            disp(['Brain extracted file already exists: ' fMemprageBrain]);
        end

        %%% interpolate memprage brain mask to TOF space
        fMemprageBrainMaskInTofSpace = fullfile(derivDir, 'memprage_brain_mask_tof_space.nii.gz');
        if force || ~exist(fMemprageBrainMaskInTofSpace, 'file')
            disp('Interpolating brain mask to TOF space...');
            cmd = {src.fs};
            cmd{end+1} = ['mri_vol2vol --mov ' fMemprageBrainMask ' --targ ' fTof ' --regheader --o ' fMemprageBrainMaskInTofSpace ' --nearest'];
            [status, result] = system(strjoin(cmd,newline),'-echo');
            if status ~= 0
                warning('Brain mask interpolation to TOF space failed with message: %s', result);
            else
                disp(['Brain mask interpolated to TOF space: ' fMemprageBrainMaskInTofSpace]);
            end
        else
            disp(['Brain mask in TOF space already exists: ' fMemprageBrainMaskInTofSpace]);
        end

        %%% tof N4 bias correction
        fTofBiasCorr = fullfile(derivDir, 'tof_n4.nii.gz');
        if force || ~exist(fTofBiasCorr, 'file')
            cmd = {src.ants};
            cmd{end+1} = ['N4BiasFieldCorrection -d 3 -i ' fTof ' -o ' fTofBiasCorr ' -x ' fMemprageBrainMaskInTofSpace];
            disp('Running ANTs N4 bias field correction on TOF image...');
            [status, result] = system(strjoin(cmd,newline),'-echo');
            if status ~= 0
                warning('N4 bias correction for TOF failed with message: %s', result);
                % Fall back to original image if bias correction fails
                fTofBiasCorr = fTof;
                disp('Using original TOF image due to N4 failure.');
            else
                disp(['Bias corrected TOF file created: ' fTofBiasCorr]);
            end
        else
            disp(['Bias corrected TOF file already exists: ' fTofBiasCorr]);
        end
    
        %%% tof brain masking
        fTofBiasCorrBrain = fullfile(derivDir, 'tof_n4_brain.nii.gz');
        if force || ~exist(fTofBiasCorrBrain, 'file')
            disp('Applying brain mask to TOF image using FreeSurfer mri_mask...');
            cmd = {src.fs};
            % Use mri_mask to apply the brain mask to the TOF image
            cmd{end+1} = ['mri_mask ' fTofBiasCorr ' ' fMemprageBrainMask ' ' fTofBiasCorrBrain];
            
            % Execute the FreeSurfer command
            [status, result] = system(strjoin(cmd,newline),'-echo');
            if status ~= 0
                warning('FreeSurfer brain masking failed with message: %s', result);
            else
                disp(['Brain-extracted TOF created: ' fTofBiasCorrBrain]);
            end
        else
            disp(['Brain-extracted TOF already exists: ' fTofBiasCorrBrain]);
        end
    else
        fTofBiasCorr = fTof;
        fTofBiasCorrBrain = fTof;
    end

    





    % %%% Create average of functional data
    % fAvCatAv = fullfile(derivDir, 'avCatAv.nii.gz');
    % fCatAv = fullfile(derivDir, 'catAv.nii.gz');
    
    % This section assumes these files are created elsewhere
    % If they don't exist, we should handle that case
    if ~exist(fAvCatAv, 'file') || ~exist(fCatAv, 'file')
        warning('Average functional data files not found. These should be created before this point.');
    end
    



    %% Draw vessel rois
    %%% find earlier vessel rois file in derivatives
    fVesselRoi_deriv = {};
    for T = 1:length(taskList)
        fVesselRoi_deriv{end+1} = dir(fullfile(rCond.(taskList{T}).dirsOrig.bidsDeriv,'vessel.nii.gz'));
        fVesselRoi_deriv{end+1} = dir(fullfile(rCond.(taskList{T}).dirsOrig.bidsDeriv,'*','vessel.nii.gz'));
        fVesselRoi_deriv{end+1} = dir(fullfile(rCond.(taskList{T}).dirsOrig.bidsDeriv,'*','*','vessel.nii.gz'));
    end
    fVesselRoi_deriv(cellfun('isempty',fVesselRoi_deriv)) = [];

    %%% copy vessel roi from derivatives or use vesselness
    fVesselRoi = fullfile(maskDir,'vessel.nii.gz');
    if ~isempty(fVesselRoi_deriv)
        fVesselRoi_deriv = fVesselRoi_deriv{1}; fVesselRoi_deriv = fullfile(fVesselRoi_deriv.folder,fVesselRoi_deriv.name);
        copyfile(fVesselRoi_deriv,fVesselRoi);
    else
        fVesselRoi_deriv  = fullfile(derivDir,'vessel.nii.gz'); 
        copyfile(fVesselness,fVesselRoi);
    end
    % fVesselRoiMaxVox  = fullfile(derivDir,'vesselMaxVox.nii.gz'); copyfile(fVesselRoi,fVesselRoiMaxVox);
    % fVesselRoiDilated = fullfile(derivDir,'vesselDilated.nii.gz'); copyfile(fVesselRoi,fVesselRoiDilated);
    


    %%% Draw in neurodesk
    fNeurodeskList = {
        fMemprage
        fTofBiasCorrBrain
        fAvMapEchoCatCnfrm
        fAvCatAv
        fCatAv
        % fVesselRoiDilated
        fVesselRoi
        % fVesselRoiMaxVox
        fVesselness
        };
    tmpDir = tempname; mkdir(tmpDir);
    for i = 1:length(fNeurodeskList)
        if ~isempty(fNeurodeskList{i})
            copyfile(fNeurodeskList{i},tmpDir);
            if exist(replace(fNeurodeskList{i},'.nii.gz','.lta'),'file')
                copyfile(replace(fNeurodeskList{i},'.nii.gz','.lta'),tmpDir);
            end
        end
    end

    if ~isempty(fMemprage)
        [~,fMemprage2,~]          = fileparts(replace(fMemprage         ,'.nii.gz','')); fMemprage2          = [fMemprage2,'.nii.gz'];
    else
        fMemprage2 = [];
    end
    [~,fTofBiasCorrBrain2,~]  = fileparts(replace(fTofBiasCorrBrain ,'.nii.gz','')); fTofBiasCorrBrain2  = [fTofBiasCorrBrain2,'.nii.gz'];
    [~,fAvMapEchoCat2,~]      = fileparts(replace(fAvMapEchoCat     ,'.nii.gz','')); fAvMapEchoCat2      = [fAvMapEchoCat2,'.nii.gz'];
    [~,fAvMapEchoCatCnfrm2,~] = fileparts(replace(fAvMapEchoCatCnfrm,'.nii.gz','')); fAvMapEchoCatCnfrm2 = [fAvMapEchoCatCnfrm2,'.nii.gz'];
    [~,fAvCatAv2,~]           = fileparts(replace(fAvCatAv          ,'.nii.gz','')); fAvCatAv2           = [fAvCatAv2,'.nii.gz'];
    [~,fCatAv2  ,~]           = fileparts(replace(fCatAv            ,'.nii.gz','')); fCatAv2             = [fCatAv2,'.nii.gz'];
    [~,fVesselRoi2,~]         = fileparts(replace(fVesselRoi        ,'.nii.gz','')); fVesselRoi2         = [fVesselRoi2,'.nii.gz'];
    [~,fVesselness2,~]        = fileparts(replace(fVesselness       ,'.nii.gz','')); fVesselness2        = [fVesselness2,'.nii.gz'];
    if ~isempty(fMemprage)
        if exist(replace(fMemprage,'.nii.gz','.lta'),'file');          fMemprage3          = [fMemprage2 ':reg=' replace(fMemprage2,'.nii.gz','.lta')];                   else fMemprage3          = fMemprage2; end
    else
        fMemprage3 = [];
    end
    if exist(replace(fTofBiasCorrBrain,'.nii.gz','.lta'),'file');  fTofBiasCorrBrain3  = [fTofBiasCorrBrain2 ':reg=' replace(fTofBiasCorrBrain2,'.nii.gz','.lta')];   else fTofBiasCorrBrain3  = fTofBiasCorrBrain2; end
    if exist(replace(fAvMapEchoCat,'.nii.gz','.lta'),'file');      fAvMapEchoCat3      = [fAvMapEchoCat2 ':reg=' replace(fAvMapEchoCat2,'.nii.gz','.lta')];           else fAvMapEchoCat3      = fAvMapEchoCat2; end            
    if exist(replace(fAvMapEchoCatCnfrm,'.nii.gz','.lta'),'file'); fAvMapEchoCatCnfrm3 = [fAvMapEchoCatCnfrm2 ':reg=' replace(fAvMapEchoCatCnfrm2,'.nii.gz','.lta')]; else fAvMapEchoCatCnfrm3 = fAvMapEchoCatCnfrm2; end
    if exist(replace(fAvCatAv,'.nii.gz','.lta'),'file');           fAvCatAv3           = [fAvCatAv2 ':reg=' replace(fAvCatAv2,'.nii.gz','.lta')];                     else fAvCatAv3           = fAvCatAv2; end
    if exist(replace(fCatAv,'.nii.gz','.lta'),'file');             fCatAv3             = [fCatAv2 ':reg=' replace(fCatAv2,'.nii.gz','.lta')];                         else fCatAv3             = fCatAv2; end
    if exist(replace(fVesselRoi,'.nii.gz','.lta'),'file');         fVesselRoi3         = [fVesselRoi2 ':reg=' replace(fVesselRoi2,'.nii.gz','.lta')];                 else fVesselRoi3         = fVesselRoi2; end
    if exist(replace(fVesselness,'.nii.gz','.lta'),'file');        fVesselness3        = [fVesselness2 ':reg=' replace(fVesselness2,'.nii.gz','.lta')];               else fVesselness3        = fVesselness2; end

    
    %%% To draw masks, copy to neurocloud, use freeview there, then copy back
    subDir = ['sub-' sub '_' acqLabel];
    cmd = {src.fs};
    cmd{end+1} = ['mkdir -p ~/' subDir];
    cmd{end+1} = ['cd ~/' subDir];
    cmd{end+1} = ['rsync sebp@takoyaki1:' tmpDir '/*{.nii.gz,.lta} .'];
    cmd{end+1} = 'echo draw VESSEL ROIs (aretery=902, vein=914, left-vessel=30, right-vessel=62)';
    cmd{end+1} = 'freeview -v \';
    cmd{end+1} = [fCatAv3             ':grayscale=0,1000 \'];
    if ~isempty(fMemprage)
        cmd{end+1} = [fMemprage3          ':resample=cubic \'];
    end
    cmd{end+1} = [fTofBiasCorrBrain3  ':resample=cubic \'];
    cmd{end+1} = [fAvMapEchoCatCnfrm3 ':resample=cubic \'];
    cmd{end+1} = [fAvCatAv3           ':grayscale=0,1000 \'];
    cmd{end+1} = [fVesselness3        ':colormap=heat \'];
    cmd{end+1} = [fVesselRoi3         ':colormap=lut'];
    cmd{end+1} = ['rsync ./*{.nii.gz,.lta,.vtk} sebp@takoyaki1:' tmpDir '/'];

    cmdFile = [tmpDir '.sh'];
    fileID = fopen(cmdFile, 'w');
    fprintf(fileID, '%s\n', cmd{:});
    fclose(fileID);
    disp('++++++++++++++++++++++++++++++++++++++++')
    disp('Commands for vessel roi drawing are in file:')
    disp(cmdFile)
    disp('Paste in freeview capable remote to transfer data, create masks and transfer back.')
    %%% Wait for user to confirm mask drawing is done
    if force || ~exist(fVesselRoi_deriv,'file')
        keyboard
        done = '';
        while ~strcmpi(done, 'done')
            disp('When done, type "done"')
            done = input('', 's');
        end
        disp('++++++++++++++++++++++++++++++++++++++++')

        % move out of tmpDir
        copyfile(fullfile(tmpDir,fVesselRoi2),fVesselRoi)
        if exist(fullfile(tmpDir,replace(fAvMapEchoCatCnfrm2,'.nii.gz','.lta')),'file')
            copyfile(fullfile(tmpDir,replace(fAvMapEchoCatCnfrm2,'.nii.gz','.lta')),replace(fAvMapEchoCatCnfrm,'.nii.gz','.lta'));
        end
        for i = 1:length(fNeurodeskList)
            if ~isempty(fNeurodeskList{i})
                [~,fName,~] = fileparts(replace(fNeurodeskList{i},'.nii.gz',''));
                if exist(fullfile(tmpDir,[fName '.lta']),'file')
                    copyfile(fullfile(tmpDir,[fName '.lta']),replace(fNeurodeskList{i},'.nii.gz','.lta'));
                end
                if exist(fullfile(tmpDir,[fName '.vtk']),'file')
                    copyfile(fullfile(tmpDir,[fName '.vtk']),replace(fNeurodeskList{i},'.nii.gz','.vtk'));
                end
            end
        end



        %%% Store to permanent bids derivatives
        copyfile(fVesselRoi,fVesselRoi_deriv)
    end
    


    %%% Output file index
    volAnat.label.calcarineVessel.f         = fVesselRoi;
    volAnat.label.calcarineVessel.label     = {'artery' 'vein' 'Left-vessel' 'Right-vessel'};
    volAnat.label.calcarineVessel.labelVal  = [902 914 30 62];
    volAnat.label.calcarineVessel.fBase     = fAvCatAv;
    if ~isempty(fMemprage)
        volAnat.label.calcarineVessel.fBaseList = {fCatAv fVesselness fAvMapEchoCatCnfrm fTofBiasCorrBrain fMemprage};
    else
        volAnat.label.calcarineVessel.fBaseList = {fCatAv fVesselness fAvMapEchoCatCnfrm fTofBiasCorrBrain};
    end
    volAnat.label.calcarineVessel.fFig      = fSegFig;
    



    

    %% Individual vessel ROIs
    if exist(fVesselRoi,'file')
        label = volAnat.label.calcarineVessel;
        imField = {'base' 'vesselness'};
        [~,b,~] = fileparts(label.fBaseList);
        im = {label.fBase label.fBaseList{contains(b,'vesselness.nii')}};
        cropSz = 10;
        volAnat.roi.vessel = getVesselRoi2(label,imField,im,cropSz);
    else
        % keyboard
        volAnat.roi.vessel = [];
    end
    

    if verbose>1 && exist('label','var')
        ax = {};

        mri = MRIread(label.f);
        mri.vol = mri.vol~=0;
        figure('WindowStyle','docked');
        imagesc(mri.vol(:,:,:,1));
        colormap(gray);
        axis image; axis off;
        ax{end+1} = gca;
        
        mri = MRIread(label.fBase);
        figure('WindowStyle','docked');
        imagesc(mri.vol(:,:,:,1));
        colormap(gray);
        axis image; axis off;
        hold on;
        ax{end+1} = gca;

        ind = ismember({volAnat.roi.vessel.class},'artery');
        hP = plot([volAnat.roi.vessel(ind).poly]);
        set(hP,'EdgeColor','r');
        set(hP,'FaceColor','none');
        ind = ismember({volAnat.roi.vessel.class},'vein');
        hP = plot([volAnat.roi.vessel(ind).poly]);
        set(hP,'EdgeColor','b');
        set(hP,'FaceColor','none');
        ind = ismember({volAnat.roi.vessel.class},'unknown');
        hP = plot([volAnat.roi.vessel(ind).poly]);
        set(hP,'EdgeColor','y');
        set(hP,'FaceColor','none');
        axis image; axis off;
        ax{end+1} = gca;
        linkaxes([ax{:}]);
    end




    %% Insert into rCond
    if nargout>1
        for T = 1:length(taskList)
            rCond.(char(taskList{T})).volAnat = volAnat;
        end
    end



    %%%%%
    disp('!!!!!!!!!!!!!!!!!!!!!!!!')
    disp('------------------------')
    disp('manually backup drawing')
    bkDir = '/local/users/Proulx-S/db/vsmDiamCenSur';
    if ~exist(fullfile(bkDir,subDir),'dir'); mkdir(fullfile(bkDir,subDir)); end
    disp(['cp ' fVesselRoi_deriv ' ' [fullfile(bkDir,subDir) filesep]])
    ls(fullfile(bkDir,'*','*'))
    disp('------------------------')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!')





function drawVesselRoi(fBaseList,fVesselRoi,fVesselMask,fBrainMask,fAvMap,fTof)



    if forceThis || ~exist(runSet{S}{rs}.fMasks.fMaskInv,'file')
        mri     = MRIread(fBase,1);
        mri.vol = ones(mri.volsize);
        MRIwrite(mri,runSet{S}{rs}.fMasks.fMaskInv);
        cmd{end+1} = ['scp sebp@takoyaki1:' runSet{S}{rs}.fMasks.fMaskInv ' sebp@takoyaki1:' fBase ' .'];
        cmd{end+1} = 'echo draw EXCLUSION mask for the BRAIN (brain=0, nonBrain=1)';
        cmd{end+1} = 'freeview -v \';
        cmd{end+1} = [replace(fBase,[fileparts(fBase) filesep],'./') ' \'];
        cmd{end+1} = [replace(runSet{S}{rs}.fMasks.fMaskInv,[fileparts(runSet{S}{rs}.fMasks.fMaskInv) filesep],'./') ':colormap=heat:opacity=0.33'];
        cmd{end+1} = ['scp ' replace(runSet{S}{rs}.fMasks.fMaskInv,[fileparts(runSet{S}{rs}.fMasks.fMaskInv) filesep],'./') ' sebp@takoyaki1:' runSet{S}{rs}.fMasks.fMaskInv ''];
    end



%     cmd = {src.afni};

% tmpDir = fullfile('.',subList{S},'venc14_cmplxMag1');
% cmd{end+1} = ['rm -fr ' tmpDir];
% cmd{end+1} = ['mkdir -p ' tmpDir];
% cmd{end+1} = ['rsync sebp@takoyaki1:{' strjoin(...
% {[char(rCond{S}.vfMRIpc_dflt_pcVenc14ap.task_50sPrd5sDur.volResp.cmplxMag1.respCat.fStat) '+orig.*']
% rCond{S}.vfMRIpc_dflt_pcVenc14ap.task_50sPrd5sDur.volResp.cmplxMag1.respCat.stats.fResp{1,1,3}},...
%     ',') '} ' tmpDir '/'];
% cmd{end+1} = ['cd ' tmpDir];


% tmpDir = fullfile('.',subList{S},'venc14_cmplx');
% cmd{end+1} = ['rm -fr ' tmpDir];
% cmd{end+1} = ['mkdir -p ' tmpDir];
% cmd{end+1} = ['rsync sebp@takoyaki1:{' strjoin(...
% {[char(rCond{S}.vfMRIpc_dflt_pcVenc14ap.task_50sPrd5sDur.volResp.cmplx.respCat.fStat) '+orig.*']
% rCond{S}.vfMRIpc_dflt_pcVenc14ap.task_50sPrd5sDur.volResp.cmplx.respCat.stats.fResp{1,1,3}},...
%     ',') '} ' tmpDir '/'];
% cmd{end+1} = ['cd ' tmpDir];

% disp(strjoin(cmd,newline))



    %% Individual vessel ROIs
    [a,b,~] = fileparts(replace(volAnat.mask.vesselRefined.f,'.nii.gz',''));
    volAnat.label.calcarineVessel.f = fullfile(a,'calcarineVessel.nii.gz');

    forceThis = forceRoi;
    % start from:
    % -'vesselMaskRefined.nii.gz'
    % do:
    % -'vesselMaskRefinedCalcarine.nii.gz'
    %   include only vessels to be analyzed (nice cross-sections in calcarine)

    % -'calcarineVessel.label'
    %  from vesselMaskRefinedCalcarine.nii.gz
    %  create new volume, select label and identify vein and arteries
    %  and unknown vessels

    % -'sesAvCat_cat_av_preproc_vesselRoi##x.nii.gz'
    %   indivudual vessels from vesselMaskRefinedCalcarine
    %   replace x with v for vein and a for artery, keep x when ambiguous
    %   replace ## with vessel number (v and a pooled in same counter)
    % -'/derivatives/set-avMap/sub-vsmDrivenP2_ses-1_acq-avMap_run-6_echo-cat_angio/manTo-vfMRI.lta'
    %   manually align avMap to functionals (optional)
    if ~contains(volAnat.mask.ulay.f,info.sub); dbstack; error('possibly wrong sub'); end
    % fieldList = fields(runCond);


    % avMap     = char(runCond.(fieldList{1}).avMap.f);
    % avMapReg  = replace(avMap,'in-vfMRI.nii.gz','manTo-vfMRI.lta');
    ulay      = volAnat.mask.ulay.f;
    vMask     = volAnat.mask.vesselRefined.f;
    % vMaskCort = replace(vMask,'.nii.gz','CorticalNoSinus.nii.gz');
    vRoiCalc  = replace(vMask,'.nii.gz','Calcarine.nii.gz');
    vRoiVesselList  = replace(vMask,'vesselMaskRefined.nii.gz','vesselRoi*.nii.gz'); vRoiVesselList = dir(vRoiVesselList); if ~isempty(vRoiVesselList); vRoiVesselList = fullfile({vRoiVesselList.folder},{vRoiVesselList.name})'; else; vRoiVesselList = {}; end
    % if ~exist(vMaskCort,'file')
    %     copyfile(vMask,vMaskCort);
    % end
    % if ~exist(vRoiCalc,'file')
    %     copyfile(vMaskCort,vRoiCalc);
    % end
    if ~exist(vRoiCalc,'file')
        copyfile(vMask,vRoiCalc);
    end


    cmd = {srcFs};
    cmd{end+1,1} = 'freeview \';
    cmd{end+1} = [ulay     ' \'];
    cmd{end+1} = [vMask    ' \'];
    % if exist(vMaskCort,'file')
    %     cmd{end+1} = [vMaskCort ' \'];
    % end
    if exist(vRoiCalc,'file')
        cmd{end+1} = [vRoiCalc ' \'];
    end
    if ~isempty(vRoiVesselList)
        cmd{end+1} = [strjoin(vRoiVesselList,' ') ' \'];
    end

    if isfield(volAnat.label,'calcarineVessel') && isfield(volAnat.label.calcarineVessel,'f') && ~isempty(volAnat.label.calcarineVessel.f)...
            && exist(volAnat.label.calcarineVessel.f,'file')
        cmd{end+1} = [volAnat.label.calcarineVessel.f ':colormap=lut \'];
    end

    for i = 1:length(avMap)
        dbstack; error('double-check that')
        if isMRI(avMap(i))
            f = avMap(i).fspec;
        else
            f = avMap(i).fList;
        end
        manReg = replace(f,'.nii.gz','.manReg.lta');

        if exist(char(manReg),'file')
            avMap(i).manReg = manReg;
            cmd{end+1} = [char(f) ':reg=' char(manReg) ':resample=cubic \'];
        else
            avMap(i).manReg = [];
            cmd{end+1} = [char(f)                      ':resample=cubic \'];
        end
    end
    cmd{end}(end) = [];
    cmd{end} = [cmd{end} '&' newline];
    clipboard('copy',strjoin(cmd,newline))
    disp(strjoin(cmd,newline))
    disp('Freesurfer command in clipboard')
    for i = 1:length(avMap)
        if isempty(avMap(i).manReg) && ~exist(char(avMap(i).manReg),'file')
            if isMRI(avMap(i))
                f = avMap(i).fspec;
            else
                f = avMap(i).fList;
            end
            [a,b,c] = fileparts(replace(char(f),'.nii.gz',''));
            disp(['Save manual registration of avMap:' newline ' ' b '.nii.gz' newline 'As:' newline ' ' char(avMap(i).manReg)])
        end
    end
    if ~exist(volAnat.label.calcarineVessel.f,'file')
        disp(['save calcarine vessel labels as:' newline char(volAnat.label.calcarineVessel.f)])
    end
    if forceThis
        dbstack;
        disp('Go into the script and draw vessel ROIs')
        keyboard
    end
    disp('FS color LUT');
    disp('902->artery');
    disp('914->vein');
    disp('30 ->Left-vessel');
    disp('62 ->Right-vessel');

    volAnat.label.cmd = cmd;


    keyboard

    %% Save to bids derivatives for posterity
    fIn = volAnat.label.calcarineVessel.f;
    [a,b,c] = fileparts(replace(fIn,'.nii.gz',''));
    fOut = fullfile(info.bidsDir,'derivatives','manual',['sub-' info.sub],[b '.nii.gz']);
    copyfile(fIn,fOut);
    volAnat.label.calcarineVessel.f = fOut;



    % %% Vessel + tissue ROI (tissue surrounding each vessels)
    % skipThis = 1;
    % fIn  = volAnat.label.calcarineVessel.f;
    % fOut = fullfile(fileparts(fIn),'calcarineVesselPlusTissueSurround.nii.gz');
    % volAnat.label.calcarineVesselPlusTissueSurround.f = fOut;
    % if ~skipThis && ( force || ~exist(fOut,'file') )
    %     dbstack; error('double-check that')
    %     mri = MRIread(fIn);
    %     mask = mri.vol>0;
    %     sz = size(mask);
    %
    %     % % % figure('WindowStyle','docked'); Lim = [129.0535  270.9403  139.3679  281.2547];
    %     % % % imagesc(mask); hold on; ax = gca; ax.DataAspectRatio = [1 1 1]; axis(Lim)
    %     % % % % get index of non-zero voxels
    %     % % % [y,x] = ind2sub(sz,find(mask));
    %     % % % % get index of these index that form the convex hull (smallest convex region containing all non-zero voxels)
    %     % % % k = convhull(y,x);
    %     % % % % convert to mask
    %     % % % mask = poly2mask(x(k),y(k),sz(1),sz(2));
    %     % % % % % % imagesc(mask);
    %     % % % % dilate mask n voxels
    %     Nmm = 2;
    %     n = round(Nmm/mean(mri.volres(1:2)));
    %     Nmm = n.*mean(mri.volres(1:2));
    %     disp(['dilating vessel mask by ' num2str(Nmm) 'mm around each vessels'])
    %     disp(fOut)
    %     mask = imdilate(mask,strel('disk',n));
    %     % % % figure('WindowStyle','docked'); Lim = [129.0535  270.9403  139.3679  281.2547];
    %     % % % imagesc(mask);
    %
    %     % write to file
    %     mri.vol = mask;
    %     MRIwrite(mri,fOut)
    % end
    %
    %
    % %% Vessel + tissue ROI (all tissue between vessels)
    % skipThis = 1;
    % fIn  = volAnat.label.calcarineVessel.f;
    % fOut = fullfile(fileparts(fIn),'calcarineVesselPlusTissueFilling.nii.gz');
    % volAnat.label.calcarineVesselPlusTissueFilling.f = fOut;
    % if ~skipThis && ( force || ~exist(fOut,'file') )
    %     dbstack; error('double-check that')
    %     mri = MRIread(fIn);
    %     mask = mri.vol>0;
    %     sz = size(mask);
    %
    %     % % % figure('WindowStyle','docked'); Lim = [129.0535  270.9403  139.3679  281.2547];
    %     % % % imagesc(mask); hold on; ax = gca; ax.DataAspectRatio = [1 1 1]; axis(Lim)
    %     % get index of non-zero voxels
    %     [y,x] = ind2sub(sz,find(mask));
    %     % get index of these index that form the convex hull (smallest convex region containing all non-zero voxels)
    %     k = convhull(y,x);
    %     % convert to mask
    %     mask = poly2mask(x(k),y(k),sz(1),sz(2));
    %     % % % imagesc(mask);
    %     % dilate mask n voxels
    %     Nmm = 2;
    %     n = round(Nmm/mean(mri.volres(1:2)));
    %     Nmm = n.*mean(mri.volres(1:2));
    %     disp(['dilating mask of vessel and in-between tissue by ' num2str(Nmm) 'mm'])
    %     disp(fOut)
    %     mask = imdilate(mask,strel('disk',n));
    %     % % % imagesc(mask);
    %
    %     % write to file
    %     mri.vol = mask;
    %     MRIwrite(mri,fOut)
    % end