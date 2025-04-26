function [volAnat,rCond] = volAnatPreproc6(rCond,force,verbose)
    global src
    if ~exist('force','var');     force = []; end
    if ~exist('verbose','var'); verbose = []; end
    if isempty(force);            force = 0; end
    if isempty(verbose);        verbose = 0; end

        

    %% Summarize across runs
    %%% Combine tasks
    taskList = fields(rCond); taskList(~contains(taskList,'task_')) = [];
    if length(taskList)>1; dbstack; error('need to combine data from multiple files'); end
    
    disp('--------------------------------')
    disp(['volAnat: ' rCond.(taskList{1}).sub '_acq-' rCond.(taskList{1}).acq '_prsc-' rCond.(taskList{1}).prsc '_venc-' rCond.(taskList{1}).vencAcq])
    disp('--------------------------------')


    %%% Combine runs
    fList = cell(size(taskList));
    for T = 1:length(taskList)
        fList{T} = rCond.(char(taskList{T})).fPreprocList;
    end
    fList = cat(1,fList{:});
    
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
    fMaskBrainInv = unique(rCond.(char(taskList)).fPreprocMaskList); fMaskBrainInv(cellfun('isempty',fMaskBrainInv)) = [];
    if length(fMaskBrainInv)>1; dbstack; error('more than one mask found'); end; fMaskBrainInv = char(fMaskBrainInv);
    %%%% conform
    fMaskBrainInv2 = replace(fAvCatAv,'_volTs.nii.gz','_brainMaskInv.nii.gz'); copyfile(fMaskBrainInv,fMaskBrainInv2);
    fMaskBrainInv = fMaskBrainInv2;
    MRIconform(fMaskBrainInv,fAvCatAv);
    fMaskBrain = replace(fMaskBrainInv,'_brainMaskInv.nii.gz','_brainMask.nii.gz');
    if force || ~exist(fMaskBrain,'file')
        cmd = {src.afni};
        cmd{end+1} = ['3dcalc -overwrite -a ' fMaskBrainInv ' -expr ''-(a-1)'' -prefix ' fMaskBrain];
        if verbose
            [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
        else
            [status,cmdout] = system(strjoin(cmd,newline)        ); if status; dbstack; error(cmdout); error('x'); end
        end
    end
    MRIconform(fMaskBrain,fAvCatAv);



    %%% Output file index
    acqLabel = strjoin({['acq-' rCond.(char(taskList)).acq] ['prsc-' rCond.(char(taskList)).prsc]},'_');
    volAnat.mask.brain.f     = fMaskBrain;
    volAnat.mask.brain.fInv  = fMaskBrainInv;
    volAnat.mask.brain.fBase = fAvCatAv;
    

    % task = char(taskList);
    % acq = rCond.(task).acq;
    % prsc = rCond.(task).prsc;
    

    %% Preprocess anat for vesselness map (for a starting point to vessel drawing)
    forceThis = force;
    %%%% correct bias field
    [fVolCorr,fVolTsCorr,fCatAvCorr,fVol,fVolField] = correctBiasField(fAvCatAv, fMaskBrain, fCatAv, [], forceThis, verbose);

    forceThis = force;
    %%%% compute vesselness
    [fComp,fNonComp,fSegFig] = computeVesselness(fVolCorr,fMaskBrain,forceThis,verbose);

    
    %% Draw vessel rois
    derivDir = fullfile(rCond.(char(taskList)).dirsOrig.bidsDeriv,acqLabel);
    fVesselRoi  = fullfile(derivDir,'vessel.nii.gz');

    if force || ~exist(fVesselRoi,'file')

        % rCond.(taskList{1})
        
        % fAvMap
        % fTof
        % drawVesselRoi([cellstr(fVolCorr) cellstr(fCatAvCorr)],fVesselRoi,cellstr{fComp{contains(fComp,'vesselSegMask.nii.gz')}},fMaskBrain,fAvMap,fTof)
        
        disp('!!!!!!!!!!')
        disp('!!!!!!!!!!')
        dbstack; warning('not implemented yet')
        disp('!!!!!!!!!!')
        disp('!!!!!!!!!!')
    end


    %%% Output file index
    volAnat.label.calcarineVessel.f        = fVesselRoi;
    volAnat.label.calcarineVessel.label    = {'artery' 'vein' 'Left-vessel' 'Right-vessel'};
    volAnat.label.calcarineVessel.labelVal = [902 914 30 62];
    volAnat.label.calcarineVessel.fBase = fVolCorr;
    volAnat.label.calcarineVessel.fFig  = fSegFig;
    if verbose
        disp('FS color LUT');
        disp('902->artery');
        disp('914->vein');
        disp('30 ->Left-vessel');
        disp('62 ->Right-vessel');
    end





    %% Individual vessel ROIs
    if exist(fVesselRoi,'file')
        label = volAnat.label.calcarineVessel;
        imField = {'base'};
        im = {label.fBase};
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
        rCond.(char(taskList)).volAnat = volAnat;
    end





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