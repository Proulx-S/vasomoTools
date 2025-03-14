function [out,avMap] = volAnatPreproc5(do,info,runCond,prcSmr,acq,avMap,force,forceRoi,verbose)
global srcFs srcAfni
if ~exist('avMap','var');     avMap = []; end
if ~exist('force','var');     force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force);            force = 0; end
if isempty(verbose);        verbose = 0; end

%% User variables
outVar = 'volAnat';
stepLabel = 'anat preprocessing';





%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if isfield(info,'outDir'); outDir = info.outDir; else, outDir = info.prcDir; end; if ~exist(outDir,'dir'); mkdir(outDir); end
if isfield(info,'ses')
    stepFile = fullfile(outDir,[strjoin({['sub-' info.sub] ['ses-' info.ses] mfilename},'_')]);
else
    stepFile = fullfile(outDir,[strjoin({['sub-' info.sub] mfilename},'_')]);
end

if exist('do','var') && ~isempty(do)
    if isfield(do,'loadIt') && ~isempty(do.loadIt); loadIt = do.loadIt; else, loadIt=0; end
    if isfield(do,'doIt')   && ~isempty(do.doIt);   doIt   = do.doIt;   else, doIt=0;   end
    if isfield(do,'saveIt') && ~isempty(do.saveIt); saveIt = do.saveIt; else, saveIt=0; end
else
    loadIt = 0; doIt = 1; saveIt = 0;
end

if loadIt && doIt
    warning(strjoin({'does not make sense to loadIt then doIt' 'Changing to' 'loadIt = 0' 'doIt   = 1' 'saveIt = 1'},newline)); loadIt = 0; doIt = 1; saveIt = 1;
end
if loadIt && saveIt
    warning(strjoin({'does not make sense to loadIt then saveIt' 'Changing to' 'loadIt = 1' 'doIt   = 0' 'saveIt = 0'},newline)); loadIt = 1; doIt = 0; saveIt = 0;
end

tic; disp(' '); disp(' '); disp(repmat('-',1,length(stepLabel))); disp(upper(stepLabel)); disp(repmat('-',1,length(stepLabel)))

if loadIt
    if exist([stepFile '.mat'],'file')
        disp(strjoin({[upper(stepLabel) ': loading from'] [stepFile '.mat']},newline)); load(stepFile); disp([upper(stepLabel) ': loaded'])
    else
        warning(strjoin({[stepFile '.mat'] 'does not exist' 'Changing to' 'loadIt = 0' 'doIt   = 1' 'saveIt = 1'},newline)); loadIt = 0; doIt = 1; saveIt = 1;
    end
end






if doIt

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% Do the processing here
    %%%%%%%%%%%%%%%%%%%%%%%%%






    if ~isfield(info,'useSynth');     info.useSynth = []; end
    if ~isfield(info,'skipHead');     info.skipHead = []; end
    if ~isfield(info,'skipVessel'); info.skipVessel = []; end
    if isempty(info.useSynth);        info.useSynth = 0; end
    if isempty(info.skipHead);        info.skipHead = 0; end
    if isempty(info.skipVessel);    info.skipVessel = 0; end



    %% %%%%%%
    % Masks %
    %%%%%% %%

    %% Image or images to use for creating mask
    % try
    if numel(prcSmr.sesCat.runAv.fList)>1; dbstack; error('X'); end
    procDir = fullfile(outDir,['sub-' info.sub],'ses-cat',['set-' acq]);
    if ~exist(procDir,'dir'); mkdir(procDir); end
    volAnat.mask.ulay.f = fullfile(procDir,'prcMeans.nii.gz');
    copyfile(char(prcSmr.sesCat.runAv.fList),volAnat.mask.ulay.f)


    % volAnat.mask.ulay.f = char(prcSmr.sesCat.runAv.fList);
    % volAnat.mask.ulay.f = fields(runCond);
    % volAnat.mask.ulay.f = runCond.(volAnat.mask.ulay.f{1}).fPreprocUnderSesCatRunCatAvList{1};


    %% Simple crop
    volAnat.mask.crop.mri = MRIread(volAnat.mask.ulay.f,1);
    %%% Aim to crop out 'crop'mm in each direction
    crop = 5; %mm
    % try
    trim = ceil(crop./volAnat.mask.crop.mri.volres); %voxel
    % catch
    %     keyboard
    % end
    trimActual = trim.*volAnat.mask.crop.mri.volres; %mm
    volAnat.mask.crop.mri.vol = true(volAnat.mask.crop.mri.volsize);
    volAnat.mask.crop.mri.vol([1:trim(1) end-trim(1)+1:end],:                            ,:) = false;
    volAnat.mask.crop.mri.vol(:                            ,[1:trim(2) end-trim(2)+1:end],:) = false;

    trimZ = trim(3);
    while volAnat.mask.crop.mri.volsize(3) - trimZ*2 < 5
        trimZ = trimZ-1;
    end
    volAnat.mask.crop.mri.vol(:,:,[1:trimZ end-trimZ+1:end]) = false;




    % %% Brain mask
    % if volAnat.mask.crop.mri.depth>5
    %     info.useSynth = 1;
    % end
    %
    % cmd = {srcFs}; cmd{end+1} = srcAfni;
    % if info.useSynth
    %     cmd{end+1} = 'echo mri_synthstrip';
    %     cmd{end+1} = 'mri_synthstrip \';
    %     cmd{end+1} = ['-i ' volAnat.mask.ulay.f ' \'];
    %     fOut = replace(volAnat.mask.ulay.f,'_volTs.nii.gz','_volSynthMask.nii.gz');
    %     cmd{end+1} = ['-m ' fOut];
    %
    %     fIn = fOut;
    %     fMask = replace(fIn,'_volSynthMask.nii.gz','_volBrainMask.nii.gz');
    %     fMaskInv = replace(fIn,'_volSynthMask.nii.gz','_volBrainMaskInv.nii.gz');
    %     cmd{end+1} = ['cp ' fIn ' ' fMask];
    %
    %     cmd{end+1} = '3dmask_tool -overwrite \';
    %     cmd{end+1} = ['-input ' fMask ' \'];
    %     cmd{end+1} = ['-prefix ' fMask ' \'];
    %     cmd{end+1} = ['-dilate_input 2'];
    %
    %     cmd{end+1} = 'echo edit brain mask if needed';
    %     cmd{end+1} = ['fslview ' volAnat.mask.ulay.f ' ' fMask];
    %
    %     if force || ~exist(fMask,'file') || ~exist(fMaskInv,'file')
    %         % disp('edit brain mask if needed')
    %         [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    %
    %         % zero first and last slice
    %         mask = MRIread(fMask);
    %         mask.vol(:,:,[1 end]) = 0;
    %         MRIwrite(mask,fMask);
    %
    %         % invert mask
    %         mask.vol = -(mask.vol-1);
    %         MRIwrite(mask,fMaskInv);
    %     end
    %
    % else
    %
    %     cmd{end+1} = ['fslview -m single ' volAnat.mask.ulay.f];
    %
    %     fIn = replace(volAnat.mask.ulay.f,'_volTs.nii.gz','_volTs-mask.nii.gz');
    %     fMask = replace(volAnat.mask.ulay.f,'_volTs.nii.gz','_volBrainMask.nii.gz');
    %     fMaskInv = replace(fMask,'_volBrainMask.nii.gz','_volBrainMaskInv.nii.gz');
    %
    %     if force || ~exist(fMask,'file') || ~exist(fMaskInv,'file')
    %         disp('draw brain mask')
    %         [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    %         movefile(fIn,fMask)
    %
    %         cmd = {srcAfni};
    %         cmd{end+1} = '3dcalc -overwrite \';
    %         cmd{end+1} = ['-prefix ' fMaskInv ' \'];
    %         cmd{end+1} = ['-a ' fMask ' \'];
    %         cmd{end+1} = '-expr ''-(a-1)''';
    %         [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    %     end
    % end
    %
    % volAnat.mask.brain.f = fMask;
    % volAnat.mask.brain.mri = MRIread(fMask);
    % volAnat.mask.brainInv.f = fMaskInv;
    % volAnat.mask.brainInv.mri = MRIread(fMaskInv);
    %
    %% Head mask
    if ~info.skipHead
        cmd = {};
        % cmd = {srcFs};
        cmd{end+1} = ['fslview -m single ' volAnat.mask.ulay.f];
        cmd = strjoin(cmd,newline); % disp(cmd)

        if ~contains(volAnat.mask.ulay.f,'.nii.gz'); dbstack; error('unexpected filename'); end
        fIn = replace(volAnat.mask.ulay.f,'.nii.gz','-mask.nii.gz');
        % fOut = replace(volAnat.mask.ulay.f,'_volTs.nii.gz','_volHeadMask.nii.gz');
        fOut = fullfile(info.bidsDir,'derivatives','manual',['sub-' info.sub],'volHeadMask.nii.gz');

        if force || ~exist(fOut,'file')
            disp('draw head mask')
            [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
            movefile(fIn,fOut)
        end

        volAnat.mask.head.f   = fOut;
        % volAnat.mask.head.mri = MRIread(fOut);
    end




    %% %%%%%
    % ROIs %
    %%%%% %%

    if ~info.skipVessel


        % -'sesAvCat_cat_av_preproc_vesselMaskRefinedCorticalNoSinus.nii.gz'
        %   only cortical vessels and exluding sinuses
        % -'sesAvCat_cat_av_preproc_vesselMaskRefinedCalcarine.nii.gz'
        %   include only vessels to be analyzed (nice cross-sections in calcarine)

        % -'calcarineVessel.label'
        %  from sesAvCat_cat_av_preproc_vesselMaskRefinedCalcarine.nii.gz
        %  create new volume, select label and identify vein and arteries
        %  and unknown vessels

        % -'sesAvCat_cat_av_preproc_vesselRoi##x.nii.gz'
        %   indivudual vessels from vesselMaskRefinedCalcarine
        %   replace x with v for vein and a for artery, keep x when ambiguous
        %   replace ## with vessel number (v and a pooled in same counter)

        %%% Mask
        volAnat.mask.vessel.f = dir(fullfile(info.bidsDir,'derivatives','manual',['sub-' info.sub],'calcarineVessel*.nii.gz'));
        volAnat.mask.vessel.f = fullfile({volAnat.mask.vessel.f.folder},{volAnat.mask.vessel.f.name})';
        if ~isempty(volAnat.mask.vessel.f) %&& exist(char(volAnat.mask.vessel.f),'file')
            %%% Label
            [~,i,~] = fileparts(replace(volAnat.mask.vessel.f,'.nii.gz',''));
            i = ismember(i,'calcarineVessel');
            if nnz(i)~=1; dbstack; error('somehting wrong'); end
            volAnat.label.vessel.f = volAnat.mask.vessel.f(i); volAnat.mask.vessel.f(i) = [];
            volAnat.label.vessel.label = {'artery' 'vein' 'Left-vessel' 'Right-vessel'};
            volAnat.label.vessel.val   = [ 902      914    30            62           ];
            %%% Roi
            volAnat.roi = [];



        else
            %% Image or images to use for creating ROIs

            mri = MRIload3(volAnat.mask.ulay.f,[],[],0);
            mask = MRIload3(volAnat.mask.head.f,[],[],0); mask = mask.vol;
            mask = mask & volAnat.mask.crop.mri.vol;
            mri = vol2vec(mri,mask);

            %% Liberal Vessel Mask
            volAnat.mask.vessel.f = fullfile(procDir,'vesselMask.nii.gz');
            volAnat.mask.vessel.fFig = replace(volAnat.mask.vessel.f,'.nii.gz','.fig');

            %%% Identify vessels based on gaussian mixture
            if force || ~exist(volAnat.mask.vessel.f,'file') || ~exist(volAnat.mask.vessel.fFig,'file')
                % gaussian mixture
                k = 2;
                X = mri.vec(:);

                hFig = figure('visible','off');
                hT = tiledlayout(1,2); hT.TileSpacing = "tight"; hT.Padding = 'tight';

                nexttile
                h = histogram(X,'Normalization','pdf'); hold on
                binCent = h.BinEdges(1:end-1) - mean(diff(h.BinEdges))/2;
                GMModel = fitgmdist(X,k);
                plot(binCent',GMModel.pdf(binCent'))
                for i = 1:k
                    n = makedist('normal',GMModel.mu(i),sqrt(GMModel.Sigma(i)));
                    p = GMModel.ComponentProportion(i);
                    plot(binCent,p.*pdf(n,binCent))
                end
                labelVal = [0 1:k];
                label = {'outerHead' '' ''};
                [~,vesselComp] = min(GMModel.ComponentProportion);
                label{labelVal==vesselComp} = 'vessel';
                label{cellfun('isempty',label)} = 'nonVessel';

                legend([{'hist' 'pdf'} label(2:end)])

                nexttile
                [idx,nlogL,P,logpdf,d2] = cluster(GMModel,X);
                seg = zeros(size(mri.vec));
                seg(:) = idx;
                mriVessel = mri;
                mriVessel.vec = seg==labelVal(ismember(label,'vessel'));
                mriVessel = vec2vol(mriVessel);
                mriVessel.vol = any(mriVessel.vol,4);
                mriVessel.nframes = 1;
                mriVessel.fspec = volAnat.mask.vessel.f;
                imagesc(mriVessel.vol)
                ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';

                MRIwrite(mriVessel,mriVessel.fspec);
                volAnat.mask.vessel.mri = mriVessel;

                fFig = volAnat.mask.vessel.fFig;
                if verbose<=1
                    set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
                end
                savefig(hFig,fFig,'compact')
                if verbose>1
                    hFig.Visible = 'on';
                    hFig.WindowStyle = 'docked';
                    drawnow
                end
            end

            %%% Manual refinement
            volAnat.mask.vesselRefined.f = replace(volAnat.mask.vessel.f,'vesselMask.nii.gz','vesselMaskRefined.nii.gz');
            if force || ~exist(volAnat.mask.vesselRefined.f,'file')
                mri = volAnat.mask.vessel.mri;
                mri.fspec = volAnat.mask.vesselRefined.f;
                MRIwrite(mri,mri.fspec);

                cmd = {srcFs};
                cmd{end+1} = 'echo "refine vessel mask (just remove false-positives)"';
                cmd{end+1} = ['fslview -m single ' volAnat.mask.ulay.f ' ' volAnat.mask.vesselRefined.f];
                [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
            end


            %%% Dilation
            f = volAnat.mask.vesselRefined.f;
            fDil = replace(f,'vesselMaskRefined.nii.gz','vesselMaskRefinedDil.nii.gz');
            volAnat.mask.vesselRefinedDil.f = fDil;
            if force || ~exist(fDil,'file')
                cmd = {srcAfni};
                cmd{end+1} = '3dmask_tool -overwrite \';
                cmd{end+1} = ['-input ' f ' \'];
                cmd{end+1} = ['-prefix ' fDil ' \'];
                cmd{end+1} = ['-dilate_input 3'];
                [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
            end
            % disp(strjoin({volAnat.mask.ulay.f
            % volAnat.mask.vesselRefined.f
            % volAnat.mask.vesselRefinedDil.f},[' \\' newline]))


            % %% Sort avMap
            % doThis = 0;
            % if isMRI(avMap)
            %     doThis = 0;
            % else
            %     for i = 1:numel(avMap)
            %         doThis = length(avMap(i).fList)>1; if doThis; break; end
            %     end
            % end
            % %%% seperate different maps
            % % each map should be a single XxYxZxE file
            % if doThis
            %     tmp.sub = char(unique({avMap.sub}));
            %     tmp.ses = cat(1,avMap.ses);
            %     if ~isfield(avMap,'label') || isempty(avMap.label) && isfield(avMap,'acq') && ~isempty(avMap.acq)
            %         tmp.acq = char(unique({avMap.acq}));
            %     else
            %         tmp.acq = char(unique({avMap.label}));
            %     end
            %     tmp.fList = {avMap.fspec}';
            %     tmp.mri = avMap;
            %     avMap = tmp; clear tmp
            %     [~,avMap.bidsList,~] = fileparts(replace(avMap.fList,'.nii.gz',''));
            %     for f = 1:length(avMap.bidsList)
            %         avMap.bidsList{f} = strsplit(avMap.bidsList{f},'_');
            %     end
            %     avMap.bidsList = cat(1,avMap.bidsList{:});
            %     un = zeros(1,size(avMap.bidsList,2));
            %     for b = 1:size(avMap.bidsList,2)
            %         un(b) = length(unique(avMap.bidsList(:,b)));
            %     end
            %     eInd = contains(avMap.bidsList(1,:),'echo-');
            %
            %     tmp = avMap.bidsList(:,~eInd);
            %     if any(un(~eInd)>1)
            %         [a,~,c] = unique(tmp(:,un(~eInd)>1));
            %         tmp = repmat(avMap,length(a),1);
            %         tmpList = {'ses' 'fList' 'mri' 'bidsList'};
            %         for i = 1:length(a)
            %             for ii = 1:length(tmpList)
            %                 tmp(i).(tmpList{ii}) = avMap.(tmpList{ii})(c==i,:);
            %             end
            %         end
            %     else
            %         tmp = avMap;
            %     end
            %     %detect un-handled cases
            %     if any([avMap.mri.nframes]>1) || length(unique(avMap.ses))>1
            %         dbstack;
            %         warning('code that')
            %         keyboard
            %     end
            %     avMap = tmp; clear tmp
            %
            %
            %     % all avMap file are the same except echo time, so catenate
            %     % that in a single file
            %     avMap = MRIload2(avMap);
            %     for i = 1:length(avMap)
            %         [~,b] = sort(avMap(i).bidsList(:,eInd));
            %         avMap(i).mri(1).vol = cat(4,avMap(i).mri(b).vol);
            %         avMap(i).mri(2:end) = [];
            %         avMap(i).ses        = avMap(i).ses(b(1),:);
            %         avMap(i).fList      = avMap(i).fList(b(1),:);
            %         avMap(i).bidsList   = avMap(i).bidsList(b(1),:);
            %
            %         avMap(i).fList = replace(avMap(i).fList,'echo-1','echo-cat');
            %         avMap(i).fList = strsplit(char(avMap(i).fList),filesep);
            %         avMap(i).fList = {strjoin([avMap(i).fList(1:end-1) {'derivatives'} avMap(i).fList(end)],filesep)};
            %         avMap(i).bidsList{eInd} = 'echo-cat';
            %         avMap(i).mri.fspec = char(avMap(i).fList);
            %         if ~exist(fileparts(avMap(i).mri.fspec),'dir'); mkdir(fileparts(avMap(i).mri.fspec)); end
            %         MRIwrite(avMap(i).mri,avMap(i).mri.fspec);
            %     end
            % end



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
        end
    end
end







%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

disp(repmat('-',1,length(stepLabel)+6)); disp([upper(stepLabel) ': DONE']); toc; disp(repmat('-',1,length(stepLabel)+6)); disp(' '); disp(' ');
eval(['out = ' outVar '; clear ' outVar]);