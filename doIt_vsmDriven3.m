% Branching from doIt_vsmDriven2.m. Here we start doing
% run-condition-specific functional analysis
clear all
close all
[outDir,pipId] = fileparts(mfilename('fullpath'));
outDir = fullfile(outDir,pipId); if ~exist(outDir,'dir'); mkdir(outDir); end



%%%%%%%%%%%%%%%%
%% Preparation %
%%%%%%%%%%%%%%%%

%%% Dependencies
% matlab
addpath(genpath(fullfile(pwd,pipId)))
addpath(genpath('/usr/local/freesurfer/stable7.4.1/matlab/'))
addpath(genpath('/space/takoyaki/1/users/proulxs/tools/chronux'))
addpath(genpath('/space/takoyaki/1/users/proulxs/tools/vasomoTools'))
addpath(genpath('/space/takoyaki/1/users/proulxs/tools/bassReg2'))
addpath(genpath('/space/takoyaki/1/users/proulxs/tools/martinosTools'))
% bash
global srcAfni srcFs
srcFs = 'source /usr/local/freesurfer/fs-stable741-env-autoselect';
srcAfni = 'export PATH=$PATH:/usr/pubsw/packages/AFNI/23.1.05';

info.dataSetLabel = 'vsmDriven'; % 'svfmri_long_short' 'satin' 'mitraShahin' 'feruJingyuan'


%%% Variables and Paths
switch info.dataSetLabel
    case 'vsmDriven'
        info.preprocDir = outDir;
        info.pipId      = pipId;
        info.bidsDir    = fullfile(info.preprocDir,'bids'); if ~exist(info.bidsDir,'dir'); mkdir(info.bidsDir); end
        %%% Data location
        info.dbDir      = '/space/takoyaki/1/users/proulxs/vasomo/source/expDb';

    otherwise
        dbstack; error('code that')
end
%% %%%%%%%%%%%%%


%% %%%%%%%%%%%%%%
% Preprocessing %
%%%%%%%%%%%%%% %%
if 0


    subList       = {};
    sesList       = {};
    bidsDirList   = {};
    sesDbList     = {};
    sesDbFileList = {};
    volTrList     = {};
    runCondList   = {};
    volTrList     = {};
    dummyList     = {};
    kList         = {};
    winLengthList = {};
    winStepList   = {};
    onsetsList    = {};
    ondursList    = {};
    dtStimList    = {};
    labelList     = {};
    extraList     = {};


    %%% Subject and session
    switch info.dataSetLabel
        case 'vsmDriven'
            K         = 5;
            winLength = 15;


            runCond  = {};
            fMap     = {};
            avMap    = {};
            memprage = {};
            sa2rage  = {};
            pcMRA    = {};
            b0       = {};
            b1       = {};


            %%%%%%%%%%%%%%
            %%%% vsmDriven
            sesDbListTmp{1} = {
                '/autofs/space/takoyaki_001/users/proulxs/vasomo/source/expDb/vsmDrivenP1/2024-07-28--bay2--vsmDrivenP1'
                '/autofs/space/takoyaki_001/users/proulxs/vasomo/source/expDb/vsmDrivenP1/2024-08-09--bay2--vsmDrivenP1'
                };
            sesDbListTmp{2} = {
                '/autofs/space/takoyaki_001/users/proulxs/vasomo/source/expDb/vsmDrivenP2/2024-07-28--bay2--vsmDrivenP2'
                '/autofs/space/takoyaki_001/users/proulxs/vasomo/source/expDb/vsmDrivenP2/2024-08-05--bay2--vsmDrivenP2'
                };
            


            for sub = 1:2
                for ses = 1:2
                    subList{end+1,1}     = ['vsmDrivenP' num2str(sub)];
                    sesList{end+1,1}     = num2str(ses);
                    bidsDirList{end+1,1} = fullfile(info.bidsDir,['sub-' subList{end}],['ses-' sesList{end}]);
                    sesDbList{end+1,1}   = sesDbListTmp{sub}{ses};
                    extraList{end+1,1}   = '';


                    %%%%% anat

                    %%%%%% avMap
                    avMap{end+1,1}.fList = dir(fullfile(bidsDirList{end,1},'anat','*acq-avMap*.nii.gz'));
                    %%%%%% memprage
                    memprage{end+1,1}.fList = dir(fullfile(bidsDirList{end,1},'anat','*_T1w.nii.gz'));
                    if ~isempty(memprage{end,1}.fList)
                        tmp = dir(fullfile(info.dbDir,subList{end},'*','fs',subList{end}));
                        memprage{end,1}.fsDir = tmp(1).folder; clear tmp
                    end
                    %%%%%% pcMRA
                    pcMRA{end+1,1}.fList = dir(fullfile(bidsDirList{end,1},'anat','*acq-venc*.nii.gz'));


                    %%%%% fmap
                    %%%%%% topup
                    b0{end+1,1}.label = 'topup';
                    b0{end,1}.fList = dir(fullfile(bidsDirList{end,1},'fmap','*_epi.nii.gz'));
                    %%%%%% sa2rage
                    b1{end+1,1}.label = 'B1';
                    b1{end,1}.fList = dir(fullfile(bidsDirList{end,1},'fmap','*_TB1SRGE.nii.gz'));


                    %%%%% func
                    runCond{end+1,1}   = {};
                    dummyList{end+1,1} = {};

                    %%%%%% vfMRI
                    dummy            = 5;
                    %%%%%%% 50sPrd5sDur
                    runCond{end,1}{1,end+1}.sub  = subList{end};
                    runCond{end,1}{1,end}.ses    = sesList{end};
                    runCond{end,1}{1,end}.label  = '50sPrd5sDur';
                    runCond{end,1}{1,end}.labelAcq = 'vfMRI';
                    runCond{end,1}{1,end}.fList = dir(fullfile(bidsDirList{end,1},'func',['*task-' runCond{end,1}{1,end}.label '*_angio.nii.gz']));
                    runCond{end,1}{1,end}.fList = fullfile({runCond{end,1}{1,end}.fList.folder},{runCond{end,1}{1,end}.fList.name})';
                    dsgn.label = runCond{end,1}{1,end}.label;
                    dsgn.dt    = 0.840;
                    initRest   = dsgn.dt*12;
                    stimPeriod = dsgn.dt*57;
                    stimDur    = dsgn.dt*6;
                    runDur     = dsgn.dt*354;
                    dsgn.onsetList = initRest:stimPeriod:(runDur-stimPeriod);
                    dsgn.ondurList = ones(size(dsgn.onsetList)).*(stimDur);
                    runCond{end,1}{1,end}.dsgn = dsgn;
                    dummyList{end,1}{1,end+1}  = repmat(dummy,size(runCond{end,1}{1,end}.fList));

                    %%%%%% 50sPrd1sDur
                    runCond{end,1}{1,end+1}.sub  = subList{end};
                    runCond{end,1}{1,end}.ses    = sesList{end};
                    runCond{end,1}{1,end}.label  = '50sPrd1sDur';
                    runCond{end,1}{1,end}.labelAcq = 'vfMRI';
                    runCond{end,1}{1,end}.fList = dir(fullfile(bidsDirList{end,1},'func',['*task-' runCond{end,1}{1,end}.label '*_angio.nii.gz']));
                    runCond{end,1}{1,end}.fList = fullfile({runCond{end,1}{1,end}.fList.folder},{runCond{end,1}{1,end}.fList.name})';
                    dsgn.label = runCond{end,1}{1,end}.label;
                    dsgn.dt    = 0.840;
                    initRest   = dsgn.dt*12;
                    stimPeriod = dsgn.dt*57;
                    stimDur    = dsgn.dt*1;
                    runDur     = dsgn.dt*354;
                    dsgn.onsetList = initRest:stimPeriod:(runDur-stimPeriod);
                    dsgn.ondurList = ones(size(dsgn.onsetList)).*(stimDur);
                    runCond{end,1}{1,end}.dsgn = dsgn;
                    dummyList{end,1}{1,end+1}  = repmat(dummy,size(runCond{end,1}{1,end}.fList));

                    %%%%%% 10sPrd1sDur
                    runCond{end,1}{1,end+1}.sub  = subList{end};
                    runCond{end,1}{1,end}.ses    = sesList{end};
                    runCond{end,1}{1,end}.label  = '10sPrd1sDur';
                    runCond{end,1}{1,end}.labelAcq = 'vfMRI';
                    runCond{end,1}{1,end}.fList = dir(fullfile(bidsDirList{end,1},'func',['*task-' runCond{end,1}{1,end}.label '*_angio.nii.gz']));
                    runCond{end,1}{1,end}.fList = fullfile({runCond{end,1}{1,end}.fList.folder},{runCond{end,1}{1,end}.fList.name})';
                    dsgn.label = runCond{end,1}{1,end}.label;
                    dsgn.dt    = 0.840;
                    initRest   = dsgn.dt*12;
                    stimPeriod = dsgn.dt*12;
                    stimDur    = dsgn.dt*1;
                    runDur     = dsgn.dt*354;
                    dsgn.onsetList = initRest:stimPeriod:(runDur-stimPeriod);
                    dsgn.nullTrial = false(size(dsgn.onsetList));
                    dsgn.nullTrial(14) = true;
                    dsgn.ondurList = ones(size(dsgn.onsetList)).*stimDur;
                    runCond{end,1}{1,end}.dsgn = dsgn;
                    dummyList{end,1}{1,end+1}  = repmat(dummy,size(runCond{end,1}{1,end}.fList));

                    %%%%%% 15sPrd1sDur
                    runCond{end,1}{1,end+1}.sub  = subList{end};
                    runCond{end,1}{1,end}.ses    = sesList{end};
                    runCond{end,1}{1,end}.label  = '15sPrd1sDur';
                    runCond{end,1}{1,end}.labelAcq = 'vfMRI';
                    runCond{end,1}{1,end}.fList = dir(fullfile(bidsDirList{end,1},'func',['*task-' runCond{end,1}{1,end}.label '*_angio.nii.gz']));
                    runCond{end,1}{1,end}.fList = fullfile({runCond{end,1}{1,end}.fList.folder},{runCond{end,1}{1,end}.fList.name})';
                    dsgn.label = runCond{end,1}{1,end}.label;
                    dsgn.dt    = 0.840;
                    initRest   = dsgn.dt*12;
                    stimPeriod = dsgn.dt*17;
                    stimDur    = dsgn.dt*1;
                    runDur     = dsgn.dt*354;
                    dsgn.onsetList = initRest:stimPeriod:(runDur-stimPeriod);
                    dsgn.nullTrial = false(size(dsgn.onsetList));
                    dsgn.nullTrial(14) = true;
                    dsgn.ondurList = ones(size(dsgn.onsetList)).*stimDur;
                    runCond{end,1}{1,end}.dsgn = dsgn;
                    dummyList{end,1}{1,end+1}  = repmat(dummy,size(runCond{end,1}{1,end}.fList));

                    %%%%%% 20sPrd1sDur
                    runCond{end,1}{1,end+1}.sub  = subList{end};
                    runCond{end,1}{1,end}.ses    = sesList{end};
                    runCond{end,1}{1,end}.label  = '20sPrd1sDur';
                    runCond{end,1}{1,end}.labelAcq = 'vfMRI';
                    runCond{end,1}{1,end}.fList = dir(fullfile(bidsDirList{end,1},'func',['*task-' runCond{end,1}{1,end}.label '*_angio.nii.gz']));
                    runCond{end,1}{1,end}.fList = fullfile({runCond{end,1}{1,end}.fList.folder},{runCond{end,1}{1,end}.fList.name})';
                    dsgn.label = runCond{end,1}{1,end}.label;
                    dsgn.dt    = 0.840;
                    initRest   = dsgn.dt*12;
                    stimPeriod = dsgn.dt*24;
                    stimDur    = dsgn.dt*1;
                    runDur     = dsgn.dt*354;
                    dsgn.onsetList = initRest:stimPeriod:(runDur-stimPeriod);
                    dsgn.nullTrial = false(size(dsgn.onsetList));
                    dsgn.nullTrial(7) = true;
                    dsgn.ondurList = ones(size(dsgn.onsetList)).*stimDur;
                    runCond{end,1}{1,end}.dsgn = dsgn;
                    dummyList{end,1}{1,end+1}  = repmat(dummy,size(runCond{end,1}{1,end}.fList));


                    %%%%%% bold
                    dummy = 5;
                    %%%%%%% 50sPrd1sDur
                    runCond{end,1}{1,end+1}.sub  = subList{end};
                    runCond{end,1}{1,end}.ses    = sesList{end};
                    runCond{end,1}{1,end}.label  = '50sPrd1sDur';
                    runCond{end,1}{1,end}.labelAcq = 'bold';
                    runCond{end,1}{1,end}.fList = dir(fullfile(bidsDirList{end,1},'func',['*task-' runCond{end,1}{1,end}.label '*_bold.nii.gz']));
                    runCond{end,1}{1,end}.fList = fullfile({runCond{end,1}{1,end}.fList.folder},{runCond{end,1}{1,end}.fList.name})';
                    dsgn.label = runCond{end,1}{1,end}.label;
                    dsgn.dt    = 0.840;
                    initRest   = dsgn.dt*12;
                    stimPeriod = dsgn.dt*57;
                    stimDur    = dsgn.dt*1;
                    runDur     = dsgn.dt*354;
                    dsgn.onsetList = initRest:stimPeriod:(runDur-stimPeriod);
                    dsgn.ondurList = ones(size(dsgn.onsetList)).*(stimDur);
                    runCond{end,1}{1,end}.dsgn = dsgn;
                    dummyList{end,1}{1,end+1}  = repmat(dummy,size(runCond{end,1}{1,end}.fList));



                end
            end
            % subList
            % runCond{1}{:}
            % runCond{2}{:}
            % pcMRA{2}.fList.name
            % memprage{4}
            % avMap{4}
            % b0{4}.fList.name
            % b1{4}.fList.name


    end



    %%%%%%%%%%%%%%
    % Session loop
    runSet  = cell(size(runCond));
    volAnat = cell(size(runCond));

    sesIndList = 1:length(subList);
    for s = 1:length(subList(sesIndList))
        S = sesIndList(s);

        % see doIt_pilotPipeline08c/stepTemplate.m
        %% Housekeeping
        info.sub          = subList{S};
        info.ses          = sesList{S};
        info.sesDb        = sesDbList{S};
        info.bidsDir      = bidsDirList{S};
        info.bidsDerivDir = fullfile(info.bidsDir,'derivatives'); if ~exist(info.bidsDerivDir,'dir'); mkdir(info.bidsDerivDir); end


        tic; disp(' '); disp(' '); disp(' ');
        tmp = {['sub-' info.sub '_ses-' info.ses]};
        tmp{end+1} = [num2str(s) '/' num2str(length(sesIndList))];
        tmp{end+1} = ['Nses=' num2str(length(subList))];
        tmp = strjoin(tmp,'; ');
        disp(repmat('+',1,length(tmp))); disp(repmat('+',1,length(tmp))); disp(repmat('+',1,length(tmp))); disp(tmp)







        %% Main analysis steps








        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Within-run preprocessing %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        do.loadIt = 0;
        do.doIt = 1;
        do.saveIt = 0;

        switch info.dataSetLabel
            case 'vsmDriven'

                %%% vfMRI
                runSet{S}{1,end+1}.sub      = subList{S};
                runSet{S}{end}.ses          = sesList{S};
                runSet{S}{end}.label        = 'vfMRI';
                runSet{S}{end}.wd           = outDir;
                runSet{S}{end}.bidsDir      = bidsDirList{S};
                runSet{S}{end}.bidsDerivDir = fullfile(bidsDirList{S},'derivatives',['set-' runSet{S}{end}.label]);
                ind = [runCond{S}{:}]; ind = {ind.labelAcq}; ind = ismember(ind,runSet{S}{end}.label);
                runSet{S}{end}.fList = [runCond{S}{ind}];
                runSet{S}{end}.fList  = cat(1,runSet{S}{end}.fList.fList);
                runSet{S}{end}.nDummy = cat(1,dummyList{S}{ind});
                if ~isempty(runSet{S}{end}.fList)
                    forceThis   = 0;
                    verboseThis = 1;
                    runSet{S}{end}.initFiles = initPreproc(runSet{S}{end},[],[],forceThis,verboseThis);

                    forceThis   = 0;
                    verboseThis = 1;
                    info.useSynth = 0;
                    fMask = volAnatPreproc2(do,info,runSet{S}{end},forceThis,verboseThis);
                    fMask = fMask.func.mask.brainInv.mri.fspec;

                    forceThis   = 0;
                    verboseThis = 1;
                    param.baseType = 'first'; % 'first' 'av' 'mcAv'
                    param.spSmFac  = 1;
                    fBase = [];
                    runSet{S}{end}.wrMocoFiles = estimMotionWR2(runSet{S}{end}.initFiles,param,fBase,fMask,forceThis,verboseThis);
                end


                %%% bold
                runSet{S}{1,end+1}.sub      = subList{S};
                runSet{S}{end}.ses          = sesList{S};
                runSet{S}{end}.label        = 'bold';
                runSet{S}{end}.wd           = outDir;
                runSet{S}{end}.bidsDir      = bidsDirList{S};
                runSet{S}{end}.bidsDerivDir = fullfile(bidsDirList{S},'derivatives',['set-' runSet{S}{end}.label]);
                ind = [runCond{S}{:}]; ind = {ind.labelAcq}; ind = ismember(ind,runSet{S}{end}.label);
                runSet{S}{end}.fList = [runCond{S}{ind}];
                runSet{S}{end}.fList  = cat(1,runSet{S}{end}.fList.fList);
                runSet{S}{end}.nDummy = cat(1,dummyList{S}{ind});
                if ~isempty(runSet{S}{end}.fList)
                    forceThis   = 0;
                    verboseThis = 1;
                    runSet{S}{end}.initFiles = initPreproc(runSet{S}{end},[],[],forceThis,verboseThis);

                    forceThis   = 0;
                    verboseThis = 1;
                    info.useSynth = 1;
                    fMask = volAnatPreproc2(do,info,runSet{S}{end},forceThis,verboseThis);
                    fMask = fMask.func.mask.brainInv.mri.fspec;

                    forceThis   = 0;
                    verboseThis = 1;
                    param.baseType = 'first'; % 'first' 'av' 'mcAv'
                    param.spSmFac  = 1;
                    fBase = [];
                    runSet{S}{end}.wrMocoFiles = estimMotionWR2(runSet{S}{end}.initFiles,param,fBase,fMask,forceThis,verboseThis);
                end


        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Between-run preprocessing %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% run or load preprocessing
        do.loadIt = 0;
        do.doIt = 1;
        do.saveIt = 0;

        switch info.dataSetLabel
            case 'vsmDriven'
                for rs = 1:length(runSet{S})
                    if isempty(runSet{S}{rs}.fList); continue; end
                    switch runSet{S}{rs}.label
                        case {'vfMRI' 'bold'}
                            forceThis   = 0;
                            verboseThis = 1;
                            param.baseType = 'firstRun_avFrame'; % 'first' 'av' 'mcAv'
                            param.spSmFac  = 1;
                            fBase = [];
                            fMask = cellstr(runSet{S}{rs}.wrMocoFiles.fMaskList);
                            if length(unique(fMask))==1
                                fMask = fMask{1}; else; dbstack; error('code that');
                            end
                            runSet{S}{rs}.brMocoFiles = estimMotionBR(runSet{S}{rs}.wrMocoFiles,fBase,fMask,param,forceThis,verboseThis);
                        otherwise
                            dbstack; error('code that');
                    end
                end
            otherwise
                dbstack; error('code that');
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Between-session preprocessing %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% run or load preprocessing
        do.loadIt = 0;
        do.doIt = 1;
        do.saveIt = 0;
        forceThis   = 0;
        verboseThis = 1;


        switch info.dataSetLabel
            case 'vsmDriven'
                for rs = 1:length(runSet{S})
                    if isempty(runSet{S}{rs}.fList); continue; end

                    if S>1
                        runSetRef = {};
                        for SS = 1:S-1
                            for rsrs = 1:length(runSet{SS})
                                if strcmp(runSet{SS}{rsrs}.sub,runSet{S}{rs}.sub) && ...
                                        str2num(runSet{SS}{rsrs}.ses)<str2num(runSet{S}{rs}.ses) && ...
                                        strcmp(runSet{SS}{rsrs}.label,runSet{S}{rs}.label) && ...
                                        ~isempty(runSet{SS}{rsrs}.fList)

                                    runSetRef{end+1} = runSet{SS}{rsrs};
                                end
                            end
                        end

                        if ~isempty(runSetRef)
                            runSetRef = runSetRef{1};

                            param.sourceType = 'avRun_avFrame';
                            param.spSmFac  = 1;
                            fBase = fullfile(runSetRef.brMocoFiles.wd,'av_cat_mcBR_av_mcWR_setPlumb_volTs.nii.gz');
                            fMask = cellstr(runSetRef.brMocoFiles.fMask);
                            if length(unique(fMask))==1
                                fMask = fMask{1}; else; dbstack; error('code that');
                            end
                            runSet{S}{rs}.bsMocoFiles = estimMotionBS2(runSet{S}{rs}.brMocoFiles,fBase,fMask,param,forceThis,verboseThis);
                            runSet{S}{rs}.bsMocoFiles.fGeomSes1 = runSetRef.initFiles.fGeom;
                        end
                    end
                    %
                    % switch runSet{S}{rs}.label
                    %     case {'vfMRI' 'bold'}
                    %         forceThis   = 0;
                    %         verboseThis = 1;
                    %         param.baseType = 'firstRun_avFrame'; % 'first' 'av' 'mcAv'
                    %         param.spSmFac  = 1;
                    %         fBase = [];
                    %         fMask = cellstr(runSet{S}{rs}.wrMocoFiles.fMaskList);
                    %         if length(unique(fMask))==1
                    %             fMask = fMask{1}; else; dbstack; error('code that');
                    %         end
                    %         runSet{S}{rs}.brMocoFiles = estimMotionBR(runSet{S}{rs}.wrMocoFiles,fBase,fMask,param,forceThis,verboseThis);
                    %     otherwise
                    %         dbstack; error('code that');
                    % end
                end
            otherwise
                dbstack; error('code that');
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Finalize preprocessing (one-step interpolation) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        do.loadIt = 0;
        do.doIt = 1;
        do.saveIt = 0;
        switch info.dataSetLabel
            case 'vsmDriven'
                for rs = 1:length(runSet{S})
                    if isempty(runSet{S}{rs}.fList); continue; end

                    forceThis   = 0;
                    verboseThis = 1;
                    initFiles    = runSet{S}{rs}.initFiles;
                    preprocFiles = {runSet{S}{rs}.wrMocoFiles runSet{S}{rs}.brMocoFiles};
                    if isfield(runSet{S}{rs},'bsMocoFiles')
                        preprocFiles = [preprocFiles {runSet{S}{rs}.bsMocoFiles}];
                    else
                        preprocFiles = [preprocFiles {[]}];
                    end

                    runSet{S}{rs}.finalFiles = finalizePreproc2(initFiles,preprocFiles,forceThis,verboseThis);
                end


                % fList = replace(runSet{S}.wrMocoFiles.fMocoList,'.nii.gz','.param.1D');
                % for f = 1:length(fList)
                %     cmd = {srcAfni};
                %     cmd{end+1} = '1dplot \';
                %     cmd{end+1} = '-yaxis "-1.5:1.5:1:1" - \';
                %     cmd{end+1} = ['-ynames `head -2 ' fList{f} ' | tail -1 | cut -c 3-` - \'];
                %     cmd{end+1} = fList{f};
                %     [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
                % end
            otherwise
                dbstack; error('code that')
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% For surface registration %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        forceThis = 0;

        %%% fs
        volAnat{S}{1,end+1}.sub      = subList{S};
        volAnat{S}{end}.ses          = sesList{S};
        volAnat{S}{end}.label        = 'fs';
        volAnat{S}{end}.dbDir        = sesDbList{S};
        fsDir = dir(fullfile(volAnat{S}{end}.dbDir,'fs',volAnat{S}{1,end}.sub));
        if isempty(memprage{S}.fList) || isempty(fsDir)
            volAnat{S}{end}.wd           = '';
            volAnat{S}{end}.bidsDir      = '';
            volAnat{S}{end}.bidsDerivDir = '';
            volAnat{S}{end}.fsDir        = '';
        else
            disp('fsDir found in dbDir')
            volAnat{S}{end}.wd           = outDir;
            volAnat{S}{end}.bidsDir      = bidsDirList{S};
            volAnat{S}{end}.bidsDerivDir = fullfile(bidsDirList{S},'derivatives',['set-' volAnat{S}{end}.label]);
            volAnat{S}{end}.fsDir        = fullfile(volAnat{S}{end}.bidsDerivDir,'fs');
            
            % copy to bids derivDir
            fsDir = fsDir(1).folder;
            if ~isFS(fsDir); dbstack; error('fsDir in dbDir does not seem to contain fs analysis'); end
            if forceThis || ~isFS(volAnat{S}{end}.fsDir)
                disp('copying fsDir from dbDir to bidsDerivDir')
                copyfile(fsDir,volAnat{S}{end}.fsDir)
            else
                disp('fsDir already in bidsDerivDir')
            end
        end

        %%% avMap
        volAnat{S}{1,end+1}.sub      = subList{S};
        volAnat{S}{end}.ses          = sesList{S};
        volAnat{S}{end}.label        = 'avMap';
        if isempty(avMap{S}.fList)
            volAnat{S}{end}.fList        = '';
            volAnat{S}{end}.dbDir        = '';
            volAnat{S}{end}.wd           = '';
            volAnat{S}{end}.bidsDir      = '';
            volAnat{S}{end}.bidsDerivDir = '';
        else
            volAnat{S}{end}.fList        = avMap{S}.fList;
            volAnat{S}{end}.dbDir        = sesDbList{S};
            volAnat{S}{end}.wd           = outDir;
            volAnat{S}{end}.bidsDir      = bidsDirList{S};
            volAnat{S}{end}.bidsDerivDir = fullfile(bidsDirList{S},'derivatives',['set-' volAnat{S}{end}.label]);
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end


    %%%%%%%%%%%%%%%%%%%%%%%
    %% avMap registration %
    %%%%%%%%%%%%%%%%%%%%%%%
    forceThis = 0;
    subList2 = unique(subList);
    avMap = cell(size(subList2));
    fs    = cell(size(subList2));
    vfMRI = cell(size(subList2));
    %%% get relevant data struct for each subject
    for S = 1:length(subList2)
        volAnat2 = volAnat(ismember(subList,subList2{S}));
        runSet2  = runSet(ismember(subList,subList2{S}));
        for s = 1:length(volAnat2)
            for ar = 1:length(volAnat2{s})
                %%% get avMap
                if strcmp(volAnat2{s}{ar}.label,'avMap') && ~isempty(volAnat2{s}{ar}.fList)
                    avMap{S}{end+1} = volAnat2{s}{ar};
                end
                %%% get fs
                if strcmp(volAnat2{s}{ar}.label,'fs') && ~isempty(volAnat2{s}{ar}.fsDir)
                    fs{S}{end+1} = volAnat2{s}{ar};
                end
            end
            for ar = 1:length(runSet2{s})
                %%% get vfMRI
                if strcmp(runSet2{s}{ar}.label,'vfMRI') && isfield(runSet2{s}{ar},'finalFiles') && strcmp(runSet2{s}{ar}.ses,'1')
                   vfMRI{S}{end+1} = runSet2{s}{ar}.finalFiles;
                end
            end
        end
        %%% make sure data struct make sense
        if length(fs{S})~=1; dbstack; error(['expecting a single data struct for fs of ' subList2{S}]); end
        fs{S} = fs{S}{1};
        if length(avMap{S})~=1; dbstack; error(['expecting a single data struct for avMap of ' subList2{S}]); end
        avMap{S} = avMap{S}{1};
        if length(vfMRI{S})~=1; dbstack; error(['expecting a single data struct for vfMRI of ' subList2{S}]); end
        vfMRI{S} = vfMRI{S}{1};
    end

    %%% write avMap into functional space without resampling
    for S = 1:length(avMap)
        fList = fullfile({avMap{S}.fList.folder},{avMap{S}.fList.name})';
        ref = fullfile(vfMRI{S}.bidsDerivDir,'sesAvCat_av_cat_av_preproc_volTs.nii.gz');
        fAvMap = fullfile(avMap{S}.bidsDerivDir,replace(avMap{S}.fList(1).name,'.nii.gz',''));
        fAvMap = replace(fAvMap,'echo-1','echo-cat'); if ~exist(fAvMap,'dir'); mkdir(fAvMap); end
        fAvMap = fullfile(fAvMap,'in-vfMRI.nii.gz');
        if forceThis || ~exist(fAvMap,'file')
            mriRef  = MRIread(ref,1);
            for i = 1:length(fList)
                mri = MRIread(fList{i});
                if i==1
                    mriRef.vol = mri.vol;
                    mriRef.volres = mri.volres;
                    mriRef.xsize  = mri.xsize;
                    mriRef.ysize  = mri.ysize;
                    mriRef.zsize  = mri.zsize;
                    mriRef.tr     = mri.tr;
                    mriRef.te     = mri.te;
                else
                    mriRef.vol = cat(4,mriRef.vol,mri.vol);
                end
            end
            MRIwrite(mriRef,fAvMap);
        end
        avMap{S}.fList_invfMRI = cellstr(fAvMap);
        for i = 1:length(volAnat2{S})
            if ~strcmp(volAnat2{S}{i}.label,'avMap'); continue; end
            volAnat2{S}{i}.fList_invfMRI = cellstr(fAvMap);
        end
    end

    % disp(strjoin([fullfile(vfMRI{S}.bidsDerivDir,'sesAvCat_av_cat_av_preproc_volTs.nii.gz')
    % fullfile({avMap{S}.fList.folder},{avMap{S}.fList.name})'
    % avMap{S}.fList_invfMRI],[' \\' newline]))
    
    % %%% Surface registration not working for avMap for a lack of WM to GM
    % %%% contrast. May be able to make it work somehow but too much work
    % [SUBJECTS_DIR,SUBJECT,~] = fileparts(fs{S}.fsDir);
    % MOV = {avMap{S}.fList.name}'; MOV = avMap{S}.fList(contains(MOV,'echo-1')); MOV = fullfile(MOV.folder,MOV.name);
    % [~,curOutDir,~] = fileparts(replace(MOV,'.nii.gz','')); curOutDir = fullfile(avMap{S}.bidsDerivDir,curOutDir);
    % surfReg(SUBJECTS_DIR,SUBJECT,MOV,curOutDir)
    for s = 1:length(runSet)
        for ar = 1:length(runSet{s})
            if ~strcmp(runSet{s}{ar}.label,'vfMRI'); continue; end
            sub = runSet{s}{ar}.sub;
            for S = 1:length(avMap)
                if ~strcmp(avMap{S}.sub,sub); continue; end
                runSet{s}{ar}.avMap.f = avMap{S}.fList_invfMRI;
            end
        end
    end
    %% %%%%%%%%%%%%%%%%%%%%%%



    %% %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Refactor from set to cond %
    %%%%%%%%%%%%%%%%%%%%%%%%%% %%
    [runCond,subList,runCondAcqList,runCondStimList] = set2cond2(runSet,runCond);
    save(mfilename,'runCond','subList','runCondAcqList','runCondStimList')
else
    load(mfilename,'runCond','subList','runCondAcqList','runCondStimList')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Anatomical processing (masks and rois) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do.loadIt = 0;
do.doIt   = 1;
do.saveIt = 0;
forceThis   = 0;
forceRoi    = 0;
verboseThis = 1;

for S = 1:size(runCond,1)
    info.sub = subList{S};
    for rca = 1:length(runCondAcqList)
        if ~isfield(runCond{S},runCondAcqList{rca}); continue; end
        out = volAnatPreproc3(do,info,runCond{S}.(runCondAcqList{rca}),forceThis,forceRoi);
        for rcs = 1:length(runCondStimList)
            if ~isfield(runCond{S}.(runCondAcqList{rca}),runCondStimList{rcs}); continue; end
            runCond{S}.(runCondAcqList{rca}).(runCondStimList{rcs}).volAnatSub = out;
        end
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%
% Subject loop
volTsAll = cell(size(subList));
volPsdAll = cell(size(subList));
allCmd = cell(size(subList));

sesIndList = 1:length(subList);
for s = 1:length(subList(sesIndList))
    S = sesIndList(s);

    % see doIt_pilotPipeline08c/stepTemplate.m
    %% Housekeeping
    info.sub          = subList{S};
    % info.ses          = sesList{S};
    % info.sesDb        = sesDbList{S};
    % info.bidsDir      = bidsDirList{S};
    % info.bidsDerivDir = fullfile(info.bidsDir,'derivatives'); if ~exist(info.bidsDerivDir,'dir'); mkdir(info.bidsDerivDir); end


    tic; disp(' '); disp(' '); disp(' ');
    tmp = {['sub-' info.sub]};
    tmp{end+1} = [num2str(s) '/' num2str(length(sesIndList))];
    tmp{end+1} = ['Nses=' num2str(length(subList))];
    tmp = strjoin(tmp,'; ');
    disp(repmat('+',1,length(tmp))); disp(repmat('+',1,length(tmp))); disp(repmat('+',1,length(tmp))); disp(tmp)
    %%






    %% Main analysis steps




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load data (and write the cross-run mean) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    forceThis = 0;

    copyFieldList = {'sub' 'ses' 'label' 'labelAcq' 'dsgn' 'wd' 'bidsDir' 'bidsDerivDir' 'ppLabelList' 'dataType' 'fOrigList' 'fPreprocList' 'fTransList' 'fTransCatList' 'bidsList' 'acqTime' 'nDummy'  'volAnat' 'volAnatSes' 'volAnatSub'};
    disp('loading volTs')
    switch info.dataSetLabel
        case 'vsmDriven'
            
            ar = 1;
            acqRunCond = 'vfMRI'; runCondAcqList{ar};
            for sr = 1:length(runCondStimList)
                stimRunCond = runCondStimList{sr};
                if ~isfield(runCond{S}.(acqRunCond),stimRunCond); continue; end
                nr = length(runCond{S}.(acqRunCond).(stimRunCond).fPreprocList);

                volTs = [];
                [nrPerSes, sesList] = groupcounts(cellstr(runCond{S}.(acqRunCond).(stimRunCond).ses));
                volTsSes = cell(size(sesList));
                volTsSub = [];

                fOrigList    = runCond{S}.(acqRunCond).(stimRunCond).fOrigList;
                fPreprocList = runCond{S}.(acqRunCond).(stimRunCond).fPreprocList;
                for rr = 1:nr
                    disp([' ' info.sub '; ' acqRunCond '; ' stimRunCond ' ' num2str(sr) '/' num2str(length(runCondStimList)) '; file ' num2str(rr) '/' num2str(nr) ' (' fPreprocList{rr} ')'])
                    volTs(rr,1).mri = MRIread(fPreprocList{rr});
                    for i = 1:length(copyFieldList)
                        if isfield(runCond{S}.(acqRunCond).(stimRunCond),copyFieldList{i})
                            if size(runCond{S}.(acqRunCond).(stimRunCond).(copyFieldList{i}),1)==1
                                volTs(rr,1).mri.(copyFieldList{i}) = runCond{S}.(acqRunCond).(stimRunCond).(copyFieldList{i});
                            else
                                volTs(rr,1).mri.(copyFieldList{i}) = runCond{S}.(acqRunCond).(stimRunCond).(copyFieldList{i})(rr);
                            end
                        end
                    end
                    nFrameOrig = MRIread(fOrigList{rr},1); nFrameOrig = nFrameOrig.nframes;
                    nFrame     = runCond{S}.(acqRunCond).(stimRunCond).nFrame(rr);
                    volTs(rr,1).mri.nDummyRemoved = nFrameOrig-nFrame;
                    if rr==1
                        volTsSub.mri     = volTs(rr).mri;
                    else
                        volTsSub.mri.vol = volTsSub.mri.vol...
                                           + volTs(rr).mri.vol;
                    end
                    % if rr==1
                    %     runCond{S}.(acqRunCond).(stimRunCond).volTsSub.mri     = runCond{S}.(acqRunCond).(stimRunCond).volTs(rr).mri;
                    % else
                    %     runCond{S}.(acqRunCond).(stimRunCond).volTsSub.mri.vol = runCond{S}.(acqRunCond).(stimRunCond).volTsSub.mri.vol...
                    %         + runCond{S}.(acqRunCond).(stimRunCond).volTs(rr).mri.vol;
                    % end


                    sesInd = ismember(sesList,runCond{S}.(acqRunCond).(stimRunCond).ses(rr));
                    if isempty(volTsSes{sesInd})
                        volTsSes{sesInd}.mri     = volTs(rr).mri;
                    else
                        volTsSes{sesInd}.mri.vol = volTsSes{sesInd}.mri.vol...
                                                   + volTs(rr).mri.vol;
                    end
                    
                end

                %divide by number of runs
                volTsSub.mri.vol = volTsSub.mri.vol ./ nr;
                % runCond{S}.(acqRunCond).(stimRunCond).volTsSub.mri.vol = runCond{S}.(acqRunCond).(stimRunCond).volTsSub.mri.vol ./ nr;
                volTsSes = cat(1,volTsSes{:});
                for ses = 1:size(volTsSes,1)
                    volTsSes(ses).mri.vol = volTsSes(ses).mri.vol ./ nrPerSes(ses);
                end
                

                %write averages
                for ses = 1:size(volTsSes,1)
                    mri = volTsSes(ses).mri;
                    f = mri.fspec; f = strsplit(f,filesep); f{end-1} = strsplit(f{end-1},'_'); f{end-1}{contains(f{end-1},'run-')} = 'run-av'; f{end-1} = strjoin(f{end-1},'_'); if ~exist(strjoin(f(1:end-1),filesep),'dir'); mkdir(strjoin(f(1:end-1),filesep)); end; f = strjoin(f,filesep);
                    if forceThis || ~exist(f,'file')
                        disp(['  writing ses-' sesList{ses} ' cross-run average (' f ')'])
                        MRIwrite(mri,f);
                    else
                        disp(['  ses-' sesList{ses} ' cross-run average already done, skipping (' f ')'])
                    end
                    mri.fspec = f;
                    volTsSes(ses).mri = mri; clear mri
                end

                mri = volTsSub.mri;
                f = mri.fspec; f = strsplit(f,filesep); f{end-1} = strsplit(f{end-1},'_'); f{end-1}{contains(f{end-1},'run-')} = 'run-av'; f{end-1}{contains(f{end-1},'ses-')} = 'ses-av'; f{end-1} = strjoin(f{end-1},'_'); if ~exist(strjoin(f(1:end-1),filesep),'dir'); mkdir(strjoin(f(1:end-1),filesep)); end; f = strjoin(f,filesep);
                f = strsplit(f,filesep); fTmp = f(1:end-2); fTmp{contains(fTmp,'ses-')} = 'ses-av'; f(1:end-2) = fTmp; if ~exist(strjoin(f(1:end-1),filesep),'dir'); mkdir(strjoin(f(1:end-1),filesep)); end; f = strjoin(f,filesep);
                if forceThis || ~exist(f,'file')
                    disp(['  writing cross-session cross-run average (' f ')'])
                    MRIwrite(mri,f);
                else
                    disp(['  cross-session cross-run average already done, skipping (' f ')'])
                end
                mri.fspec = f;
                volTsSub.mri = mri; clear mri

                
                %compile
                runCond{S}.(acqRunCond).(stimRunCond).volTs    = volTs; clear volTs
                runCond{S}.(acqRunCond).(stimRunCond).volTsSes = volTsSes; clear volTsSes
                runCond{S}.(acqRunCond).(stimRunCond).volTsSub = volTsSub; clear volTsSub

            end

        otherwise
            dbstack; error('code that');
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Event-related response %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    do.loadIt  = 0;
    do.doIt    = 1;
    do.saveIt  = 0;
    do.writeIt = 1;
    forceThis   = 0;
    verboseThis = 1;
    switch info.dataSetLabel
        case 'vsmDriven'
            ar = 1;
            acqRunCond = 'vfMRI'; runCondAcqList{ar};
            for sr = 1:length(runCondStimList)
                stimRunCond = runCondStimList{sr};
                if ~isfield(runCond{S}.(acqRunCond),stimRunCond); continue; end
                volAnat = runCond{S}.(acqRunCond).(stimRunCond).volAnatSub;

                %individual runs
                [volResp, ~, info] = volTsGetResp2(do,info,runCond{S}.(acqRunCond).(stimRunCond).volTs,[],volAnat,forceThis,verboseThis);
                for rr = 1:size(volResp,1)
                    runCond{S}.(acqRunCond).(stimRunCond).volTs(rr,1).volResp = volResp(rr,1);
                end
                clear volResp

                %within session average across runs
                for ses = 1:size(runCond{S}.(acqRunCond).(stimRunCond).volTsSes,1)
                    [runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses).volResp, ~, info] = volTsGetResp2(do,info,runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses),[],volAnat,forceThis,verboseThis);
                end

                %cross-session and cross-run average
                [runCond{S}.(acqRunCond).(stimRunCond).volTsSub.volResp, ~, info] = volTsGetResp2(do,info,runCond{S}.(acqRunCond).(stimRunCond).volTsSub,[],volAnat,forceThis,verboseThis);
                clear volAnat
            end

        otherwise
            dbstack; error('code that')
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%

    
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Visualize Event-related response %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 0
        % condLabel = 'task_10sPrd1sDur';
        % condLabel = 'task_15sPrd1sDur';
        % condLabel = 'task_50sPrd1sDur';
        % condLabel = 'task_50sPrd5sDur';
        % condLabel = 'task_20sPrd1sDur';

        S = 2;
        acqCondLabel  = 'vfMRI'; fields(runCond{S});

        stimCondLabel = ''; fields(runCond{S}.(acqCondLabel));
        cmd = fsCommand2(runCond{S}.(acqCondLabel),stimCondLabel);

        stimCondLabel = 'task_20sPrd1sDur';
        cmd = fsCommand2(runCond{S}.(acqCondLabel),stimCondLabel);


        S = 2;
        acqCondLabel  = 'vfMRI'; fields(runCond{S});

        stimCondLabel = ''; fields(runCond{S}.(acqCondLabel));
        cmd = fsCommand2(runCond{S},acqCondLabel,stimCondLabel);

        stimCondLabel = 'task_20sPrd1sDur';
        cmd = fsCommand(runCond{S},acqCondLabel,stimCondLabel);

        runCond{S}.(acqCondLabel).(stimCondLabel).volTs
        % mri = [runCond{S}.(acqCondLabel).(stimCondLabel).volTs.mri];
        % tmp = [runCond{S}.(acqCondLabel).(stimCondLabel).volTs.volResp];
        % tmp(1).ts.fspec
        % {mri.fspec}'



        %
        % tmp = MRIread(char(fRespFav));
        % figure('WindowStyle','docked');
        % imagesc(tmp.vol(:,:,:,1)); colorbar
        % tmp = MRIread(char(fRespF));
        % figure('WindowStyle','docked');
        % imagesc(tmp.vol(:,:,:,1)); colorbar
        %
        % tmp = MRIread(char(fRespAv));
        % figure('WindowStyle','docked');
        % imagesc(tmp.vol(:,:,:,3)); colorbar
    end




    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Full (space-time-spectral) multitaper analysis %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    do.loadIt  = 0;
    do.doIt    = 1;
    do.saveIt  = 0;
    do.writeIt = 1;
    switch info.dataSetLabel
        case 'vsmDriven'
            % info.K = kList{S};
            
            ar = 1;
            acqRunCond = 'vfMRI'; runCondAcqList{ar};
            
            for sr = 1:length(runCondStimList)
                stimRunCond = runCondStimList{sr};
                if ~isfield(runCond{S}.(acqRunCond),stimRunCond) || ~isfield(runCond{S}.(acqRunCond).(stimRunCond),'volTsSes') || isempty(runCond{S}.(acqRunCond).(stimRunCond).volTsSes)
                    continue; end
                volAnat = runCond{S}.(acqRunCond).(stimRunCond).volAnatSub;
                roiLabelList = [volAnat.roi{:}]; roiLabelList = {roiLabelList.label}';

                for ses = 1:size(runCond{S}.(acqRunCond).(stimRunCond).volTsSes,1)
                    volTs    = rmfield(runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1),'volResp');
                    info.ses = volTs.mri.ses;
                    if ~isfield(volTs,'dsgn'); volTs.dsgn = runCond{S}.(acqRunCond).(stimRunCond).dsgn; end
                    
                    switch stimRunCond
                        case {'task_10sPrd1sDur'
                                'task_15sPrd1sDur'}
                            info.win = inf;
                            info.K = 10;
                        case {'task_20sPrd1sDur'}
                            info.win = round(9 / (volTs.mri.tr/1000)); %info.win.*(volTs.mri.tr/1000)
                            info.win(2) = 1;
                            info.K = 3;
                        case {'task_50sPrd1sDur'
                                'task_50sPrd5sDur'}
                            info.win = round(15 / (volTs.mri.tr/1000)); %info.win.*(volTs.mri.tr/1000)
                            % info.win(2) = floor(info.win/2);
                            info.win(2) = 1;
                            info.K = 5;
                    end

                    % runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1).volPsd = volPsdFullMt3(do,info,volTs.mri,volAnat.mask{2}.vesselRefinedDil.f);
                    % runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1).volPsd = volPsdFullMt3(do,info,volTs.mri,volAnat.roi{ismember(roiLabel,'vesselCalcarine')}.f);
                    runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1).volPsd = volPsdFullMt3(do,info,volTs.mri,volAnat.roi{ismember(roiLabelList,'vesselCalcarine')});
                end
                clear volTs volAnat volResp volPsd roiLabelList
            end

            % %plot
            % for c = 1:size(volTs.volTs,2)
            %     plotSpecAll(volTs.volTs(c).volPsd,volTs.volTs(c),[])
            %     tryL2svd(volTs.volTs(c).volPsd)
            % end
            % 
            % % I = 1;
            % % load(subList{I},'volTs')
            % % volTs = rmfield(volTs,{'volTs' 'volResp'});
            % % save(subList{I},'volTs','-v7.3')
            % %
            
            continue

        otherwise

%%
            % volPsd = volPsdFullMtAna(do,info,volTs,volAnat);
            info.onsetList = onsets;
            info.ondurList = durs;
            info.skipSvd  = 0;
            info.K        = 3;
            info.win      = 7; % should be at least K*2
            info.perm     = 0;

            if ~isfield(info.funcGroup,'condLabel')
                volPsd = volPsdFullMt2(do,info,volTs,volAnat);
                volPsd = volPsdFullMt2(do,info,volTs,volAnat,volPsd);
                % allCmd{I} = strjoin(viewResp(volTs.resp,volAnat,{volPsd.psd.fspec volPsd.svd.fspec.spSVmag},{':colormap=turbo' ':colormap=turbo'}),newline);
                % % viewResp(volTs.resp,volAnat,{volPsd.psd.fspec volPsd.svd.fspec.spSVmag volPsd.svd.fspec.spSVmag_fdrThresh volPsd.svd.fspec.spSVphase_fdrThresh},{':colormap=turbo' ':colormap=turbo' ':colormap=turbo' ':colormap=turbo'});
                % tryL2svd(volPsd)
                plotSpecAll(volPsd,volTs)
                drawnow

                continue

                % I=3;
                % clipboard('copy',allCmd{I})
                x = 245; %cursor_info.Position(1);
                y = 170; %cursor_info.Position(2);
                % x = 149; %cursor_info.Position(1);
                % y = 193; %cursor_info.Position(2);

                mriResp = MRIread(fullfile(info.dbDir,[info.dataSetLabel '_' info.sub],'mri','fmri_long_ushoot','L_resamp_e0.nii.gz'));

                figure('WindowStyle','docked');
                imagesc(std(volPsd.resp.tsOnBase.vol,[],4)); colorbar
                % imagesc(mean(volPsd.resp.tsOnBase.vol,4)); colorbar

                figure('WindowStyle','docked');
                imagesc(std(mriResp.vol,[],4)); colorbar
                % imagesc(mean(mriResp.vol,4)); colorbar

                figure('WindowStyle','docked');
                tr = volPsd.resp.tsOnBase.tr/1000;
                t = 0:tr:(volPsd.resp.tsOnBase.nframes-1).*tr;
                plot(t,smooth(squeeze(volPsd.resp.tsOnBase.vol(y,x,:,:)),4)); hold on
                tr = 0.2;
                t = -1:tr:((mriResp.nframes-1).*tr-1);
                plot(t/1.5,squeeze(mriResp.vol(y,x,:,:)));


            else
                for condInd = 1:length(info.funcGroup.condLabel)
                    tmpInd = info.funcGroup.condInd==condInd;
                    volTsMean(condInd) = volTs(find(tmpInd,1));
                    volTsMean(condInd).vol = mean(cat(5,volTs(info.funcGroup.condInd==condInd).vol),5);
                    volPsdMean(condInd) = volPsdFullMt2(do,info,volTsMean(condInd),volAnat(1));
                    tryL2svd(volPsdMean(condInd))
                end


            end



            continue


            info.win = 8;
            volPsd = volPsdFullMt2(do,info,volTs,volAnat);




            volPsd = volPsdFullMt2(do,info,volTs,volAnat);
            tryL2svd(volPsd(1))



            viewResp(volTs.resp,volAnat,{volPsd.psd.fspec volPsd.svd.fspec.spSVmag volPsd.svd.fspec.spSVmag_fdrThresh volPsd.svd.fspec.spSVphase_fdrThresh},{':colormap=turbo' ':colormap=turbo' ':colormap=turbo' ':colormap=turbo'});


            F = figure('WindowStyle','docked');
            Ht = tiledlayout(8,1); Ht.TileSpacing = 'tight'; Ht.Padding = 'tight';
            ax = {};
            ax{end+1} = plotTs2(volTs,Ht);
            ax{end+1} = plotCoh3(volPsd,Ht);


            F = figure('WindowStyle','docked');
            Ht = tiledlayout(12,1); Ht.TileSpacing = 'tight'; Ht.Padding = 'tight';
            ax = {};
            for roiInd = 1:12
                ax{end+1} = plotTs2(volTs,Ht,[],[],volAnat,roiInd);
            end

            midResp = 35;
            F = figure('WindowStyle','docked');
            Ht = tiledlayout(8,1); Ht.TileSpacing = 'tight'; Ht.Padding = 'tight';
            ax = {};
            ax{end+1} = plotTs2(volTs,Ht,[],[],volAnat.fun.roi.vesselCenter);
            ax{end+1} = plotPsd2(volPsd,Ht,midResp);
            ax{end+1} = plotPsdGram2(volPsd,Ht);
            ax{end+1} = plotPsdTrialGram2(volPsd,Ht,'av');
            ax{end+1} = plotPsdTrialGram2(volPsd,Ht,'pc');
            nexttile
            ax{end+1} = plotPsdTrialGram2(volPsd,Ht,'MD');
            ax{end+1} = plotPsdTrialGram2(volPsd,Ht,'MDpc');
            xline(ax{end},midResp)

            F = figure('WindowStyle','docked');
            Ht = tiledlayout(8,1); Ht.TileSpacing = 'tight'; Ht.Padding = 'tight';
            ax = {};
            ax{end+1} = plotTs2(volTs,Ht);
            ax{end+1} = plotCoh2(volPsd,Ht);
            ax{end+1} = plotCohGram2(volPsd,Ht);
            ax{end+1} = plotCohTrialGram2(volPsd,Ht,'av');
            ax{end+1} = plotCohTrialGram2(volPsd,Ht,'pc');
            ax{end+1} = plotCohTrialGram2(volPsd,Ht,'ek');
            ax{end+1} = plotCohTrialGram2(volPsd,Ht,'MD');
            ax{end+1} = plotCohTrialGram2(volPsd,Ht,'MDpc');



            %
            % F = figure('WindowStyle','docked');
            % Ht = tiledlayout(7,1); Ht.TileSpacing = 'tight'; Ht.Padding = 'tight';
            % ax = {};
            % ax{end+1} = plotTs2(volTs,Ht);
            % ax{end+1} = plotPsd2(volPsd,Ht);
            % ax{end+1} = plotPsdGram2(volPsd,Ht);
            % ax{end+1} = plotPsdTrialGram2(volPsd,Ht,'av');
            % ax{end+1} = plotPsdTrialGram2(volPsd,Ht,'pc');
            % ax{end+1} = plotPsdTrialGram2(volPsd,Ht,'MD');
            % ax{end+1} = plotPsdTrialGram2(volPsd,Ht,'MDpc');
            %
            % ax{end+1} = plotCoh2(volPsd,Ht);
            % ax{end+1} = plotCohGram2(volPsd,Ht);
            % ax{end+1} = plotCohTrialGram2(volPsd,Ht,'av');
            % % ax{end+1} = plotCohTrialGram2(volPsd,Ht,'pc');
            % ax{end+1} = plotCohTrialGram2(volPsd,Ht,'ek');
            % % ax{end+1} = plotCohTrialGram2(volPsd,Ht,'MD');
            % ax{end+1} = plotCohTrialGram2(volPsd,Ht,'MDpc');
            %
            % plotTs2(volTs,[]);
            % plotPsdTrialGram2(volPsd,[],'MD');
            % plotPsdTrialGram2(volPsd,[],'MDpc');
            % plotCohTrialGram2(volPsd,[],'MD');
            % plotCohTrialGram2(volPsd,[],'MDpc');

            return


            % Fs = 10;
            % tEnd = 5*60;
            % t = 0:1/Fs:tEnd;
            % dur = 40;
            % fStim1 = 1/(dur+20);
            % trialOnsets = 10:1/fStim1:tEnd; trialOnsets = trialOnsets(1:4);
            % missingInd1 = true(size(t));
            % for i = 1:length(trialOnsets)
            %     missingInd1(t>=trialOnsets(i) & t<trialOnsets(i)+dur) = false;
            % end
            % fStim2 = 1/(dur+1);
            % trialOnsets = 10:1/fStim2:tEnd; trialOnsets = trialOnsets(1:4);
            % missingInd2 = true(size(t));
            % for i = 1:length(trialOnsets)
            %     missingInd2(t>=trialOnsets(i) & t<trialOnsets(i)+dur) = false;
            % end
            % figure('WindowStyle','docked');
            % plot(t,missingInd1-0.5); hold on
            % plot(t,(missingInd2-0.5)*2);
            % ylim([-2 2])
            %
            % length(t(~missingInd1)) * 1/Fs == length(t(~missingInd2)) * 1/Fs
            % T = length(t(~missingInd1));
            % K = 5;
            % [TW,W,K] = K2W(T,K);
            %
            % [lambda1,u1] = MDslepian(W,K,t(~missingInd1),Fs);
            % [lambda2,u2] = MDslepian(W,K,t(~missingInd2),Fs);
            % figure('WindowStyle','docked');
            % subplot(2,1,1)
            % plot(t(~missingInd2),u1)
            % subplot(2,1,2)
            % plot(t(~missingInd2),u2)

            %%

    end







    %% Housekeeping
    disp(' '); tmp = ['sub-' info.sub '_ses-' info.ses ': DONE']; toc;
    disp(tmp); disp(repmat('+',1,length(tmp))); disp(repmat('+',1,length(tmp))); disp(repmat('+',1,length(tmp))); disp(' ');




    return

    %% Simulate stimulus response
    if ~isempty(onsetsList)
        volPsdSim = simPsd(onsetsList{S},ondursList{S},volTs,info,0.001);

        axes(axSingleVox);
        yyaxis right
        plot(squeeze(volPsdSim.f),volPsdSim.vec)
        axSingleVox.YScale = 'log';

        axes(axAv);
        yyaxis right
        plot(squeeze(volPsdSim.f),volPsdSim.vec)
        axAv.YScale = 'log';

    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Arbitrary simulations %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    figure('WindowStyle','docked');
    dursListList = [];

    for s = 1:1:10
        TR = 3;
        fStim = 1/(TR*20);
        fDur  = 1/(TR*s);
        nStim = 5;
        T = 1/fStim * (nStim+1);
        onsetsList = (0:1/fStim:1/fStim*nStim)';
        ondursList = ones(size(onsetsList)) .* 1/fDur;

        volTsTmp = vec2vol(volTs); volTsTmp.vol = [];
        volTsTmp.vol = [];
        volTsTmp.tr = TR*1000;
        volTsTmp.nframes = T/TR;
        volPsdSim = simPsd(onsetsList,ondursList,volTsTmp,info);

        plot3(squeeze(volPsdSim.f),volPsdSim.vec,ones(size(volPsdSim.vec)).*ondursList(1)); hold on

        dursListList(end+1) = ondursList(1);
    end
    legend(num2str(dursListList'))
    ax = gca;
    ax.YScale = 'log';
    %% %%%%%%%%%%%%%%%%%%%%%%%

end

return


S = 2;
acqRunCond  = 'vfMRI'; %ar = 1; runCondAcqList{ar};
stimRunCond = 'task_50sPrd5sDur'; %sr = 1; runCondStimList{sr};
% {'task_10sPrd1sDur'}
% {'task_15sPrd1sDur'}
% {'task_20sPrd1sDur'}
% {'task_50sPrd1sDur'}
% {'task_50sPrd5sDur'}
% for sr = 1:length(runCondStimList)
% stimRunCond = runCondStimList{sr};
volAnat = runCond{S}.(acqRunCond).(stimRunCond).volAnatSub;
roiLabelList = [volAnat.roi{:}]; roiLabelList = {roiLabelList.label}';
roiLabel = 'vesselCalcarine';
roiInd = ismember(roiLabelList,roiLabel);


for ses = 1:size(runCond{S}.(acqRunCond).(stimRunCond).volTsSes,1)
    volTs    = rmfield(runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1),'volResp');
    info.ses = volTs.mri.ses;
    if ~isfield(volTs,'dsgn'); volTs.dsgn = runCond{S}.(acqRunCond).(stimRunCond).dsgn; end

    switch stimRunCond
        case {'task_10sPrd1sDur'
                'task_15sPrd1sDur'}
            info.win = inf;
            info.K = 10;
        case {'task_20sPrd1sDur'}
            info.win = round(9 / (volTs.mri.tr/1000)); %info.win.*(volTs.mri.tr/1000)
            info.win(2) = 1;
            info.K = 3;
        case {'task_50sPrd1sDur'
                'task_50sPrd5sDur'}
            info.win = round(15 / (volTs.mri.tr/1000)); %info.win.*(volTs.mri.tr/1000)
            % info.win(2) = floor(info.win/2);
            info.win(2) = 1;
            info.K = 5;
    end

    % runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1).volPsd = volPsdFullMt3(do,info,volTs.mri,volAnat.mask{2}.vesselRefinedDil.f);
    % runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1).volPsd = volPsdFullMt3(do,info,volTs.mri,volAnat.roi{ismember(roiLabel,'vesselCalcarine')}.f);
    runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1).volPsd = volPsdFullMt3(do,info,volTs.mri,volAnat.roi{roiInd});
end
% cmd = fsCommand2(runCond{S}.(acqRunCond),stimRunCond);
plotSpecAll2(...
    runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1).volPsd,...
    volTs.mri,...
    runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1).volResp,...
    [],...
    runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1).volPsd.volAnat.f,...
    runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses,1).volPsd.volAnat.funcUlay,1);




%%%%%%%%%%%%%%%%%%
%% Visualize all %
%%%%%%%%%%%%%%%%%%

% mri1 = MRIread(runCond{1}.(acqRunCond).(runCondStim).volAnatSub.mask{2}.vesselRefinedDil.f);
% mri2 = MRIread(runCond{2}.(acqRunCond).(runCondStim).volAnatSub.mask{2}.vesselRefinedDil.f);
% figure('WindowStyle','docked');
% imagesc(mri2.vol)

ar = 1;
acqRunCond = 'vfMRI'; runCondAcqList{ar};
stimRunCond = 'task_50sPrd1sDur';
for S = 1:2
    for ses = 1:length(runCond{S}.(acqRunCond).(stimRunCond).volTsSes)
        volAnat = runCond{S}.(acqRunCond).(stimRunCond).volAnatSub;
        volTs   = runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses).mri;
        volPsd  = runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses).volPsd;
        volResp = runCond{S}.(acqRunCond).(stimRunCond).volTsSes(ses).volResp;
        ulay = strsplit(volTs.fspec,filesep); ulay(end-1) = []; ulay{end} = ['task-' volPsd.label '_av_cat_av_preproc_volTs.nii.gz']; ulay = strjoin(ulay,filesep);
        roiLabelList = [volAnat.roi{:}]; roiLabelList = {roiLabelList.label}';
        

        Q = 1;
        fMask = volAnat.roi{ismember(roiLabelList,'vesselCalcarine')}.f
        plotSpecAll2(volPsd,volTs,volResp,[],volAnat.mask{2}.vesselRefinedDil.f,ulay,Q)


        % plotSpecAll(volTs.volTs(c).volPsd,volTs.volTs(c),[])
        % tryL2svd(volTs.volTs(c).volPsd)
    end
end

for sr = 1:length(runCondStimList)
    stimRunCond = runCondStimList{sr};
    if ~isfield(runCond{S}.(acqRunCond),stimRunCond) || ~isfield(runCond{S}.(acqRunCond).(stimRunCond),'volTsSes') || isempty(runCond{S}.(acqRunCond).(stimRunCond).volTsSes)
        continue; end
end

runCond{S}



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individual vessel ROIs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 1;
info.sub = subList{S}; info.ses = '';
ar = 1; acqRunCond = 'vfMRI'; runCondAcqList{ar};
out = volAnatPreproc3(do,info,runCond{S}.(runCondAcqList{rca}),[],[]);




%% %%%%%%%%%%%%%%%



return

%% Write to nii
if 0
    volPsd.tr = volPsd.tr/1000;
    origFile = fullfile(info.dbDir,info.sesDb,info.dbFile);

    MRIwrite(vec2vol(volPsd),replace(origFile,'_ravg.nii.gz','_ravgPsd.nii.gz'));
    disp(['psd in: ' replace(origFile,'_ravg.nii.gz','_ravgPsd.nii.gz')])
    bidsFile = volPsd.fspec;
    MRIwrite(vec2vol(volPsd),replace(bidsFile,'_angio.nii.gz','_angioPsd.nii.gz'));
    disp(['psd in: ' replace(bidsFile,'_angio.nii.gz','_angioPsd.nii.gz')])

    volPsd.vec = log(volPsd.vec);
    MRIwrite(vec2vol(volPsd),replace(origFile,'_ravg.nii.gz','_ravgLogPsd.nii.gz'));
    disp(['logPsd in: ' replace(origFile,'_ravg.nii.gz','_ravgLogPsd.nii.gz')])
    bidsFile = volPsd.fspec;
    MRIwrite(vec2vol(volPsd),replace(bidsFile,'_angio.nii.gz','_angioLogPsd.nii.gz'));
    disp(['logPsd in: ' replace(bidsFile,'_angio.nii.gz','_angioLogPsd.nii.gz')])
end

% %% Agregate across subjects
% volTsAll{I} = volTs;
% volAnatAll{I} = volAnat;
% volPsdAll{I} = volPsd;


%% Subject loop: DONE
%%%%%%%%%%%%%%%%%%%%%





%% Housekeeping
saveIt = 0;
outVar = 'volTs';
stepLabel = 'preprocessing';
% WARNING: outVar = outVarAll; clear outVarAll
if exist([outVar 'All'],'var'); eval([outVar ' = ' outVar 'All; clear ' outVar 'All']); end; warning([newline outVar ' = ' outVar 'All; clear ' outVar 'All']);
stepLabel = [stepLabel '; Nsubj=' num2str(length(sesIndList)) '/' num2str(length(subList))];
stepFile = fullfile(outDir,['group-subj_' outVar]);
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

saveIt = 0;
outVar = 'volAnat';
stepLabel = 'anat processing';
% WARNING: outVar = outVarAll; clear outVarAll
if exist([outVar 'All'],'var'); eval([outVar ' = ' outVar 'All; clear ' outVar 'All']); end; warning([newline outVar ' = ' outVar 'All; clear ' outVar 'All']);
stepLabel = [stepLabel '; Nsubj=' num2str(length(sesIndList)) '/' num2str(length(subList))];
stepFile = fullfile(outDir,['group-subj_' outVar]);
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

saveIt = 0;
outVar = 'volPsd';
stepLabel = 'full mt analysis';
% WARNING: outVar = outVarAll; clear outVarAll
if exist([outVar 'All'],'var'); eval([outVar ' = ' outVar 'All; clear ' outVar 'All']); end; warning([newline outVar ' = ' outVar 'All; clear ' outVar 'All']);
stepLabel = [stepLabel '; Nsubj=' num2str(length(sesIndList)) '/' num2str(length(subList))];
stepFile = fullfile(outDir,['group-subj_' outVar]);
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

stepLabel = 'all data loaded';
disp(' '); disp(' '); disp(repmat('x',1,length(stepLabel))); disp(upper(stepLabel)); toc; disp(repmat('x',1,length(stepLabel))); disp(' '); disp(' ');



return
%%%%%%%%%%%%%%%%
%% Visualization
%% Plot All
subList
sesIndList
S = 3;
cropMask = volAnat{S}.fun.mask.head;

% IAll = {}; tfAll = {}; voxAll = {};

% volTs{I}.t(:,:,:,:,:,2)
[tf,vox] = plotAll(volPsd{S},volTs{S},[],[],cropMask,pipId);
IAll{end+1} = S; tfAll{end+1} = tf; voxAll{end+1} = vox;

s = 1; %length(IAll);
S = IAll{s}; tf = tfAll{s}; vox = voxAll{s};
plotAll(volPsd{S},volTs{S},tf,vox,cropMask,pipId);
