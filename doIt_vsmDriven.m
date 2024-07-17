clear all
close all
[outDir,pipId] = fileparts(mfilename('fullpath'));
outDir = fullfile(outDir,pipId); if ~exist(outDir,'dir'); mkdir(outDir); end
%%%%%%%%%%%%%%
%% Preparation
%%%%%%%%%%%%%%
%% Dependencies
% matlab
addpath(genpath(fullfile(pwd,pipId)))
addpath(genpath('/usr/local/freesurfer/stable7.4.1/matlab/'))
% addpath(genpath('/space/takoyaki/1/users/proulxs/tools/BrewerMap'))
addpath(genpath('/space/takoyaki/1/users/proulxs/tools/chronux'))
addpath(genpath('/space/takoyaki/1/users/proulxs/tools/vasomoTools'))
addpath(genpath('/space/takoyaki/1/users/proulxs/tools/bassReg'))
addpath(genpath('/space/takoyaki/1/users/proulxs/tools/martinosTools'))
% addpath(genpath('/space/takoyaki/1/users/proulxs/tools/spm12'))
% addpath(genpath('/cluster/freesurfer/unwarp/gradient_nonlin_unwarp')) % Polimeni's gradient distortion correction routine
% bash
global srcAfni srcFs
% srcFs = 'source /usr/local/freesurfer/nmr-stable741-env-bash';
srcFs = 'source /usr/local/freesurfer/fs-stable741-env-autoselect';
srcAfni = 'export PATH=$PATH:/usr/pubsw/packages/AFNI/23.1.05';



info.dataSetLabel = 'svfmri_long_short'; % 'svfmri_long_short' 'satin'



%% Variables and Paths
switch info.dataSetLabel
    
    % Seb satin
    case 'satin'
        info.preprocDir = fullfile('/autofs/space/takoyaki_001/users/proulxs/others/div/vasomo',pipId);
        % info.preprocDir = '/autofs/space/takoyaki_001/users/proulxs/vasomo/vasomoAna/vasomoInflow/pilotPipeline07';
        info.pipId   = pipId;
        info.bidsDir = fullfile('/space/takoyaki/1/users/proulxs/seqDev/satin/bids');
        info.bidsFile = [];
        %%% Data location
        info.dbDir   = '/space/takoyaki/1/users/proulxs/seqDev/satin/doIt_satinPip01/sub-pilot01/ses-1/func';
        info.dbFile  = 'sub-pilot01_ses-1_echo-1_run-1_task-vis_acq-satOFFg2t10l6f26_bold/preproc/preproc_volTs.nii.gz'; %'L_resamp_e0.nii.gz';
        info.runCond = {'set-f0b0' 'set-f1b1'};
        % tmp = dir(fullfile(info.dbDir)); disp({tmp.name}')
        % tmp = dir(fullfile(info.dbDir,'set-f0b0')); disp({tmp.name}')
        % tmp = dir(fullfile(info.dbDir,'set-f0b0',info.dbFile)); disp({tmp.name}')
        % tmp = dir(fullfile(info.dbDir,info.runCond{2},info.dbFile)); disp({tmp.name}')

    
    % Divya Long and Short
    case 'svfmri_long_short'
        info.preprocDir = fullfile('/autofs/space/takoyaki_001/users/proulxs/others/div/vasomo',pipId);
        % info.preprocDir = '/autofs/space/takoyaki_001/users/proulxs/vasomo/vasomoAna/vasomoInflow/pilotPipeline07';
        info.pipId   = pipId;
        info.bidsDir = fullfile('/autofs/space/takoyaki_001/users/proulxs/others/div/vasomo/bids');
        %%% Data location
        info.dbDir   = '/autofs/cluster/meso/users/divya/single_vessel_fMRI_long_short';
        info.dbFile  = 'L_echo0_ravg.nii.gz'; %'L_resamp_e0.nii.gz';
        tmp = dir(fullfile(info.dbDir,'svfmri_long_short_*')); disp({tmp.name}')
        % Timing information - 0.2 secs resampled response.
        % tvec = -1:0.2:56;
        % tvec = 0 is the stimulus on for 16 secs, then off.
end



subList       = {};
sesList       = {};
sesDbList     = {};
sesDbFileList = {};
volTrList        = {};
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

    % Seb satin
    case 'satin'
        info.dataSetLabel = 'satin';
        disp('/space/takoyaki/1/users/proulxs/seqDev/satin/doIt_satinPip01/sub-pilot01/ses-1/func/set-f0b0/sub-pilot01_ses-1_echo-1_run-1_task-vis_acq-satOFFg2t10l6f26_bold/preproc')
        disp(fullfile(info.dbDir,info.runCond{2}))
        dir(fullfile(info.dbDir,info.runCond{2}))
        bidsStr = 'sub-pilot01_ses-1_echo-1_run-1_task-vis_acq-satONg2t10l6f26_bold';
        info.dbFile = fullfile('preproc','preproc_volTs.nii.gz');
        dir(fullfile(info.dbDir,info.runCond{2},bidsStr,info.dbFile));
        % /space/takoyaki/1/users/proulxs/seqDev/satin/doIt_satinPip01/sub-pilot01/ses-1/func/set-f0b0
        subList{end+1}       = 'pilot01';
        sesList{end+1}       = '1';
        %/autofs/space/takoyaki_001/users/proulxs/vasomo/source/expDb/satInflow_pilot01/2023-10-29--bay2--satInflow_pilot01/bids/ses-1/func
        % sesDbList{end+1}     = '/space/takoyaki/1/users/proulxs/seqDev/satin/doIt_satinPip01/sub-pilot01/ses-1/func/set-f0b0';
        sesDbList{end+1}     = fullfile(info.dbDir,info.runCond{2},bidsStr);
        sesDbFileList{end+1} = fullfile('preproc','preproc_volTs.nii.gz');
        runCondList{end+1}   = info.runCond;
        exampleFile = fullfile(info.dbDir,info.runCond{2},bidsStr,info.dbFile);
        volTs = MRIread(exampleFile,1);
        matSz = 400; PAT = 4; refLine = 0;  rfTR = [];
        volTrList{end+1} = volTs.tr/1000; %flashVolTr(rfTR,matSz,PAT,refLine);
        % volTs = MRIread(exampleFile,1);
        dummyList{end+1}     = 1;
        kList{end+1}         = 6;
        disp([volTs.nframes volTs.tr/1000]); clear volTs
        winLengthList{end+1} = 8; % [number of frames]
        winStepList{end+1}   = 1;   % [number of frames]
        disp([winLengthList{end} winStepList{end}])
        disp([winLengthList{end} winStepList{end}].*volTrList{end})
        dir(fullfile(info.dbDir,info.runCond{2},bidsStr,'preproc/cond-visOn_model-custom_startTime.1D'))
        onsetsList{end+1}    = [21.000 78.000 135.000 192.000]';
        ondursList{end+1}    = [18.000 18.000  18.000  18.000]';
        dtStimList{end+1}    = 3;
        labelList{end+1}     = {'vis'};
        extraList{end+1}     = '';

    % Divya Long and Short
    case 'svfmri_long_short'
        info.dataSetLabel    = 'svfmri_long_short';
        subList{end+1}       = 'svfmri_long_short_pilot07';
        sesList{end+1}       = '1';
        sesDbList{end+1}     = 'svfmri_long_short_pilot07/mri/fmri_long_ushoot';
        matSz = 400; refLine = 48; PAT = 4; rfTR = 0.029; volTs = MRIread(fullfile(info.dbDir,sesDbList{end},info.dbFile),1);
        volTrList{end+1}        = (refLine+(matSz-refLine)/PAT)*rfTR;      if isempty(volTrList{end}); volTrList{end} = volTs.tr/1000; else; volTs.tr = volTrList{end}*1000; end
        dummyList{end+1}     = 1;
        kList{end+1}         = 6;
        disp([volTs.nframes volTs.tr/1000])
        winLengthList{end+1} = 8; % [number of frames]
        winStepList{end+1}   = 1;   % [number of frames]
        disp([winLengthList{end} winStepList{end}])
        disp([winLengthList{end} winStepList{end}].*volTrList{end})
        onsetsList{end+1}    = [24.0 80.8 137.6 194.4 251.2]';
        ondursList{end+1}    = [16.0 16.0  16.0  16.0  16.0]';
        dtStimList{end+1}    = 4;
        labelList{end+1}     = {'vis'};
        extraList{end+1}     = dir('/autofs/cluster/meso/users/divya/single_vessel_fMRI_long_short/svfmri_long_short_pilot07/mri/fmri_long/lbl/*.nii.gz');
        % /autofs/cluster/meso/users/divya/single_vessel_fMRI_long_short/svfmri_long_short_pilot07/mri/fmri_long/lbl


        subList{end+1}       = 'svfmri_long_short_pilot06';
        sesList{end+1}       = '1';
        sesDbList{end+1}     = 'svfmri_long_short_pilot06/mri/fmri_long_ushoot';
        matSz = 400; refLine = 48; PAT = 4; rfTR = 0.029; volTs = MRIread(fullfile(info.dbDir,sesDbList{end},info.dbFile),1);
        volTrList{end+1}        = (refLine+(matSz-refLine)/PAT)*rfTR;      if isempty(volTrList{end}); volTrList{end} = volTs.tr/1000; else; volTs.tr = volTrList{end}*1000; end
        dummyList{end+1}     = 1;
        kList{end+1}         = 6;
        disp([volTs.nframes volTs.tr/1000])
        winLengthList{end+1} = 8; % [number of frames]
        winStepList{end+1}   = 1;   % [number of frames]
        disp([winLengthList{end} winStepList{end}])
        disp([winLengthList{end} winStepList{end}].*volTrList{end})
        onsetsList{end+1}    = [24.0 80.8 137.6 194.4 251.2]';
        ondursList{end+1}    = [16.0 16.0  16.0  16.0  16.0]';
        dtStimList{end+1}    = 4;
        labelList{end+1}     = {'vis'};
        extraList{end+1}     = dir('/autofs/cluster/meso/users/divya/single_vessel_fMRI_long_short/svfmri_long_short_pilot06/mri/fmri_long/lbl/*.nii.gz');

        subList{end+1}       = 'svfmri_long_short_pilot04';
        sesList{end+1}       = '1';
        sesDbList{end+1}     = 'svfmri_long_short_pilot04/mri/fmri_long_ushoot';
        matSz = 400; refLine = 48; PAT = 4; rfTR = 0.029; volTs = MRIread(fullfile(info.dbDir,sesDbList{end},info.dbFile),1);
        volTrList{end+1}        = (refLine+(matSz-refLine)/PAT)*rfTR;      if isempty(volTrList{end}); volTrList{end} = volTs.tr/1000; else; volTs.tr = volTrList{end}*1000; end
        dummyList{end+1}     = 1;
        kList{end+1}         = 6;
        disp([volTs.nframes volTs.tr/1000])
        winLengthList{end+1} = 8; % [number of frames]
        winStepList{end+1}   = 1;   % [number of frames]
        disp([winLengthList{end} winStepList{end}])
        disp([winLengthList{end} winStepList{end}].*volTrList{end})
        onsetsList{end+1}    = [24.0 80.8 137.6 194.4 251.2]';
        ondursList{end+1}    = [16.0 16.0  16.0  16.0  16.0]';
        dtStimList{end+1}    = 4;
        labelList{end+1}     = {'vis'};
        extraList{end+1}     = dir('/autofs/cluster/meso/users/divya/single_vessel_fMRI_long_short/svfmri_long_short_pilot04/mri/fmri_long/lbl/*.nii.gz');

        subList{end+1}       = 'svfmri_long_short_pilot03';
        sesList{end+1}       = '1';
        sesDbList{end+1}     = 'svfmri_long_short_pilot03/mri/fmri_long_ushoot';
        matSz = 400; refLine = 48; PAT = 4; rfTR = 0.029; volTs = MRIread(fullfile(info.dbDir,sesDbList{end},info.dbFile),1);
        volTrList{end+1}        = (refLine+(matSz-refLine)/PAT)*rfTR;      if isempty(volTrList{end}); volTrList{end} = volTs.tr/1000; else; volTs.tr = volTrList{end}*1000; end
        dummyList{end+1}     = 1;
        kList{end+1}         = 6;
        disp([volTs.nframes volTs.tr/1000])
        winLengthList{end+1} = 8; % [number of frames]
        winStepList{end+1}   = 1;   % [number of frames]
        disp([winLengthList{end} winStepList{end}])
        disp([winLengthList{end} winStepList{end}].*volTrList{end})
        onsetsList{end+1}    = [24.0 80.8 137.6 194.4 251.2]';
        ondursList{end+1}    = [16.0 16.0  16.0  16.0  16.0]';
        dtStimList{end+1}    = 4;
        labelList{end+1}     = {'vis'};
        extraList{end+1}     = dir('/autofs/cluster/meso/users/divya/single_vessel_fMRI_long_short/svfmri_long_short_pilot03/mri/fmri_long/lbl/*.nii.gz');


        % subList{end+1} = 'svfmri_long_short_pilot06';
        % sesList{end+1} = '1';
        % sesDbList{end+1} = 'svfmri_long_short_pilot06/mri/fmri_long_ushoot';
        %
        % subList{end+1} = 'svfmri_long_short_pilot04';
        % sesList{end+1} = '1';
        % sesDbList{end+1} = 'svfmri_long_short_pilot04/mri/fmri_long_ushoot';
        %
        % subList{end+1} = 'svfmri_long_short_pilot03';
        % sesList{end+1} = '1';
        % sesDbList{end+1} = 'svfmri_long_short_pilot03/mri/fmri_long_ushoot';
end

%%%%%%%%%%%%%%%
%% Subject loop
volTsAll = cell(size(subList));
volPsdAll = cell(size(subList));
allCmd = cell(size(subList));

subIndList = 1:length(subList);
for i = 1:length(subList(subIndList))
    I = subIndList(i);

    % see doIt_pilotPipeline08c/stepTemplate.m
    %% Housekeeping
    info.sub   = subList{I};
    info.ses   = sesList{I};
    info.sesDb = sesDbList{I};
    info.tr    = volTrList{I};
    info.dummy = dummyList{I};
    info.K     = kList{I};
    info.win   = [winLengthList{I} winStepList{I}];

    dsgn.dt     = dtStimList{I};
    dsgn.onsets = onsetsList{I};
    dsgn.ondurs = ondursList{I};
    dsgn.label  = labelList{I};
    onsets = onsetsList{I};
    durs = ondursList{I};

    tic; disp(' '); disp(' '); disp(' ');
    tmp = {['sub-' info.sub '_ses-' info.ses]};
    tmp{end+1} = [num2str(i) '/' num2str(length(subIndList))];
    tmp{end+1} = ['Nsub=' num2str(length(subList))];
    tmp = strjoin(tmp,'; ');
    disp(repmat('+',1,length(tmp))); disp(repmat('+',1,length(tmp))); disp(repmat('+',1,length(tmp))); disp(tmp)



    %% Main analysis steps
    %%% run or load preprocessing
    do.loadIt = 0;
    do.doIt = 1;
    do.saveIt = 0;
    
    switch info.dataSetLabel
        
        case 'satin'
            switch info.sub
                case 'pilot01'
                    % info.funcGroup.id
                    info.funcGroup.id        = 'f0b0Nf1b1';
                    info.funcGroup.condLabel = info.runCond;
                    info.funcGroup.list      = {};
                    info.funcGroup.condInd   = [];
                    % volTsFileList = cell()
                    for runCond = 1:length(info.runCond)
                        tmpdir = dir(fullfile(info.dbDir,info.funcGroup.condLabel{runCond},'*_echo-1*run-*_*'));
                        tmpInd = ~contains({tmpdir.name}','-cat');
                        tmpFile = fullfile({tmpdir(tmpInd).folder},{tmpdir(tmpInd).name},info.dbFile)';


                        info.funcGroup.condInd = cat(1,info.funcGroup.condInd,ones(size(tmpFile)).*runCond);
                        info.funcGroup.list    = cat(1,info.funcGroup.list,fullfile({tmpdir(tmpInd).folder},{tmpdir(tmpInd).name},info.dbFile)');
                        % volTsFileList = fullfile({tmpdir(tmpInd).folder},{tmpdir(tmpInd).name},info.dbFile)';
                        % % bidsStrList = bidsStrList(tmpInd);
                        % for runInd = 1:size(volTsFileList,1)
                        %     volTsFileList{runInd}
                        % end
                    end
                    for runInd = 1:size(info.funcGroup.list,1)
                        disp(['reading run ' num2str(runInd) '/' num2str(size(info.funcGroup.list,1))])
                        curVolTs = MRIread(info.funcGroup.list{runInd});
                        curVolTs = vol2vec(curVolTs,[],-1);
                        % discard initial transient
                        curVolTs.vol(:,:,:,1:info.dummy) = [];
                        curVolTs.t(1:info.dummy,:) = [];
                        curVolTs.nframes = curVolTs.nframes - info.dummy;
                        curVolTs.dummy = info.dummy;
                        % add functional design
                        curVolTs.dsgn = dsgn;
                        curVolTs.runCondLabel = info.funcGroup.condLabel{info.funcGroup.condInd(runInd)};
                        volTs(runInd,1) = curVolTs;
                    end
            end



        case 'svfmri_long_short'
            info.sub = strsplit(info.sub,'_'); info.sub = info.sub{end};
            curBidsDir = fullfile(info.bidsDir,['sub-' info.sub],['ses-' info.ses]);
            if ~exist(curBidsDir,'dir'); mkdir(fullfile(curBidsDir,'func')); end


            curBidsFile = ['sub-' info.sub '_ses-' info.ses '_run-avg_angio'];

            %%%%%%%%%%%%%%
            %%%%%%%%%%%%%%
            %
            %%%%%%%%%%%%%%
            %%%%%%%%%%%%%%



            % copyfile(fullfile(info.dbDir,info.sesDb,info.dbFile),fullfile(curBidsDir,'func',[curBidsFile '.nii.gz']))

            info.funcGroup.id = 'runAvg';
            info.funcGroup.list = {curBidsFile};

            volTs = MRIread(fullfile(curBidsDir,'func',[info.funcGroup.list{1} '.nii.gz']));
            % correct tr if needed
            if isfield(info,'tr') && ~isempty(info.tr) && info.tr~=volTs.tr/1000
                volTs.tr = info.tr*1000;
                cmd = {srcFs};
                cmd{end+1} = ['mri_convert -tr ' num2str(info.tr*1000) ' ' volTs.fspec ' ' volTs.fspec];
                [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
            end
            volTs = vol2vec(volTs,[],-1);
            % discard initial transient
            volTs.vol(:,:,:,1:info.dummy) = [];
            volTs.t(1:info.dummy,:) = [];
            volTs.nframes = volTs.nframes - info.dummy;
            volTs.dummy = info.dummy;
            % add functional design
            volTs.dsgn = dsgn;
    end

    %%% run or load anatomical processing (masks and rois)
    do.loadIt = 0;
    do.doIt   = 1;
    do.saveIt = 0;
    % info.brainFile = volTs.mask.brain.fspec;
    volAnat = volAnatPreproc(do,info,extraList{I});


    % % % % % % % %%% simulate psd
    % % % % % % % % % % % onsets2 = onsets*2;
    % % % % % % % % % % % onsets2(end) = [];
    % % % % % % % % % % % durs2 = durs;
    % % % % % % % % % % % durs2(end) = [];
    % % % % % % % % % [volPsdSim, volTsSim] = simPsd(onsets2,durs2,volTs,info,0.01,0);
    % % % % % % % % % [volPsdSim, volTsSim] = simPsd(onsets,durs,volTs,info,0.01,0);
    % % % % % % % % % [volPsdSim, volTsSim] = simPsd2(onsets2,durs2,volTs,info,0.01,0);
    % % % % % % % % [volPsdSim, volTsSim] = simPsd2(onsets,durs,volTs,info,0.01,0);
    % % % % % % %
    % % % % % % % info3 = rmfield(info,{'onsetList' 'ondurList' 'funcGroup'});
    % % % % % % % info3.tr = 1;
    % % % % % % % info3.K  = 3;
    % % % % % % % dsgn3.dt = info3.tr;
    % % % % % % % dsgn3.onsetFreq = 1/100;
    % % % % % % % dsgn3.onFreq = 1/5;
    % % % % % % % dsgn3.dummy = 0;
    % % % % % % % dsgn3.nframes = round(15*60/dsgn3.dt);
    % % % % % % % dsgn3.winLength = round(50/dsgn3.dt)*dsgn3.dt; % [sec]
    % % % % % % % [volPsdSim, volTsSim] = simPsd3(dsgn3,info3,0);
    % % % % % % % midResp = volPsdSim.dsgn.onsets(1) + volPsdSim.dsgn.ondurs(1)*2/3;
    % % % % % % %
    % % % % % % % F = figure('WindowStyle','docked');
    % % % % % % % Ht = tiledlayout(8,1); Ht.TileSpacing = 'tight'; Ht.Padding = 'tight';
    % % % % % % % ax = {};
    % % % % % % % ax{end+1} = plotTs2(volTsSim,Ht);
    % % % % % % % ax{end+1} = plotPsd2(volPsdSim,Ht,midResp);
    % % % % % % % ax{end+1} = plotPsdGram2(volPsdSim,Ht);
    % % % % % % % ax{end+1} = plotPsdTrialGram2(volPsdSim,Ht,'av');
    % % % % % % % ax{end+1} = plotPsdTrialGram2(volPsdSim,Ht,'pc');
    % % % % % % % nexttile
    % % % % % % % ax{end+1} = plotPsdTrialGram2(volPsdSim,Ht,'MD');
    % % % % % % % ax{end+1} = plotPsdTrialGram2(volPsdSim,Ht,'MDpc');
    % % % % % % % xline(ax{end},midResp)
    % % % % % % %
    % % % % % % %
    % % % % % % %
    % % % % % % %
    % % % % % % % figure('WindowStyle','docked');
    % % % % % % % t   = volPsdSim.psdTrialGramMD.t(1,1,1,1,1,1,:,1);
    % % % % % % % f   = volPsdSim.psdTrialGramMD.f(1,1,1,1,:,1,1,1);
    % % % % % % % [~,wInd] = min(abs(t-66));
    % % % % % % % psd = volPsdSim.psdTrialGramMD.vec.psdPC(:,:,:,:,:,1,wInd);
    % % % % % % % plot(squeeze(f),squeeze(psd))
    % % % % % % %
    % % % % % % %
    % % % % % % % F = figure('WindowStyle','docked');
    % % % % % % % Ht = tiledlayout(8,1); Ht.TileSpacing = 'tight'; Ht.Padding = 'tight';
    % % % % % % % ax = {};
    % % % % % % % ax{end+1} = plotTs2(volTsSim,Ht);
    % % % % % % % ax{end+1} = plotCoh2(volPsdSim,Ht);
    % % % % % % % ax{end+1} = plotCohGram2(volPsdSim,Ht);
    % % % % % % % ax{end+1} = plotCohTrialGram2(volPsdSim,Ht,'av');
    % % % % % % % ax{end+1} = plotCohTrialGram2(volPsdSim,Ht,'pc');
    % % % % % % % ax{end+1} = plotCohTrialGram2(volPsdSim,Ht,'ek');
    % % % % % % % ax{end+1} = plotCohTrialGram2(volPsdSim,Ht,'MD');
    % % % % % % % ax{end+1} = plotCohTrialGram2(volPsdSim,Ht,'MDpc');



    %%% run or load full stimulus-response
    do.loadIt  = 0;
    do.doIt    = 1;
    do.saveIt  = 0;
    do.writeIt = 1;
    % info.onsetList = onsets;
    % info.ondurList = durs;
    for runInd = 1:size(volTs,1)
        resp(runInd) = volTsGetResp(do,info,volTs(runInd),volAnat(runInd));
    end
    [volTs.resp] = deal(resp(:)); clear resp
    % viewResp(volTs.resp,volAnat);



    %%% run or load full (spatio-temporal-spectral) multitaper analysis
    do.loadIt  = 0;
    do.doIt    = 1;
    do.saveIt  = 0;
    do.writeIt = 1;
    % volPsd = volPsdFullMtAna(do,info,volTs,volAnat);
    info.onsetList = onsets;
    info.ondurList = durs;
    info.skipSvd = 0;
    info.win = inf;
    info.perm = 0;


    if ~isfield(info.funcGroup,'condLabel')
        volPsd = volPsdFullMt2(do,info,volTs,volAnat);
        allCmd{I} = strjoin(viewResp(volTs.resp,volAnat,{volPsd.psd.fspec volPsd.svd.fspec.spSVmag},{':colormap=turbo' ':colormap=turbo'}),newline);
        % viewResp(volTs.resp,volAnat,{volPsd.psd.fspec volPsd.svd.fspec.spSVmag volPsd.svd.fspec.spSVmag_fdrThresh volPsd.svd.fspec.spSVphase_fdrThresh},{':colormap=turbo' ':colormap=turbo' ':colormap=turbo' ':colormap=turbo'});
        tryL2svd(volPsd)

        % I=3;
        % clipboard('copy',allCmd{I})
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




    %% Housekeeping
    disp(' '); tmp = ['sub-' info.sub '_ses-' info.ses ': DONE']; toc;
    disp(tmp); disp(repmat('+',1,length(tmp))); disp(repmat('+',1,length(tmp))); disp(repmat('+',1,length(tmp))); disp(' ');

    return

    if 1
        %% Plot
        x = 180;%75;%289;%146;%210;%190;%220;%180;
        y = 221;%135;%195;%263;%251;%179;%178;%220;
        % f0 = 0.11;
        f0 = 1/sum([40.8000   16.0000]);


        volTsTmp = vec2vol(vol2vec(volTs));
        % figure('WindowStyle','docked');
        % ht = tiledlayout(2,1);
        % dt = volTsTmp.tr/1000;
        % t = -1:dt:((size(volTsTmp.vol,4)-1)*dt-1);
        % nexttile
        % imagesc(volTsTmp.imMean)
        % hold on
        % xline(x,'r'); yline(y,'r');
        % ax = gca; ax.DataAspectRatio = [1 1 1]; ax.Colormap = gray; ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off';
        % nexttile
        % plot(t,squeeze(volTsTmp.vol(y,x,:,:)),'k')
        % grid on
        % axis tight
        % xlabel('time (s)')
        % xline([0 16],'r')
        % ylabel('demeaned response')
        % title(ht,subList{1},'interpreter','none')
        %
        %

        figure('WindowStyle','docked');
        ht = tiledlayout(5,4); ht.Padding = 'tight'; ht.TileSpacing = "tight";

        %%% background
        axIm = {};
        axIm{end+1} = nexttile([2 1]);
        imagesc(volTsTmp.imMean)
        hold on
        xline(x,'r'); yline(y,'r');
        ax = gca; ax.DataAspectRatio = [1 1 1]; ax.Colormap = gray; ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off';
        ylabel(colorbar,'signal (a.u.)')



        %%% power distribution
        axIm{end+1} = nexttile([2 1]);
        volPsd = vec2vol(volPsd);
        [~,b] = min(abs(volPsd.f-f0));
        imagesc(volPsd.vol(:,:,1,b))
        ax = gca; ax.DataAspectRatio = [1 1 1]; ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off';
        ylabel(colorbar,'psd')
        % ax.ColorScale = 'log';
        ax.Colormap = jet;
        hold on
        xline(x,'r'); yline(y,'r');
        ax = gca; ax.DataAspectRatio = [1 1 1]; ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off';



        %%% spatial component mag
        axIm{end+1} = nexttile([2 1]);
        [~,b] = min(abs(volPsd.f-f0));
        tmp = zeros(size(volPsd.vol2vec));
        tmp(volPsd.vol2vec) = volPsd.svd.u(:,1,b);
        imagesc(abs(tmp))
        ax = gca; ax.DataAspectRatio = [1 1 1]; ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off';
        ylabel(colorbar,'singular vector mag')
        % ax.ColorScale = 'log';
        ax.Colormap = jet;
        hold on
        xline(x,'r'); yline(y,'r');
        ax = gca; ax.DataAspectRatio = [1 1 1]; ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off';

        %%% spatial component phase
        axIm{end+1} = nexttile([2 1]);
        imagesc(angle(tmp))
        ax = gca; ax.DataAspectRatio = [1 1 1]; ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off';
        ylabel(colorbar,'singular vector phase')
        % ax.ColorScale = 'log';
        ax.Colormap = hsv;
        hold on
        xline(x,'r'); yline(y,'r');
        ax = gca; ax.DataAspectRatio = [1 1 1]; ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off';



        linkaxes([axIm{:}])

        K = volPsd.K;
        T = volPsd.nframes*volTs.tr/1000;
        [TW,W,K] = K2W(T,K);

        %%% single-vox spectrum
        axSingleVox = nexttile([1 4]);
        plot(squeeze(volPsd.f),squeeze(volPsd.vol(y,x,1,:)),'k')
        grid on
        axis tight
        xlabel('f (Hz)')
        ylabel('power')
        ax = gca; ax.YScale = 'log';
        title('single-vox')

        xline(f0,'r')
        yLim = ylim; hold on
        plot(f0+[-1 1]*W,yLim(2)*[1 1],'r','LineWidth',3)

        f0x = [1/16.0000 1/sum([40.8000   16.0000])];
        hFund = xline(f0x,'g');
        % f0x = 1/16.0000 + [-1 1].*1/sum([40.8000   16.0000]);
        % hMod = xline(f0x,'m');



        %%% average spectrum
        axAv = nexttile([1 4]);
        volPsd = vol2vec(volPsd);
        plot(squeeze(volPsd.f),squeeze(mean(volPsd.vec,2)),'k')
        grid on
        axis tight
        xlabel('f (Hz)')
        ylabel('psd')
        ax = gca; ax.YScale = 'log';
        title('cross-voxel average')
        yLim = ylim;
        yLim(1) = 23;
        ylim(yLim)
        xlim(volPsd.f([1 end]))

        xline(f0,'r')
        yLim = ylim; hold on
        plot(f0+[-1 1]*W,yLim(2)*[1 1],'r','LineWidth',3)

        f0x = [1/16.0000 1/sum([40.8000   16.0000])];
        hFund = xline(f0x,'g');
        % f0x = 1/16.0000 + [-1 1].*1/sum([40.8000   16.0000]);
        % hMod = xline(f0x,'m');




        %%% Coherence
        nexttile([1 4])
        f = squeeze(volPsd.f);
        coh = squeeze(volPsd.svd.coh(1,1,:));
        plot(f,coh,'k')
        grid on
        axis tight
        xlabel('f (Hz)')
        ylabel('coherence')
        ylim([1/volPsd.K 1])

        title(ht,subList{1},'interpreter','none')

        xline(f0,'r')
        yLim = ylim; hold on
        plot(f0+[-1 1]*W,yLim(2)*[1 1],'r','LineWidth',3)
        f0x = [1/16.0000 1/sum([40.8000   16.0000])];
        hFund = xline(f0x,'g');
        % f0x = 1/16.0000 + [-1 1].*1/sum([40.8000   16.0000]);
        % hMod = xline(f0x,'m');

        legend([hFund(1)],{'stim fundamentals'})


        % [40.8000   16.0000]






        %% Simulate stimulus response
        if ~isempty(onsetsList)
            volPsdSim = simPsd(onsetsList{I},ondursList{I},volTs,info,0.001);

            axes(axSingleVox);
            yyaxis right
            plot(squeeze(volPsdSim.f),volPsdSim.vec)
            axSingleVox.YScale = 'log';

            axes(axAv);
            yyaxis right
            plot(squeeze(volPsdSim.f),volPsdSim.vec)
            axAv.YScale = 'log';

        end



        %% %%%%%%%%%%%%%%%%%%%%%%
        % Arbitrary simulations %
        %%%%%
        close all
        figure('WindowStyle','docked');
        dursListList = [];

        for i = 1:1:10
            TR = 3;
            fStim = 1/(TR*20);
            fDur  = 1/(TR*i);
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
        %%%%%

        %%%%%
        close all
        for i = 1:1:60
            TR = 3; %0.01;
            fStim = 1/(60);
            fDur  = 1/(1*i);
            nStim = 1;
            T = 1/fStim * (nStim+2);
            onsetsList = (1/fStim:1/fStim:1/fStim*(nStim+1))';
            onsetsList(end) = [];
            ondursList = ones(size(onsetsList)) .* 1/fDur;

            startOffset = 1;
            onsetsList = onsetsList-startOffset;
            % T = T - 60 - startOffset;

            volTsTmp = vec2vol(volTs);
            volTsTmp.vol = [];
            volTsTmp.tr = TR*1000;
            volTsTmp.nframes = T/TR;
            info.dtrndOrder = -1;
            [volPsdSim, volTsSim] = simPsd(onsetsList,ondursList,volTsTmp,info);


            figure('WindowStyle','docked');

            ht = tiledlayout(2,1);

            ax = nexttile([1 1]);
            t = squeeze(volTsSim.t)';
            s = squeeze(volTsSim.vol);
            plot(t,s); xlabel('time (s)'); ylabel('MR signal')

            s = dtrnd2(s,t,[],info.dtrndOrder);
            hold on
            plot(t,s);
            ylim([-0.8 1.2])


            ax = nexttile([1 1]);
            f   = squeeze(volPsdSim.f);
            psd = squeeze(volPsdSim.vec);
            plot(f,psd); xlabel('Hz'); ylabel('psd')
            xlim([0 0.5])
            ylim([1e-5 1e1])
            ax.YScale = 'log';
            % hold on
            % n = 2^nextpow2(size(s,1));
            % psd2 = fftshift(fft(s-mean(s),n));
            % psd2 = psd2(n/2:end);
            % plot(f,abs(psd2))
            grid on

            title(ht,num2str(1/fDur))



            % % plot3(squeeze(volTsSim.t)',squeeze(volTsSim.vol),ones(size(squeeze(volTsSim.vol))).*dursList(1)); hold on
            % % plot3(squeeze(volPsdSim.f),volPsdSim.vec,ones(size(volPsdSim.vec)).*dursList(1));
            % plot(squeeze(volPsdSim.f),volPsdSim.vec);
            % % hold on
            % xlim([0 1])
            % ylim([1e-5 1e1])
            %
            % dursListList(end+1) = dursList(1);
            %
            % ax = gca;
            % ax.YScale = 'log';
            %
            % title(num2str(dursList(1)))
            %
            % keyboard
        end
        %%%%%


        %%%%%
        close all
        for i = 0:0.1:3
            TR = 3; %0.01;
            fStim = 1/(60);
            fDur  = 1/10;
            nStim = 1;
            T = 1/fStim * (nStim+2);
            onsetsList = (1/fStim:1/fStim:1/fStim*(nStim+1))';
            onsetsList(end) = [];
            ondursList = ones(size(onsetsList)) .* 1/fDur;

            startOffset = i;
            onsetsList = onsetsList-startOffset;
            % T = T - 60 - startOffset;

            volTsTmp = vec2vol(volTs);
            volTsTmp.vol = [];
            volTsTmp.tr = TR*1000;
            volTsTmp.nframes = T/TR;
            info.dtrndOrder = -1;
            [volPsdSim, volTsSim] = simPsd(onsetsList,ondursList,volTsTmp,info,0.01);


            figure('WindowStyle','docked');

            ht = tiledlayout(2,1);

            ax = nexttile([1 1]);
            t = squeeze(volTsSim.t)';
            s = squeeze(volTsSim.vol);
            plot(t,s); xlabel('time (s)'); ylabel('MR signal')

            s = dtrnd2(s,t,[],info.dtrndOrder);
            hold on
            plot(t,s);
            ylim([-0.8 1.2])


            ax = nexttile([1 1]);
            f   = squeeze(volPsdSim.f);
            psd = squeeze(volPsdSim.vec);
            plot(f,psd); xlabel('Hz'); ylabel('psd')
            xlim([0 0.5])
            ylim([1e-5 1e1])
            ax.YScale = 'log';
            % hold on
            % n = 2^nextpow2(size(s,1));
            % psd2 = fftshift(fft(s-mean(s),n));
            % psd2 = psd2(n/2:end);
            % plot(f,abs(psd2))
            grid on

            title(ht,num2str(1/fDur))
        end
        %%%%%




        %%%%%
        close all
        % info.K = 4;
        % info.win = [inf 0];

        offset = 10;
        TR = 3;
        info.win = [40 TR];
        fStim = 1/80;
        fDur  = 1/12;
        nStim = 4;
        T = 1/fStim * (nStim) + offset;
        onsetsList = (offset:1/fStim:(offset+1/fStim*(nStim)))';
        onsetsList(end) = [];
        ondursList = ones(size(onsetsList)) .* 1/fDur;

        volTsTmp = vec2vol(volTs);
        volTsTmp.vol = [];
        volTsTmp.tr = TR*1000;
        volTsTmp.nframes = T/TR;
        info.dtrndOrder = 2;
        [volPsdSim, volTsSim] = simPsd(onsetsList,ondursList,volTsTmp,info,0.01);



        figure('WindowStyle','docked');

        ht = tiledlayout(3,1);

        ax = nexttile([1 1]);
        t = squeeze(volTsSim.t)';
        s = squeeze(volTsSim.vol);
        plot(t,s); xlabel('time (s)'); ylabel('MR signal')

        s = dtrnd2(s,t,[],info.dtrndOrder);
        hold on
        plot(t,s);
        ylim([-0.8 1.2])
        xlim([0 330])

        xline(onsetsList)
        % axis tight
        grid on


        ax = nexttile([1 1]);
        f   = squeeze(volPsdSim.f);
        psd = squeeze(volPsdSim.vec);
        plot(f,psd); xlabel('Hz'); ylabel('psd')
        xlim([0 0.5])
        % ylim([1e-5 1e1])
        ax.YScale = 'log';
        % hold on
        % n = 2^nextpow2(size(s,1));
        % psd2 = fftshift(fft(s-mean(s),n));
        % psd2 = psd2(n/2:end);
        % plot(f,abs(psd2))
        grid on

        title(ht,num2str(1/fDur))

        xline(fStim)



        ax = nexttile([1 1]);
        t   = squeeze(volPsdSim.psdGram.tWin)';
        f   = squeeze(volPsdSim.psdGram.f);
        psd = squeeze(volPsdSim.psdGram.vec);
        imagesc(t,f,psd)
        ax.ColorScale = 'log';
        grid on


        yline(fStim)
        ylim([0 0.5])
        xlim([0 T])




        %%%%%%



        % % plot3(squeeze(volTsSim.t)',squeeze(volTsSim.vol),ones(size(squeeze(volTsSim.vol))).*dursList(1)); hold on
        % % plot3(squeeze(volPsdSim.f),volPsdSim.vec,ones(size(volPsdSim.vec)).*dursList(1));
        % plot(squeeze(volPsdSim.f),volPsdSim.vec);
        % % hold on
        % xlim([0 1])
        % ylim([1e-5 1e1])
        %
        % dursListList(end+1) = dursList(1);
        %
        % ax = gca;
        % ax.YScale = 'log';
        %
        % title(num2str(dursList(1)))
        %
        % keyboard
    end
    %%%%%




    %%%%%
    %% %%%




    % tmp = splitOverlapping(volPsd);
    % plotAllGram(tmp)
    % plotAll3(volPsd,volTs,[],[],[],pipId);
    % tf = pickTimeFreq;
    %
    % [~,~,hFig1,hFig2] = plotAll2(volPsd,volTs,tf,vox,cropMask,pipId);
    %
    % tag = 'vlf4';
    % [tf,vox] = plotAll(volPsd,volTs,[],[],cropMask,pipId);
    % save(fullfile(outDir,[strjoin({['sub-' info.sub] ['ses-' info.ses] ['tfVox-' tag]},'_')]),'tf','vox','I')
    %
    % plotAll(volPsd,volTs,tf,vox,cropMask,pipId,1);

end

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
stepLabel = [stepLabel '; Nsubj=' num2str(length(subIndList)) '/' num2str(length(subList))];
stepFile = fullfile(outDir,['group-subj_' outVar]);
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

saveIt = 0;
outVar = 'volAnat';
stepLabel = 'anat processing';
% WARNING: outVar = outVarAll; clear outVarAll
if exist([outVar 'All'],'var'); eval([outVar ' = ' outVar 'All; clear ' outVar 'All']); end; warning([newline outVar ' = ' outVar 'All; clear ' outVar 'All']);
stepLabel = [stepLabel '; Nsubj=' num2str(length(subIndList)) '/' num2str(length(subList))];
stepFile = fullfile(outDir,['group-subj_' outVar]);
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

saveIt = 0;
outVar = 'volPsd';
stepLabel = 'full mt analysis';
% WARNING: outVar = outVarAll; clear outVarAll
if exist([outVar 'All'],'var'); eval([outVar ' = ' outVar 'All; clear ' outVar 'All']); end; warning([newline outVar ' = ' outVar 'All; clear ' outVar 'All']);
stepLabel = [stepLabel '; Nsubj=' num2str(length(subIndList)) '/' num2str(length(subList))];
stepFile = fullfile(outDir,['group-subj_' outVar]);
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

stepLabel = 'all data loaded';
disp(' '); disp(' '); disp(repmat('x',1,length(stepLabel))); disp(upper(stepLabel)); toc; disp(repmat('x',1,length(stepLabel))); disp(' '); disp(' ');



return
%%%%%%%%%%%%%%%%
%% Visualization
%% Plot All
subList
subIndList
I = 3;
cropMask = volAnat{I}.fun.mask.head;

% IAll = {}; tfAll = {}; voxAll = {};

% volTs{I}.t(:,:,:,:,:,2)
[tf,vox] = plotAll(volPsd{I},volTs{I},[],[],cropMask,pipId);
IAll{end+1} = I; tfAll{end+1} = tf; voxAll{end+1} = vox;

i = 1; %length(IAll);
I = IAll{i}; tf = tfAll{i}; vox = voxAll{i};
plotAll(volPsd{I},volTs{I},tf,vox,cropMask,pipId);
