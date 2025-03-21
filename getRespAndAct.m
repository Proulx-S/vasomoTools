function [fRespCat,fRespRun,fActCat,fActRun] = getRespAndAct(volTs,dsgn,fMask,param,force,verbose)
if ~exist('force','var');     force = []; end
if ~exist('verbose','var'); verbose = []; end
if ~exist('dsgn','var');       dsgn = []; end
if ~exist('fMask','var');     fMask = []; end
if isempty(force);     force = 0; end
if isempty(verbose); verbose = 0; end
if isempty(dsgn)
    if isfield(volTs,'dsgn') && isequal(dsgn)
        dsgn = dsgn;
    else
        error('badly specified dsgn')
    end
end
% if ~isfield(param,'skipMov'); param.skipMov = []; end
% if ~isfield(param,'skipCat'); param.skipCat = []; end
% if ~isfield(param,'skipRun'); param.skipRun = []; end
if ~isfield(param,'dryRun'); param.dryRun = []; end
% if isempty(param.skipMov); param.skipMov = 0; end
% if isempty(param.skipCat); param.skipCat = 0; end
% if isempty(param.skipRun); param.skipRun = 0; end
if isempty(param.dryRun); param.dryRun = 0; end

%% Data files
if ~isempty(volTs)
    if isfield(volTs,'fspec')
        fVolTs = {volTs.fspec}';
    else
        fVolTs = cell(size(volTs));
        for r = 1:length(volTs)
            fVolTs{r,1} = volTs(r).mri.fspec;
        end
    end
    for i = 1:length(fVolTs)
        if ~exist(fVolTs{i},'file'); dbstack; error('write file to disk aka code lazy bum'); end
    end
else
    fVolTs = [];
end


%% Mask files
if isempty(fMask)
    dbstack; error('double-check that')
    if isfield(volTs,'vol2vec')
        fMask = volTs.vol2vec;
    end
end
if ~all(diff([volTs.tr])<0.01); dbstack; error('runs have different tr'); end

%% Functional design
if ~isa(dsgn,'runDsgn')
    dbstack; error('old version, reconciliate')
end
if ~isempty(dsgn.cond)
    k = sort(unique(dsgn.cond)); if any(diff(k)-1); dbstack; error('cond indices are not monotonically increasing'); end
    dsgn.condLabel = repmat({'stim'},size(k));
    if any(k==0)
        dsgn.condLabel{k==0} = 'catch';
    end
    dsgn.condK = length(k);
else
    dbstack; error('double-check')
    param.funDsgn.condSeq   = ones(size(param.funDsgn.startSeq));
    param.funDsgn.condLabel = {'stim'};
    if isfield(dsgn,'nullTrial') && ~isempty(dsgn.nullTrial)
        param.funDsgn.condSeq(dsgn.nullTrial) = 2;
        param.funDsgn.condLabel{end+1} = 'catch';
        param.funDsgn.k = 2;
    end
end
param.tr = [volTs.tr]'./1000;
param.dsgn = dsgn;


%% Run afni's 3dDeconvolve for response timecourse estimation
param.model = 'TENTzero';
% On each run
R = size(fVolTs,1);
clear fRespRun
for r = 1:R
    fRespRun(r,:) = runAfni(fVolTs(r,:),[r R],param,fMask,force,verbose); % analysis performed on each echoe within that function
end
% On catenated runs
if R>1
    fRespCat      = runAfni(fVolTs,     [0 R],param,fMask,force,verbose); % analysis performed on each echoe within that function
end

%% Run afni's 3dDeconvolve for double-gamma response amplitude (and delay)
param.model = 'SPMG2';
% On each run
R = size(fVolTs,1);
clear fActRun
for r = 1:R
    fActRun(r,:) = runAfni(fVolTs(r,:),[r R],param,fMask,force,verbose); % analysis performed on each echoe within that function
end
% On catenated runs
if R>1
    fActCat      = runAfni(fVolTs,     [0 R],param,fMask,force,verbose); % analysis performed on each echoe within that function
end




%% Plot design matrices
verboseThis = verbose;
forceThis   = force;
for r = 1:R
    fRespRun(r,1).xMat = plotDsgnMat(fRespRun(r,1),forceThis,verboseThis);
end
if R>1
    fRespCat.xMat      = plotDsgnMat(fRespCat     ,forceThis,verboseThis);
end
for r = 1:R
    fActRun(r,1).xMat = plotDsgnMat(fActRun(r,1),forceThis,verboseThis);
end
if R>1
    fActCat.xMat      = plotDsgnMat(fActCat     ,forceThis,verboseThis);
end



%% Refactor
verboseThis = verbose;
forceThis   = force;
fRespRun = unpackAfni(fRespRun,[],forceThis,verboseThis);
fActRun  = unpackAfni(fActRun, [],forceThis,verboseThis);
if R>1
    fRespCat = unpackAfni(fRespCat,[],forceThis,verboseThis);
    try
        fActCat  = unpackAfni(fActCat, [],forceThis,verboseThis);
    catch
        keyboard
    end
else
    fRespCat = fRespRun; fRespCat.r = 0;
    fActCat  = fActRun;  fActCat.r = 0;
end




return




%%%%%%%%%%%%%%
%% Make videos
%%%%%%%%%%%%%%
switch param.model
    case {'SPMG2' 'SPMG3'}
    case {'TENT' 'TENTzero'}

        %%%% Individual runs
        if ~param.skipRun
            for i = 1:numel(fRun)
                fIn = fRun(i).fs.fRespTsOnBaseAv;
                fOut = replace(fIn,'.nii.gz','');
                % fRun(i).mov.fRespOnBaseAvMovie = [fOut '.avi'];
                % fRun(i).mov.fRespOnBaseAvMovieHighBit = [fOut '.mj2'];
                fRun(i).mov.fRespOnBaseAvMovie = [];
                fRun(i).mov.fRespOnBaseAvMovieHighBit = [];
            end
        end
        %%%% Whole session
        if ~param.skipCat
            for i = 1:numel(fSes)
                fIn = fSes(i).fs.fRespTsOnBaseAv;
                fOut = replace(fIn,'.nii.gz','');
                fSes(i).mov.fRespOnBaseAvMovie = [fOut '.avi'];
                fSes(i).mov.fRespOnBaseAvMovieHighBit = [fOut '.mj2'];
            end
        end



        if ~param.skipMov && volTs(1).depth==1
            forceThis = force;
            disp('Making movies')
            nLoop = 4;
            if ~isempty(fMask)
                mask = MRIread(fMask); mask = logical(mask.vol);
            else
                mask = true(volTs.height,volTs.width,volTs.depth);
            end

            % %%%% Individual runs
            % for R = 1:nRun
            %     for E = 1:nEcho+1
            %         if E>nEcho && nEcho==1; break; end
            %         disp([' run' num2str(R) '/' num2str(nRun)])
            %         disp(['  file' num2str(E) '/' num2str(size(fSes,2)+1)])
            %         if E>nEcho
            %             fIn = fRun_echoRms(R).fRespOnBase;
            %             fOut = replace(fIn,'.nii.gz','');
            %             fRun_echoRms(R).fRespOnBaseMovie = [fOut '.avi'];
            %             fRun_echoRms(R).fRespOnBaseMovieHighBit = [fOut '.mj2'];
            %         else
            %             fIn = fRun(R,E).fRespOnBase;
            %             fOut = replace(fIn,'.nii.gz','');
            %             fRun(R,E).fRespOnBaseMovie = [fOut '.avi'];
            %             fRun(R,E).fRespOnBaseMovieHighBit = [fOut '.mj2'];
            %         end
            %         if forceThis || ~exist([fOut '.avi'],'file')
            %             vOut = VideoWriter(fOut,'Uncompressed AVI');
            %         end
            %         if forceThis || ~exist([fOut '.mj2'],'file')
            %             vOutHighBit = VideoWriter(fOut,'Archival');
            %         end
            %
            %         if forceThis || ~exist([fOut '.avi'],'file') || ~exist([fOut '.mj2'],'file')
            %             resp = MRIread(fIn);
            %             % Scale
            %             mask = repmat(mask,[1 1 resp.depth resp.nframes]);
            %             resp.vol(mask) = resp.vol(mask) - min(resp.vol(mask));
            %             resp.vol(mask) = resp.vol(mask) ./ max(resp.vol(mask));
            %             mask = mask(:,:,1,1);
            %             % Crop
            %             resp.vol(all(~mask,2),:,:,:) = [];
            %             resp.vol(:,all(~mask,1),:,:) = [];
            %             % Upsample
            %             resp.vol = imresize(resp.vol,3,'nearest');
            %             % Set frame rate to 1cycle/sec
            %             vOutHighBit.FrameRate = resp.nframes;
            %             vOut.FrameRate = resp.nframes;
            %
            %             % Write
            %             open(vOut)
            %             open(vOutHighBit)
            %             for L = 1:nLoop
            %                 for f = 1:resp.nframes
            %                     writeVideo(vOut,resp.vol(:,:,1,f));
            %                     writeVideo(vOutHighBit,uint16(resp.vol(:,:,1,f)*(2^16-1)));
            %                 end
            %             end
            %             close(vOut)
            %             close(vOutHighBit)
            %
            %             disp('  done')
            %         else
            %             disp('  already done,skipping')
            %         end
            %     end
            % end


            if ~param.skipCat
                %%%% Whole session
                disp(' runCat')
                for E = 1:nEcho+1
                    if E>nEcho && nEcho==1; break; end
                    disp(['  file' num2str(E) '/' num2str(size(fSes,2)+1)])
                    if E>nEcho && nEcho>1
                        dbstack; error('code that');
                        fIn = fSes_echoRms.fRespOnBase;
                        fOut = replace(fIn,'.nii.gz','');
                        fSes_echoRms.fRespOnBaseMovie = [fOut '.avi'];
                        fSes_echoRms.fRespOnBaseMovieHighBit = [fOut '.mj2'];
                    else
                        fIn = fSes(1,E).fs.fRespTsOnBaseAv;
                        fOut = replace(fIn,'.nii.gz','');
                        % fSes(1,E).mov.fRespOnBaseAvMovie = [fOut '.avi'];
                        % fSes(1,E).mov.fRespOnBaseAvMovieHighBit = [fOut '.mj2'];
                    end
                    if forceThis || ~exist([fOut '.avi'],'file')
                        vOut = VideoWriter(fOut,'Uncompressed AVI');
                    end
                    if forceThis || ~exist([fOut '.mj2'],'file')
                        vOutHighBit = VideoWriter(fOut,'Archival');
                    end

                    if forceThis || ~exist([fOut '.avi'],'file') || ~exist([fOut '.mj2'],'file')
                        resp = MRIread(fIn);
                        mask = repmat(mask,[1 1 resp.depth resp.nframes]);
                        resp.vol(~mask) = 0;
                        % Scale
                        resp.vol(mask) = resp.vol(mask) - min(resp.vol(mask));
                        resp.vol(mask) = resp.vol(mask) ./ max(resp.vol(mask));
                        mask = mask(:,:,1,1);
                        % Crop
                        resp.vol(all(~mask,2),:,:,:) = [];
                        resp.vol(:,all(~mask,1),:,:) = [];
                        % Upsample
                        resp.vol = imresize(resp.vol,3,'nearest');
                        % Set frame rate to 1cycle/sec
                        vOutHighBit.FrameRate = resp.nframes;
                        vOut.FrameRate = resp.nframes;

                        % Write
                        open(vOut)
                        open(vOutHighBit)
                        for L = 1:nLoop
                            for f = 1:resp.nframes
                                writeVideo(vOut,resp.vol(:,:,1,f));
                                writeVideo(vOutHighBit,uint16(resp.vol(:,:,1,f)*(2^16-1)));
                            end
                        end
                        close(vOut)
                        close(vOutHighBit)

                        disp('  done')
                    else
                        disp('  already done,skipping')
                    end
                end
            end
        end
    otherwise
        dbstack; error('X');
end





function fRes = runAfni(fList,rR,param,fMask,force,verbose)
    global src
    if ~exist('rR','var');                     rR = []; end
    if ~exist('fMask','var');               fMask = []; end
    if ~exist('force','var');               force = []; end
    if ~exist('verbose','var');           verbose = []; end
    if ~isfield(param,'getResid'); param.getResid = []; end
    if isempty(force);                   force = 0; end
    if isempty(verbose);               verbose = 0; end
    if isempty(param.getResid); param.getResid = 0; end

    r = rR(1); % r=0 for catenated runs
    R = rR(2); % R is the number of runs, whether catenated or not
    
    % Handle multiple runs
    if r
        param.tr            = param.tr(r,:);
        param.nFrame        = param.nFrame(r,:);
        param.nDummyRemoved = param.nDummyRemoved(r,:);
    end
    


    %%%%%%%%%%%%%%%%%%%%
    %% Functional design
    %%%%%%%%%%%%%%%%%%%%
    if isfield(param,'dsgn') && isa(param.dsgn,'runDsgn')
        trStim   = param.dsgn.dt;
        durSeq   = param.dsgn.ondurList;
        condSeq  = param.dsgn.cond;
        startSeq = param.dsgn.onsetList;
        HRmodel  = param.model;
    else
        dbstack; error('double-check')
        k = param.funDsgn.k;
        trStim = param.funDsgn.trStim;
        durSeq =  param.funDsgn.durSeq;
        condSeq = param.funDsgn.condSeq;
        startSeq = param.funDsgn.startSeq;
        HRmodel = param.model;
    end
    if isfield(param,'nDummy') && ~isempty(param.nDummy)
        dbstack; error('old convention, double-check'); %param.nDummyIgnore = param.nDummy;
    else
        param.nDummyIgnore = 0;
    end


    sz = size(fList);
    multiRun  = sz(1)>1;
    multiEcho = sz(2)>1;
    

    %%%%%%%
    %% Mask
    %%%%%%%
    mriMask = MRIread(fMask);
    mriMask.vol([1:5 end-4:end],:              ) = 0;
    mriMask.vol(:              ,[1:5 end-4:end]) = 0;




    cmd = {src.afni};
    if size(fList,2)>1; dbstack; error('double-check that'); end
    for E = 1:size(fList,2)

        %% Define files
        fIn = fList(:,E);
        if multiRun
            fOut = char(fIn(1)); fOut = strsplit(fOut,'_');
            fOut{contains(fOut,'run-')} = 'run-cat'; fOut = strjoin(fOut,'_'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        else
            fOut = char(fIn);
        end
        

        fStat  = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats']); % let's let afni use its native file format, it is sometimes glitchy otherwise
        fFit   = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_fit.nii.gz']);
        if param.getResid
            fResid = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_resid.nii.gz']);
        else
            fResid = '';
        end
        fMask  = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_mask.nii.gz' ]);
        mriMask.fspec = fMask; MRIwrite(mriMask,fMask);
        fResp    = repmat({''},size(param.dsgn.condLabel));
        fRespStd = repmat({''},size(param.dsgn.condLabel));
        switch HRmodel
            case {'TENT' 'TENTzero'}
                for k = 1:param.dsgn.condK
                    fResp{k}    = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_respAv.nii.gz']);
                    fRespStd{k} = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_respSd.nii.gz']);
                end
            case {'SPMG2' 'SPMG3'}
            otherwise
                dbstack; error('figure that out')
        end
        if param.dryRun
            dbstack; error('double-check that')
            tmpName = tempname;
            fStim   = [tmpName '_startTime.1D'  ];
            fMat    = [tmpName '_stats.xmat.1D' ];
            fMatFig = [tmpName '_stats.xmat.fig'];
        else
            % Make sure durations are the same across trials within the same condition
            kList = sort(unique(param.dsgn.cond));
            durList = [];
            for ik = 1:param.dsgn.condK
                durList = [durList diff(param.dsgn.ondurList(param.dsgn.cond==kList(ik)))];
            end
            if any(durList); dbstack; error('different durations across trials within the same condition, code that (hint: -stim_times_AM1)'); end
            % Move on assuming same duration across trials within the same condition
            fStim = cell(size(param.dsgn.condLabel));
            for k = 1:param.dsgn.condK
                fStim{k} = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_startTime.1D']);
            end
            fMat    = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats.xmat.1D']);
            fMatFig = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats.xmat.fig']);
        end

   

        fRes(1,E).fIn      = fIn;
        fRes(1,E).fFit     = fFit;
        fRes(1,E).fResid   = fResid;
        fRes(1,E).fStat    = fStat;
        fRes(1,E).fMask    = fMask;
        fRes(1,E).fMat     = fMat;
        fRes(1,E).fMatFig  = fMatFig;
        fRes(1,E).fResp    = fResp;
        fRes(1,E).fRespStd = fRespStd;
        fRes(1,E).r        = r;
        fRes(1,E).R        = R;
        


        %% Contruct afni command
        cmdTmp = {};
        if multiRun
            cmdTmp{end+1} = 'echo catenated runs';
        end
        if multiEcho
            cmdTmp{end+1} = ['echo ' num2str(E) '/' num2str(size(fList,2))];
        end

        [cmdTmpTmp,param.dsgn.nReg] = afniCmd(fIn,fStim,fMask,param,fResp,fRespStd,fFit,fResid,fMat,fStat,verbose,param.dryRun);
        
        if param.dryRun
            system(strjoin([{srcAfni} cmdTmpTmp],newline))
        end

        if force || anyDontExist([fResp fRespStd [fStat '+orig.BRIK']])
            cmdTmp = [cmdTmp cmdTmpTmp];
        end
        cmdTmp{end+1} = ['echo ''   ''' strjoin(cellstr(fResp)   ,' ')];
        cmdTmp{end+1} = ['echo ''   ''' strjoin(cellstr(fRespStd),' ')];
        cmdTmp{end+1} = ['echo ''   '''                [fStat '+orig']         ];
        cmdTmp{end+1} = ['echo ''   '''                 fMat          ];
        if ~force && ~anyDontExist([fResp fRespStd [fStat '+orig.BRIK']])
            cmdTmp{end+1} = 'echo ''   ''already done, skipping';
        end

        fRes(1,E).param = param;
        fRes(1,E).cmd =  strjoin(cmdTmp,newline);

        cmd = [cmd cmdTmp];
    end

    %% Run afni command
    [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end




function [cmd,nReg] = afniCmd(fIn,fStim,fMask,param,fResp,fRespStd,fFit,fResid,fMat,fStat,verbose,dryRun)
    % function [cmd,nReg] = afniCmd2(fIn,fMask,fStim,nDummyIgnore,tr,startSeq,durSeq,condSeq,HRmodel,label,param,fResp,fFit,fResid,fMat,fStat,verbose,nDummyRemoved,trDecon,dryRun,nFrame)
    % param.nDummyRemoved [int]: number of initial frames that are already removed from the
    % timeseries. The stimulus timeseries must therefore be adjusted
    % accordingly.
    % param.nDummyIgnore [int]: number of initial frames to ignore from the input
    % timeseries. Only frames after these will be feed to 3dDeconvolve using
    % the [nDummyIgnore..$] notation.
    % nDummy [int]: total number of dummy initial frames
    tr      = mean(param.tr);
    nFrame  = max(param.nFrame);
    trDecon = param.trDecon;
    if any(diff(param.nDummyRemoved)>1); dbstack; error('nDummyRemoved should be the same across runs'); end
    nDummy = mode(param.nDummyIgnore + param.nDummyRemoved);
    fIn = cellstr(fIn);
    cmd = {'3dDeconvolve -overwrite \'};
    if ~dryRun
        cmd{end+1} = ['-input ' sprintf(['%s[' num2str(param.nDummyIgnore) '..$] '],fIn{:}) ' \'];
        if ~isempty(fMask)
            cmd{end+1} = ['-mask ' fMask ' \'];
        end
    else
        dbstack; error('code that')
        if isempty(nFrame)
            nFrame = MRIread(fIn{1},1); nFrame = nFrame.nframes;
        end
        cmd{end+1} = ['-nodata ' num2str(nFrame) ' ' num2str(tr,'%0.16f') ' \'];
    end
    cmd{end+1} = '-polort A -local_times \';
    cmd{end+1} = ['-stim_times_subtract ' num2str(mean(tr.*nDummy),'%f') ' \'];
    
    % Set design
    dsgn = param.dsgn;
    nRegAll = [];
    cmd{end+1} = ['-num_stimts ' num2str(dsgn.condK) ' \'];
    kList = sort(unique(dsgn.cond));
    cmd{end+1} = ['-TR_times ' num2str(trDecon,'%f') ' \'];
    for k = 1:dsgn.condK
        cmd{end+1} = ['-stim_label ' num2str(k) ' ' [char(dsgn.task) '_' dsgn.condLabel{k}] ' \'];

        % write design to file
        if dryRun
            fStim = [tempname '_startTime.1D' ];
        end
        fido = fopen(fStim{k}, 'w');
        if ~iscell(fIn); dbstack; error('fIn must be type cell'); end
        for i = 1:size(fIn,1)
            fprintf(fido,'%.3f ',dsgn.onsetList(kList(k)==dsgn.cond));
            fprintf(fido,'\n');
        end
        fclose(fido);

        dur = dsgn.ondurList(kList(k)==dsgn.cond); if ~isempty(dur) && any(diff(dur)); dbstack; error('stim duration cannot be different across trials'); end
        dur = dur(1);
        
        % Set model
        switch param.model
            case 'SPMG2'
                nReg = 2;
                cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim{k} ' ''' param.model '(' num2str(dur,'%0.3f') ')'' \'];
                nRegAll(k) = nReg;
            case 'SPMG3'
                dbstack; error('double-check that')
                nReg = 3;
                if max(abs(diff(durSeq)))/max(durSeq) > 0.0001; dbstack; error('stim duration cannot be different across trials'); end
                cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' ''' HRmodel '(' num2str(mean(durSeq),'%0.3f') ')'' \'];
            case 'TENT'
                dbstack; error('code that')
            case 'TENTzero'
                % set the deconvolution window to the maximum (all the way up to the next stimulus or the end of the run)
                eTime     = dsgn.onsetList(dsgn.cond==kList(k));
                eTimeNext = find(dsgn.cond==kList(k))+1;
                if eTimeNext(end) > length(dsgn.onsetList)
                    eTimeNext(end) = [];
                    eTimeNext = dsgn.onsetList(eTimeNext);
                    eTimeNext(end+1) = (nFrame + mode(param.nDummyRemoved)) * tr;
                else
                    eTimeNext = dsgn.onsetList(eTimeNext);
                end
                deconWin = min(eTimeNext - eTime);
                if (deconWin/trDecon)/ceil(deconWin/trDecon)>0.9
                    deconWin = ceil(deconWin/trDecon)*trDecon;
                else
                    deconWin = floor(deconWin/trDecon)*trDecon;
                end
                b = 0;
                c = round((deconWin-trDecon)/trDecon)*trDecon;
                nReg = round( (c-b)/trDecon + 1 );
                % (c-b)/(nReg-1)
                cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim{k} ' ''TENTzero(' num2str(b) ',' num2str(c) ',' num2str(nReg) ')'' \'];
                nReg = nReg - 2;
                if ~dryRun
                    cmd{end+1} = ['-iresp ' num2str(k) ' ' fResp{k}    ' \'];
                    cmd{end+1} = ['-sresp ' num2str(k) ' ' fRespStd{k} ' \'];
                end
                nRegAll(k) = nReg;
            otherwise
                dbstak; error('X');
        end
    end
    nReg = nRegAll;


    if dryRun
        dbstack; error('code that')
        fMat = replace(fMat,'_stats.xmat.1D','_stats.xmatPerTrial.1D');
        fMat  = [tempname '_stats.xmat.1D'];
    end



    % Set outputs
    if ~dryRun
        cmd{end+1} = ['-fitts ' fFit ' \'];
        if ~isempty(fResid)
            cmd{end+1} = ['-errts ' fResid ' \'];
        end
        cmd{end+1} = '-bout -fout \';
    end
    cmd{end+1} = ['-x1D ' fMat ' \'];
    if ~dryRun
        if verbose>0
            cmd{end+1} = ['-bucket ' fStat];
        else
            cmd{end+1} = ['-bucket ' fStat ' 2>/dev/null'];
        end
    else
        cmd{end}(end-1:end) = [];
    end
