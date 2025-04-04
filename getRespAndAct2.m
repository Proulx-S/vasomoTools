function [fRespCat,fRespRun,fActCat,fActRun] = getRespAndAct2(fVolTs,dsgn,fMask,param,force,verbose,passDown)
if ~exist('force','var');     force   = []; end
if ~exist('verbose','var'); verbose   = []; end
if ~exist('dsgn','var');       dsgn   = []; end
if ~exist('fMask','var');     fMask   = []; end
if ~exist('passDown','var'); passDown = {}; end
if isempty(force);     force = 0; end
if isempty(verbose); verbose = 0; end
% if isempty(dsgn)
%     if isfield(fVolTs,'dsgn') && isequal(dsgn)
%         dsgn = dsgn;
%     else
%         error('badly specified dsgn')
%     end
% end
% if ~isfield(param,'skipMov'); param.skipMov = []; end
% if ~isfield(param,'skipCat'); param.skipCat = []; end
% if ~isfield(param,'skipRun'); param.skipRun = []; end
if ~isfield(param,'dryRun'); param.dryRun = []; end
if ~isfield(param,'PCflag'); param.PCflag = []; end
        % if isempty(param.skipMov); param.skipMov = 0; end
% if isempty(param.skipCat); param.skipCat = 0; end
% if isempty(param.skipRun); param.skipRun = 0; end
if isempty(param.dryRun); param.dryRun = 0; end
if isempty(param.PCflag); param.PCflag = 0; end

%% Assert data files
if ~iscell(fVolTs) && ~ischar(fVolTs{1})
    dbstack; error('old version, reconciliate')
end

%% Assert dsgn
if ~isa(dsgn,'runDsgn')
    dbstack; error('old version, reconciliate')
end

%% Assert mask files
if ~iscell(fMask) && ~ischar(fMask)
    dbstack; error('old version, reconciliate')
end

%% Assert param
param.dsgn = dsgn; clear dsgn
if any((param.nFrameOrig-param.nFrame).*param.tr > param.dsgn.onsetList(1))
    disp('!!!!!!!!!!!!!!!!!!!!!')
    disp('WARNING: too many dummy removed->timeseries begins after first event')
    disp('WARNING: discarding first event, but you really should reprocess with less dummies')
    disp('!!!!!!!!!!!!!!!!!!!!!')
    param.dsgn.onsetList(1) = [];
    param.dsgn.ondurList(1) = [];
    param.dsgn.cond(1)      = [];
    % dbstack; error('too many dummy removed->timeseries begins after first event');
end


%% Run afni's 3dDeconvolve for response timecourse estimation
param.model = 'TENTzero';
R = size(fVolTs,1);
% On each run
clear fRespRun
for r = 1:R
    if param.PCflag
        fRespRun(r,:) = runAfni(fVolTs(r,:,:),[r R],param,fMask(r,1),force,verbose); % analysis performed on each echoe within that function
    else
        fRespRun(r,:) = runAfni(fVolTs(r,1),[r R],param,fMask(r,1),force,verbose); % analysis performed on each echoe within that function
    end
end
% On catenated runs
if R>1
    fRespCat      = runAfni(fVolTs,     [0 R],param,fMask,force,verbose); % analysis performed on each echoe within that function
end




% %% Rerun special case
% if any(diff(param.nFrame))
%     forceThis = 1;
%     switch fVolTs(1).fOrigList{1}
%         case '/local/users/sebp/martinos/vsmDriven/doIt_generalPreproc/vsmDriven/bids/sub-vsmDrivenP5/ses-1/func/sub-vsmDrivenP5_ses-1_task-10sPrd1sDur_acq-vfMRIinflow_run-1_angio.nii.gz'
%             rBad = 3;
%             rRef = 1;
%             % Extract parameters to fix
%             cmd = strsplit(fRespRun(rRef,:).cmd,newline)';
%             cmd = strsplit(cmd{contains(cmd,'TENTzero(')});
%             passDown{end+1}.pr = cmd{contains(cmd,'TENTzero(')};
%             passDown{end  }.id = 'TENTzeroParam';
%             % Rerun, passing down the fixed parameters
%             fRespRun(rBad,:) = runAfni(fVolTs(rBad,:),[rBad R],param,fMask,forceThis,verbose,passDown); % analysis performed on each echoe within that function
%         otherwise
%             error('need to define parameters for this special case')
%     end
% end


%% Run afni's 3dDeconvolve for double-gamma response amplitude (and delay)
param.model = 'SPMG2';
% On each run
R = size(fVolTs,1);
clear fActRun
for r = 1:R
    if param.PCflag
        fActRun(r,:) = struct;
    else
        fActRun(r,:) = runAfni(fVolTs(r,:),[r R],param,fMask,force,verbose); % analysis performed on each echoe within that function
    end
end
% On catenated runs
if R>1
    if param.PCflag
        fActCat      = struct;
    else
        fActCat      = runAfni(fVolTs,     [0 R],param,fMask,force,verbose); % analysis performed on each echoe within that function
    end
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
    if numel(fields(fActRun(r,1)))
        fActRun(r,1).xMat = plotDsgnMat(fActRun(r,1),forceThis,verboseThis);
    end
end
if R>1
    if numel(fields(fActRun(r,1)))
        fActCat.xMat      = plotDsgnMat(fActCat     ,forceThis,verboseThis);
    end
end



%% Refactor
verboseThis = verbose;
forceThis   = force;
fRespRun = unpackAfni(fRespRun,[],forceThis,verboseThis);
if numel(fields(fActRun))
    fActRun  = unpackAfni(fActRun, [],forceThis,verboseThis);
end
if R>1
    fRespCat = unpackAfni(fRespCat,[],forceThis,verboseThis);
    if numel(fields(fActCat))
        fActCat  = unpackAfni(fActCat, [],forceThis,verboseThis);
    end
else
    fRespCat = fRespRun; fRespCat.r = 0;
    if numel(fields(fActCat))
        fActCat  = fActRun;  fActCat.r = 0;
    end
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



        if ~param.skipMov && fVolTs(1).depth==1
            forceThis = force;
            disp('Making movies')
            nLoop = 4;
            if ~isempty(fMask)
                mask = MRIread(fMask); mask = logical(mask.vol);
            else
                mask = true(fVolTs.height,fVolTs.width,fVolTs.depth);
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





function fRes = runAfni(fList,rR,param,fMask,force,verbose,passDown)
    global src
    if ~exist('rR','var');                     rR = []; end
    if ~exist('fMask','var');               fMask = []; end
    if ~exist('force','var');               force = []; end
    if ~exist('verbose','var');           verbose = []; end
    if ~exist('passDown','var');         passDown = []; end
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
        param.nFrameOrig    = param.nFrameOrig(r,:);
        % param.nDummyRemoved = param.nDummyRemoved(r,:);
    end
    


    %%%%%%%%%%%%%%%%%%%%
    %% Functional design
    %%%%%%%%%%%%%%%%%%%%
    if isfield(param,'dsgn') && isa(param.dsgn,'runDsgn')

        if isempty(param.dsgn.condK); param.dsgn.condK = length(param.dsgn.condLabel); end
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


    sz = size(fList,[1 2 3]);
    multiRun  = sz(1)>1;
    multiEcho = sz(2)>1;
    if diff([sz(3)>1 param.PCflag]); dbstack; error('expecting PCflag or 2 volTs in the 3rd dimension'); end
    

    %%%%%%%
    %% Mask
    %%%%%%%
    if iscell(fMask) && size(fMask,1)>1
        if length(unique(fMask))>1; dbstack; error('multiple masks not supported'); end
        fMask = fMask{1};
    end
    mriMask = MRIread(char(fMask));
    mriMask.vol([1:5 end-4:end],:              ) = 0;
    mriMask.vol(:              ,[1:5 end-4:end]) = 0;




    cmd = {src.afni};
    for E = 1:size(fList,2)

        %% Define files
        fIn = fList(:,E,:);
        if multiRun
            fOut = fIn(1,:,:);
            for iii = 1:size(fOut,3)
                fOut{iii} = strsplit(fOut{iii},'_');
                fOut{iii}{contains(fOut{iii},'run-')} = 'run-cat';
                fOut{iii} = strjoin(fOut{iii},'_'); if ~exist(fileparts(fOut{iii}),'dir'); mkdir(fileparts(fOut{iii})); end
            end
        else
            fOut = fIn;
        end
        

        if param.PCflag
            fStat  = fullfile(fileparts(replace(replace(fOut(:,:,1),'part-real','part-realImag'),'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats']); % let's let afni use its native file format, it is sometimes glitchy otherwise
            fFit   = fullfile(fileparts(replace(replace(fOut(:,:,1),'part-real','part-realImag'),'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_fit.nii.gz']);
        else
            fStat  = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats']); % let's let afni use its native file format, it is sometimes glitchy otherwise
            fFit   = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_fit.nii.gz']);
        end
        if param.getResid
            dbstack; error('double-check that')
            fResid = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_resid.nii.gz']);
        else
            fResid = '';
        end
        fMask  = fullfile(fileparts(replace(replace(fOut(:,:,1),'part-real','part-realImag'),'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_mask.nii.gz' ]);
        if ~exist(fileparts(char(fMask)),'dir'); mkdir(fileparts(char(fMask))); end
        mriMask.fspec = fMask; MRIwrite(mriMask,char(fMask));
        switch HRmodel
            case {'TENT' 'TENTzero'}
                if param.PCflag
                    fResp    = repmat({''},[size(param.dsgn.condLabel,[1 2]) 2]);
                    fRespStd = repmat({''},[size(param.dsgn.condLabel,[1 2]) 2]);
                    for k = 1:param.dsgn.condK
                        fResp(1,k,:)    = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_respAv.nii.gz']);
                        fRespStd(1,k,:) = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_respSd.nii.gz']);
                    end
                else
                    fResp    = repmat({''},size(param.dsgn.condLabel));
                    fRespStd = repmat({''},size(param.dsgn.condLabel));
                    for k = 1:param.dsgn.condK
                        fResp(k)    = cellstr(fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_respAv.nii.gz']));
                        fRespStd(k) = cellstr(fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_respSd.nii.gz']));
                    end
                end
            case {'SPMG2' 'SPMG3'}
                fResp    = [];
                fRespStd = [];
                if param.PCflag; dbstack; error('should not attempt to run SPMG models on PC data'); end
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
            if param.PCflag
                fStim = repmat({''},[size(param.dsgn.condLabel,[1 2]) 2]);
                for k = 1:param.dsgn.condK
                    fStim(:,k,:) = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_startTime.1D']);
                end
                fMat    = fullfile(fileparts(replace(replace(fOut(:,:,1),'part-real','part-realImag'),'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats.xmat.1D']);
                fMatFig = fullfile(fileparts(replace(replace(fOut(:,:,1),'part-real','part-realImag'),'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats.xmat.fig']);
            else
                fStim = cell(size(param.dsgn.condLabel));
                for k = 1:param.dsgn.condK
                    fStim(k) = cellstr(fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_startTime.1D']));
                end
                fMat    = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats.xmat.1D']);
                fMatFig = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats.xmat.fig']);
            end
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
        
        %SPECIAL CASE
        if ~isempty(passDown)
            if strcmp(passDown{end}.id,'TENTzeroParam')
                tmp = cmdTmpTmp(contains(cmdTmpTmp,'TENTzero('));
                if numel(tmp)>1; dbstack; error('code that for more than one stimulus type'); end
                tmp = strsplit(char(tmp),' '); tmp{contains(tmp,'TENTzero(')} = passDown{end}.pr; tmp = {strjoin(tmp,' ')};
                cmdTmpTmp(contains(cmdTmpTmp,'TENTzero(')) = tmp;
            else
                dbstack; error('double-check that')
            end
        end


        if param.dryRun
            system(strjoin([{srcAfni} cmdTmpTmp],newline))
        end

        switch param.model
            case {'TENTzero' 'TENT'}
                if force || anyDontExist([fResp(:)' fRespStd(:)' [char(fStat) '+orig.BRIK']])
                    cmdTmp = [cmdTmp cmdTmpTmp];
                    cmdTmp{end+1} = ['echo ''   ''' strjoin(cellstr(fResp)   ,' ')];
                    cmdTmp{end+1} = ['echo ''   ''' strjoin(cellstr(fRespStd),' ')];
                end
            case {'SPMG2' 'SPMG3'}
                if force || ~exist([char(fStat) '+orig.BRIK'],'file')
                    cmdTmp = [cmdTmp cmdTmpTmp];
                end
            otherwise
                dbstack; error('figure that out')
        end
        cmdTmp{end+1} = ['echo ''   '''                [char(fStat) '+orig']         ];
        cmdTmp{end+1} = ['echo ''   '''                 char(fMat)          ];
        if ~force && ~anyDontExist([fResp(:)' fRespStd(:)' [fStat '+orig.BRIK']])
            cmdTmp{end+1} = 'echo ''   ''already done, skipping';
        end

        fRes(1,E).param = param;
        fRes(1,E).cmd =  strjoin(cmdTmp,newline);

        cmd = [cmd cmdTmp];
    end

    %% Run afni command
    [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || isempty(cmdout) || contains(cmdout,'ERROR','IgnoreCase',false); dbstack; error(cmdout); end




function [cmd,nReg] = afniCmd(fIn,fStim,fMask,param,fResp,fRespStd,fFit,fResid,fMat,fStat,verbose,dryRun)
    % function [cmd,nReg] = afniCmd2(fIn,fMask,fStim,nDummyIgnore,tr,startSeq,durSeq,condSeq,HRmodel,label,param,fResp,fFit,fResid,fMat,fStat,verbose,nDummyRemoved,trDecon,dryRun,nFrame)
    % param.nDummyRemoved [int]: number of initial frames that are already removed from the
    % timeseries. The stimulus timeseries must therefore be adjusted
    % accordingly.
    % param.nDummyIgnore [int]: number of initial frames to ignore from the input
    % timeseries. Only frames after these will be feed to 3dDeconvolve using
    % the [nDummyIgnore..$] notation.
    % nDummy [int]: total number of dummy initial frames
    if ~isfield(param,'forceDeconWin'); param.forceDeconWin = []; end
    tr      = mean(param.tr);
    nFrame  = max(param.nFrame);
    trDecon = param.trDecon;
    nDummyRemoved = param.nFrameOrig - param.nFrame;
    if any(diff(nDummyRemoved)>1); dbstack; error('nDummyRemoved should be the same across runs'); end
    nDummy = mode(param.nDummyIgnore + nDummyRemoved);
    fIn = cellstr(fIn);
    cmd = {'3dDeconvolve -overwrite \'};
    if ~dryRun
        cmd{end+1} = ['-input ' sprintf(['%s[' num2str(param.nDummyIgnore) '..$] '],fIn{:}) ' \'];
        if ~isempty(fMask)
            cmd{end+1} = ['-mask ' char(fMask) ' \'];
        end
    else
        dbstack; error('code that')
        if isempty(nFrame)
            nFrame = MRIread(fIn{1},1); nFrame = nFrame.nframes;
        end
        cmd{end+1} = ['-nodata ' num2str(nFrame) ' ' num2str(tr,'%0.16f') ' \'];
    end
    cmd{end+1} = '-polort A \';
    cmd{end+1} = ['-local_times -stim_times_subtract ' num2str(mean(tr.*nDummy),'%f') ' \'];
    
    % Set design
    dsgn = param.dsgn;
    if param.PCflag
        cmd{end+1} = ['-num_stimts ' num2str(dsgn.condK*2) ' \'];
    else
        cmd{end+1} = ['-num_stimts ' num2str(dsgn.condK) ' \'];
    end
    cmd{end+1} = ['-TR_times ' num2str(trDecon,'%f') ' \'];
    kList = sort(unique(dsgn.cond));
    nReg = zeros(1,dsgn.condK);
    for k = 1:dsgn.condK
        if param.PCflag
            [cmdTmp,nReg(k)] = setAfniStimFileAndCmd(k,kList,dsgn,fIn,fStim(:,k,:),param,dryRun,fResp,fRespStd);
        else
            [cmdTmp,nReg(k)] = setAfniStimFileAndCmd(k,kList,dsgn,fIn,fStim,param,dryRun,fResp,fRespStd);
        end
        cmd = [cmd cmdTmp];
    end


    if dryRun
        dbstack; error('code that')
        fMat = replace(fMat,'_stats.xmat.1D','_stats.xmatPerTrial.1D');
        fMat  = [tempname '_stats.xmat.1D'];
    end



    % Set outputs
    if ~dryRun
        cmd{end+1} = ['-fitts ' char(fFit) ' \'];
        if ~isempty(fResid)
            dbstack; error('code that')
            cmd{end+1} = ['-errts ' char(fResid) ' \'];
        end
        cmd{end+1} = '-bout -fout \';
    end
    %%%%%%%%%%%%%%
    %%% GET AROUND LINUX PATH LENGTH SOFT LIMITATION
    fMatTmp = [tempname '.xmat.1D'];
    cmd{end+1} = ['-x1D ' char(fMatTmp) ' \'];
    %%%%%%%%%%%%%%
    if ~dryRun
        if verbose>0
            cmd{end+1} = ['-bucket ' char(fStat)];
        else
            cmd{end+1} = ['-bucket ' char(fStat) ' 2>/dev/null'];
        end
    else
        cmd{end}(end-1:end) = [];
    end
    %%%%%%%%%%%%%%
    %%% GET AROUND LINUX PATH LENGTH SOFT LIMITATION
    cmd{end+1} = ['cp ' char(fMatTmp) ' ' char(fMat)];
    %%%%%%%%%%%%%%



function [cmd,nReg] = setAfniStimFileAndCmd(k,kList,dsgn,fIn,fStim,param,dryRun,fResp,fRespStd)
    cmd = {};
    if param.PCflag
        cmd{end+1} = ['-stim_label ' num2str(k)            ' ' [char(dsgn.task) '_' dsgn.condLabel{k} 'Real'] ' \'];
        cmd{end+1} = ['-stim_label ' num2str(dsgn.condK+k) ' ' [char(dsgn.task) '_' dsgn.condLabel{k} 'Imag'] ' \'];
    else
        cmd{end+1} = ['-stim_label ' num2str(k) ' ' [char(dsgn.task) '_' dsgn.condLabel{k}] ' \'];
    end

    % write design to file
    if dryRun
        fStim = [tempname '_startTime.1D' ];
    end

    if param.PCflag
        %real
        fido = fopen(fStim{:,:,1}, 'w');
        for i = 1:size(fIn,1)
            fprintf(fido,'%.3f ',dsgn.onsetList(kList(k)==dsgn.cond));
            fprintf(fido,'\n');
        end
        for i = 1:size(fIn,1)
            fprintf(fido,'*');
            fprintf(fido,'\n');
        end
        fclose(fido);
        %imag
        fido = fopen(fStim{:,:,2}, 'w');
        for i = 1:size(fIn,1)
            fprintf(fido,'*');
            fprintf(fido,'\n');
        end
        for i = 1:size(fIn,1)
            fprintf(fido,'%.3f ',dsgn.onsetList(kList(k)==dsgn.cond));
            fprintf(fido,'\n');
        end
        fclose(fido);
    else
        fido = fopen(char(fStim), 'w');
        if ~iscell(fIn); dbstack; error('fIn must be type cell'); end
        for i = 1:size(fIn,1)
            fprintf(fido,'%.3f ',dsgn.onsetList(kList(k)==dsgn.cond));
            fprintf(fido,'\n');
        end
        fclose(fido);
    end
    

    
    % Set model
    switch param.model
        case 'SPMG2'
            dur = dsgn.ondurList(kList(k)==dsgn.cond); if ~isempty(dur) && any(diff(dur)); dbstack; error('stim duration cannot be different across trials'); end
            dur = dur(1);
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
                eTimeNext(end+1) = max(param.nFrame) * mean(param.tr); % use max nFrame in case of a run interruption
                % eTimeNext(end+1) = (param.nFrame + mode(param.nFrameOrig - param.nFrame)) * mean(param.tr);
            else
                eTimeNext = dsgn.onsetList(eTimeNext);
            end
            deconWin = min(eTimeNext - eTime);
            % deconWin = deconWin - 3*tr; % ensure at least one acquisition tr (not trDecon) of baseline between each stimulus
            if (deconWin/param.trDecon)/ceil(deconWin/param.trDecon)>0.9
                deconWin = ceil(deconWin/param.trDecon)*param.trDecon;
            else
                deconWin = floor(deconWin/param.trDecon)*param.trDecon;
            end
            b = 0;
            c = round((deconWin-param.trDecon)/param.trDecon)*param.trDecon;
            nReg = round( (c-b)/param.trDecon + 1 );
            % (c-b)/(nReg-1)
            if param.PCflag
                cmd{end+1} = ['-stim_times ' num2str(k)            ' ' fStim{:,:,1} ' ''TENTzero(' num2str(b) ',' num2str(c) ',' num2str(nReg) ')'' \'];
                cmd{end+1} = ['-stim_times ' num2str(dsgn.condK+k) ' ' fStim{:,:,2} ' ''TENTzero(' num2str(b) ',' num2str(c) ',' num2str(nReg) ')'' \'];
            else
                cmd{end+1} = ['-stim_times ' num2str(k) ' ' char(fStim) ' ''TENTzero(' num2str(b) ',' num2str(c) ',' num2str(nReg) ')'' \'];
            end
            nReg = nReg - 2;
            if ~dryRun
                if param.PCflag
                    cmd{end+1} = ['-iresp ' num2str(k)            ' ' fResp{:,:,1}    ' \'];
                    cmd{end+1} = ['-sresp ' num2str(k)            ' ' fRespStd{:,:,1} ' \'];
                    cmd{end+1} = ['-iresp ' num2str(dsgn.condK+k) ' ' fResp{:,:,2}    ' \'];
                    cmd{end+1} = ['-sresp ' num2str(dsgn.condK+k) ' ' fRespStd{:,:,2} ' \'];
                else
                    cmd{end+1} = ['-iresp ' num2str(k) ' ' char(fResp)    ' \'];
                    cmd{end+1} = ['-sresp ' num2str(k) ' ' char(fRespStd) ' \'];
                end
            end
        otherwise
            dbstak; error('X');
    end