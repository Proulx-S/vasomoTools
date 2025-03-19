function [out1, out2, out3, out4, info] = volTsGetResp5(do,info,volTs,dsgn,volAnat,force,verbose)
% global srcAfni
% Assuming same dsgn and volAnat for every run
if isempty(do)
    do.loadIt = 0;
    do.doIt   = 1;
    do.saveIt = 0;
end

readFlag = 0;
if ~isfield(info,'doCat'); info.doCat = []; end
if ~isfield(info,'doRun'); info.doRun = []; end
if ~isfield(info,'doMov'); info.doMov = []; end
if isempty(info.doCat);    info.doCat = 1; end
if isempty(info.doRun);    info.doRun = 0; end
if isempty(info.doMov);    info.doMov = 1; end

if ~isMRI(volTs) && isfield(volTs,'mri')
    if isfield(volTs,'dryRun')
        dryRun = [volTs.dryRun];
    else
        dryRun = false(size(volTs));
    end
    volTs = [volTs.mri]';
    [volTs.dryRun] = deal(dryRun);
end

if ~exist('dsgn','var');       dsgn = []; end
if ~exist('volAnat','var'); volAnat = []; end
if ~exist('verbose','var'); verbose = []; end
if ~exist('force','var');     force = []; end
if isempty(force);     force = 0; end
if isempty(verbose); verbose = 0; end

if isempty(dsgn)
    if isfield(volTs,'dsgn')
        if numel(volTs)>1 && ~isequal(volTs.dsgn); dbstack; error('different dsgn for different volTs runs. Code that'); end
        dsgn = volTs(1).dsgn;
    else
        dbstack; error('must provide either ''volTs.dsgn'' or ''dsgn''');
    end
end
if isempty(volAnat)
    if isfield(volTs,'volAnat')
        if ~isequaln(volTs.dsgn); dbstack; error('different volAnat for different volTs runs. Code that'); end
        volAnat = volTs(1).volAnat;
    end
end



% if ~isfield(info,'K');                   info.K = []        ; end
% if ~isfield(info,'win');               info.win = zeros(0,2); end
% if ~isfield(info,'skipSvd');       info.skipSvd = 0         ; end
% if ~isfield(info,'dtrndOrder'); info.dtrndOrder = []        ; end
% if ~isfield(info,'onsets');   info.onsets = []        ; end
% if isempty(info.onsets);   info.onsets = []        ; end
% if ~isfield(info,'ondurList');   info.ondurList = []        ; end


% if ~isfield(volTs,'dsgn');   volTs.dsgn = []        ; end
% if ~isempty(volTs.dsgn)
%     if isfield(volTs.dsgn,'onsets')
%         onsets = volTs.dsgn.onsetList;
%     else
%         onsets = volTs.dsgn.onsetList;
%     end
% end
% if ~isempty(volTs.dsgn)
%     if isfield(volTs.dsgn,'ondurs')
%         ondurs = volTs.dsgn.ondurs;
%     else
%         ondurs = volTs.dsgn.ondurList;
%     end
% end




%% User variables
outVar1 = 'volRespCat';
outVar2 = 'volActCat';
outVar3 = 'volResp';
outVar4 = 'volAct';
% if info.doCat
%     outVar2 = 'volRespSes';
% end
stepLabel = 'event-related response processing';



%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if exist('do','var') && ~isempty(do)
    if isfield(do,'loadIt') && ~isempty(do.loadIt); loadIt = do.loadIt; else, loadIt=0; end
    if isfield(do,'doIt')   && ~isempty(do.doIt);   doIt   = do.doIt;   else, doIt=0;   end
    if isfield(do,'saveIt') && ~isempty(do.saveIt); saveIt = do.saveIt; else, saveIt=0; end
else
    loadIt = 0; doIt = 1; saveIt = 0;
end

if saveIt
    if isfield(info,'outDir'); outDir = info.outDir; else, outDir = info.preprocDir; end; if ~exist(outDir,'dir'); mkdir(outDir); end
    if isfield(info,'ses')
        stepFile = fullfile(outDir,[strjoin({['sub-' info.sub] ['ses-' info.ses] mfilename},'_')]);
    else
        stepFile = fullfile(outDir,[strjoin({['sub-' info.sub] mfilename},'_')]);
    end
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

% % %% Add time
% % volTs = addTime(volTs);
% 
% %% Normalize to thermal noise
% % thermalNoiseRange = [0.5 inf];
% % modeToRemove = 1:5;
% % funTs = normPSD3(funTs,thermalNoiseRange,modeToRemove);
% 
% %% Detrend run-by-run
% [volTs,~,info.dtrndOrder] = dtrnd2(volTs,[],[],info.dtrndOrder);
% % volTsTmp = vec2vol(volTs);
% % volTsTmp.vol = volTsTmp.vol - volTsTmp.imMean;
% 
% 
% % volTs = vol2vec(volTs);

if any(diff([volTs.nFrame])) || any(diff([volTs.nFrameOrig]))
    error('not all runs have the same number of frames')
end
param.nDummyRemoved = volTs(1).nFrameOrig - volTs(1).nFrame;


% if length(volTs)==1 && isfield(volTs,'nDummy') && ~isempty(volTs.nDummy)
%     param.nDummy = volTs.nDummy;
% else
%     param.nDummy = info.dummy;
% end
% param.nDummyRemoved = param.nDummy;

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



switch info.dataSetLabel
    case {'vsmDriven' 'vsmRing' 'zhangxuanPC'}
        if isfield(dsgn,'nullTrial') && ~isempty(dsgn.nullTrial)
            dsgn.condList = dsgn.nullTrial + 1;
            dsgn.condLabelList = {'stim' 'catch'}';
        elseif isempty(dsgn.cond)
            dsgn.cond = ones(size(dsgn.onsetList));
            dsgn.condLabel = {'stim'}';
        end
        % % nullTrialControl = [floor(find(dsgn.nullTrial)/2) find(dsgn.nullTrial) + ceil((length(dsgn.nullTrial) - find(dsgn.nullTrial))/2)];
        % % dsgn.condList(nullTrialControl) = 3;
        % % dsgn.condLabelList = {'stim' 'catch' 'catchCtrl'}';
        %
        % % roiInd = [volAnat.roi{:}]; roiInd = {roiInd.label}';
        % % roiInd = ismember(roiInd,'vesselCalcarine');
        % % [files,fSes,~,~,param_getResp2] = getResp3(volTs,dsgn,volAnat.roi{roiInd}.f,param,forceThis,verboseThis);
        % [files,fSes,~,~,param_getResp2] = getResp3(volTs,dsgn,mask,param,forceThis,verboseThis);
        % % [files,fSes,~,~,param_getResp2] = getResp3(volTs,dsgn,volAnat.mask.head.f,param,forceThis,verboseThis);
        %
        % % % % % % % % % [files,fRun,fSes,fSes_echoCat,param] =
        % % % % % % % % % getRespTmp(volTs,dsgn,volAnat.mask.head.f,param,forceThis,verboseThis);
        % % % % % % % % % % getRespTmp.m implements multiple FIR (more than one response
        % % % % % % % % % time course, one per stimulus condition, e.g. here one for
        % % % % % % % % % regular trials and one for the catch trial). The difficulty here
        % % % % % % % % % is that it is hard to get a F test for individual stimulus
        % % % % % % % % % condition--by default there is only the omnibus F test. The
        % % % % % % % % % solution might be to use glt but I'm not well verse on that.
        % % % % % % % %
        % % % % % % % % % %% Response from each runs and concatenated runs
        % % % % % % % % % for rr = 1:length(volTs)
        % % % % % % % % %     [files,fRun,fSes,fSes_echoCat,param] = getResp2(volTs(rr),dsgn,volAnat.mask.head.f,param,forceThis);
        % % % % % % % % % end
        % % % % % % % % % %% Response from concatenated runs
        % % % % % % % % % [files,fRun,fSes,fSes_echoCat,param] = getResp2(volTs,dsgn,volAnat.mask.head.f,param,forceThis,verboseThis);




        %% Response timecourse
        param.skipCat = ~info.doCat;
        param.skipRun = ~info.doRun;
        param.skipMov = ~info.doMov;
        param.model = 'TENTzero';
        forceThis = force;
        verboseThis = verbose;
        if isfield(info,'dryRun')
            param.dryRun = info.dryRun;
        end
        param.nFrame = [volTs.nFrame]';
        [volResp,volRespCat,~,param_getResp] = getAct4(volTs,dsgn,fMask,param,forceThis,verboseThis);
        volRespCat.afni.param = param_getResp;
        for R = 1:size(volResp,1)
            volResp(R,1).afni.param = param_getResp;
        end

        %%% Extract trial error
        % This redistribute the timeseries squared error to each regressors
        % on a trial-by-trial basis.
        % !!!Note that I'm not sure how this would behave if regressors
        % !!!overlap--more than one regressor fitting a given time point in
        % !!!the timeseries. For that one would probably need to take into
        % !!!account the estimated coefficients for an appropriately weighted
        % !!!redistribution of error--the higher regressor should eat up more
        % !!!of the error.
        forceThis   = force;
        verboseThis = verbose;
        disp('Computing response error')
        if ~param.skipCat && ~param.skipRun
            fCatTrialEr = replace(volRespCat.fs.fRespTs,'_resp.nii.gz','_respTrialSSE.nii.gz');
            fCatTrialN  = replace(volRespCat.fs.fRespTs,'_resp.nii.gz','_respTrialN.mat');
            erCat = [];
            nCat  = [];
        end
        if ~param.skipRun
            mustDoIt = false;
            fTrialEr = cell(size(volResp));
            fTrialN  = cell(size(volResp));
            for R = 1:size(volResp,1)
                nReg = volResp(R).afni.param.perTrialXmat(R).nReg;
                fTrialN{R} = replace(volResp(R,1).fs.fRespTs,'_resp.nii.gz','_respTrialN.mat');
                for e = 1:length(nReg)
                    fTrialEr{R}{1,e} = replace(volResp(R,1).fs.fRespTs,'_resp.nii.gz',['_respTrial' num2str(e,'%03i') 'Er.nii.gz']);
                    if ~exist(fTrialEr{R}{1,e},'file')
                        mustDoIt = true;
                    end
                end
                if ~exist(fTrialN{R},'file')
                    mustDoIt = true;
                end
            end
            if forceThis || mustDoIt || ~exist(fCatTrialEr,'file') || ~exist(fCatTrialN,'file')
                for R = 1:size(volResp,1)
                    disp([' ' num2str(R) '/' num2str(size(volResp,1))])
                    %%% get per trial design matrix
                    xMat = volResp(R).afni.param.perTrialXmat(R).mat;
                    nReg = volResp(R).afni.param.perTrialXmat(R).nReg;
                    xMat = xMat(:,end-sum(nReg)+1:end);
                    xMat2 = cell(size(nReg));
                    for e = 1:length(nReg)
                        xMat2{e} = xMat(:,1:nReg(e));
                        xMat(:,1:nReg(e)) = [];
                    end
                    xMat = xMat2; clear xMat2

                    %%% get residual
                    mriResid = MRIload3(volResp(R).afni.fResid,volResp(R).fs.fMask,[],0);
                    resid = mriResid.vec;

                    %%% get residual per trial
                    for e = 1:length(nReg)
                        if e==1
                            n  = zeros([size(xMat{e},2) 1             length(nReg)]);
                            er = zeros([size(xMat{e},2) size(resid,2) length(nReg)]);
                        end
                        n (:,1,e) = sum(xMat{e},1)';
                        er(:,:,e) = xMat{e}'*resid;
                    end
                    er = cat(1,zeros([1 size(er,[2 3])]),er,zeros([1 size(er,[2 3])]));
                    n  = cat(1,zeros([1,size(n ,[2 3])]),n ,zeros([1 size(n ,[2 3])]));

                    for e = 1:length(nReg)
                        mriResid.vec = [];
                        mriEr = mriResid;
                        mriEr.vec = er(:,:,e);
                        mriEr.nframes = size(er(:,:,e),1);
                        MRIwrite(vec2vol(mriEr),fTrialEr{R}{e});
                    end
                    save(fTrialN{R},'n');

                    % if ~param.skipCat
                    %     if R==1
                    %         erCat = er;
                    %         nCat  = n  ;
                    %     else
                    %         erCat = erCat + er;
                    %         nCat   = nCat + n  ;
                    %     end
                    % end
                    disp('done')
                end
            else
                disp('already done, skipping')
            end
            for R = 1:size(volResp,1)
                volResp(R,1).fs.fRespEr = fTrialEr{R}';
                volResp(R,1).fs.fRespN  = fTrialN{R};
            end
        end

        if ~param.skipCat && ~param.skipRun
            % if forceThis || ~exist(fCatTrialEr,'file') || ~exist(fCatTrialN,'file')
            %     mriResid.vec = [];
            %     mriEr = mriResid;
            %     mriEr.vec = erCat;
            %     mriEr.nframes = size(erCat,1);
            %     MRIwrite(vec2vol(mriEr),fCatTrialEr);
            %     n = nCat;
            %     save(fCatTrialN,'n');
            % end
            volRespCat.fs.fRespEr = [fTrialEr{:}]';
            volRespCat.fs.fRespN  = fTrialN ;
        end



        %% SPM double gamma fit (gamma variate + d/dt derivative)
        param.skipCat = ~info.doCat;
        param.skipRun = ~info.doRun;
        param.skipMov = ~info.doMov;
        param.model = 'SPMG2';
        forceThis   = force;
        verboseThis = verbose;
        if isfield(info,'dryRun')
            param.dryRun = info.dryRun;
        end
        [volAct,volActCat,~,param_getAct] = getAct2(volTs,dsgn,mask.fspec,param,forceThis,verboseThis);
        volActCat.afni.param = param_getAct;
        for R = 1:size(volResp,1)
            volAct(R,1).afni.param = param_getAct;
        end

        % [filesAct,fRunAct,fSesAct,~,param_getAct] = getAct(volTs,dsgn,volAnat.mask.head.f,param,forceThis,verboseThis);
    otherwise
        dbstack; error('double-check that');
        forceThis = force;
        verboseThis = verbose;
        [files,fRun,fSes,fSes_echoCat,param] = getResp(volTs,volAnat,param,forceThis);
end


% volResp.resp = fSesResp;
% volResp.act  = fSesAct;
% % % 
% % % 
% % % %% Summarize output
% % % for c = 1:size(files.resp.f,2)
% % %     volResp.ts(1,c) = MRIread(files.resp.f{1,c},~readFlag);
% % %     volResp.tsOnBase(1,c) = MRIread(files.respOnBase.f{1,c},~readFlag);
% % %     volResp.vid{1,c} = files.respOnBaseMovieHighBit.f{1,c};
% % % end
% % % volResp.base = MRIread(char(files.base.f));
% % % volResp.F  = MRIread(char(files.respF.f));
% % % volResp.Fq = MRIread(char(files.respF_fdr.f));
% % % 
% % % %design matrix for reconstruction
% % % cmd = {srcAfni};
% % % cmd{end+1} = '3dinfo \';
% % % cmd{end+1} = ['-label ' char(files.stat.f) ' \'];
% % % [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
% % % coefLabel = strsplit(cmdout,newline); coefLabel = strsplit(coefLabel{1},'|');
% % % coefLabel = coefLabel(endsWith(coefLabel,'_Coef'));
% % % coefLabel2 = {};
% % % coefLabel2{1,end+1} = coefLabel(contains(coefLabel,'Pol#'));
% % % coefLabelOrig2 = {};
% % % coefLabelOrig2{1,end+1} = coefLabel(contains(coefLabel,'Pol#'));
% % % coefLabel = coefLabel(~contains(coefLabel,'Pol#'));
% % % coefLabelOrig = coefLabel;
% % % for i = 1:length(coefLabel)
% % %     coefLabel{i} = strsplit(coefLabel{i},'#');
% % %     coefLabel{i} = coefLabel{i}{1};
% % % end
% % % while ~isempty(coefLabel)
% % %     if strcmp(coefLabel2{1,end}{1},coefLabel{1})
% % %         coefLabel2{1,end}(end+1) = coefLabel(1); coefLabel(1) = [];
% % %         coefLabelOrig2{1,end}(end+1) = coefLabelOrig(1); coefLabelOrig(1) = [];
% % %     else
% % %         coefLabel2{1,end+1}(1) = coefLabel(1); coefLabel(1) = [];
% % %         coefLabelOrig2{1,end+1}(1) = coefLabelOrig(1); coefLabelOrig(1) = [];
% % %     end
% % % end
% % % 
% % % dsgn.dsgnMat = param_getResp2.funDsgn.mat;
% % % dsgn.dsgnMatLabel = cat(1, cat(2,coefLabel2{:}), cat(2,coefLabelOrig2{:}) );
% % % volResp.dsgn = dsgn;
% % % 
% % % %coef for reconstruction
% % % in = char(files.stat.f);
% % % out = replace(in,'_stats.nii.gz','_coefs.nii.gz');
% % % cmd = {srcAfni};
% % % cmd{end+1} = '3dcalc -overwrite \';
% % % cmd{end+1} = ['-prefix ' out ' \'];
% % % cmd{end+1} = ['-a ' in '[' strjoin(dsgn.dsgnMatLabel(end,:),',') '] \'];
% % % cmd{end+1} = '-expr a ';
% % % [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
% % % 
% % % volResp.coef = MRIread(out,1);
% % % 
% % % 
% % % 
% % % % for rr = 1:size(files.resp.f,1)
% % % %     volRespRun(rr,1).ts = MRIread(files.resp.f{rr,1},~readFlag);
% % % %     volRespRun(rr,1).tsOnBase = MRIread(files.respOnBase.f{rr,1},~readFlag);
% % % %     volRespRun(rr,1).base = MRIread(files.base.f{rr,1});
% % % %     volRespRun(rr,1).F  = MRIread(files.respF.f{rr,1});
% % % %     volRespRun(rr,1).Fq = MRIread(files.respF_fdr.f{rr,1});
% % % %     volRespRun(rr,1).vid = files.respOnBaseMovieHighBit.f{rr,1};
% % % % end
% % % 
% % % 
% % % % for rr = 1:size(filesAct.coef.f,1)
% % % %     volRespRun(rr,1).(param_getAct.model).coef    = MRIread(filesAct.coef.f{rr,1},~readFlag);
% % % %     volRespRun(rr,1).(param_getAct.model).coefPol = MRIread(filesAct.coefPol.f{rr,1},~readFlag);
% % % %     volRespRun(rr,1).(param_getAct.model).F       = MRIread(filesAct.F.f{rr,1});
% % % %     volRespRun(rr,1).(param_getAct.model).Fq      = MRIread(filesAct.F_fdr.f{rr,1});
% % % % end
% % % 
% % % 
% % % % if ~isempty(fSes)
% % % %     volRespSes.ts = MRIread(char(files.resp.fSes),~readFlag);
% % % %     volRespSes.tsOnBase = MRIread(char(files.respOnBase.fSes),~readFlag);
% % % %     volRespSes.base = MRIread(char(files.base.fSes));
% % % %     volRespSes.baseCat = MRIread(char(files.base.fCat));
% % % %     volRespSes.F  = MRIread(char(files.respF.fSes));
% % % %     volRespSes.Fcat  = MRIread(char(files.respF.fCat));
% % % %     volRespSes.Fq = MRIread(char(files.respF_fdr.fSes));
% % % %     volRespSes.FqCat = MRIread(char(files.respF_fdr.fCat));
% % % %     volRespSes.vid = char(files.respOnBaseMovieHighBit.fSes);
% % % % else
% % % %     volRespSes = [];
% % % % end


end




%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

disp(repmat('-',1,length(stepLabel)+6)); disp([upper(stepLabel) ': DONE']); toc; disp(repmat('-',1,length(stepLabel)+6)); disp(' '); disp(' ');
% eval(['out1 = ' outVar1 '; clear ' outVar1]);
% eval(['out2 = ' outVar2 '; clear ' outVar2]);
if exist(outVar1,'var')
    eval(['out1 = ' outVar1 '; clear ' outVar1]);
else
    eval('out1 = [];');
end
if exist(outVar2,'var')
    eval(['out2 = ' outVar2 '; clear ' outVar2]);
else
    eval('out2 = [];');
end
if exist(outVar3,'var')
    eval(['out3 = ' outVar3 '; clear ' outVar3]);
else
    eval('out3 = [];');
end
if exist(outVar4,'var')
    eval(['out4 = ' outVar4 '; clear ' outVar4]);
else
    eval('out4 = [];');
end
