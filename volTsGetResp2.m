function [out1, out2, info] = volTsGetResp2(do,info,volTs,dsgn,volAnat,force,verbose)
% Assuming same dsgn and volAnat for every run
if isempty(do)
    do.loadIt = 0;
    do.doIt   = 1;
    do.saveIt = 0;
end

readFlag = 0;
if ~isfield(info,'doCat'); info.doCat = []; end
if isempty(info.doCat);    info.doCat = 0; end

if ~isMRI(volTs) && isfield(volTs,'mri')
    volTs = [volTs.mri]';
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
outVar1 = 'volRespRun';
if info.doCat
    outVar2 = 'volRespSes';
end
stepLabel = 'event-related response processing';




%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if isfield(info,'outDir'); outDir = info.outDir; else, outDir = info.preprocDir; end; if ~exist(outDir,'dir'); mkdir(outDir); end
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

%% Mask
if ~isempty(volAnat)
    if length(volAnat)==1
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

if all(diff([volTs.nDummy])==0)
    param.nDummy = volTs.nDummy;
end
if all(diff([volTs.nDummyRemoved])==0)
    param.nDummyRemoved = volTs.nDummyRemoved;
end

% if length(volTs)==1 && isfield(volTs,'nDummy') && ~isempty(volTs.nDummy)
%     param.nDummy = volTs.nDummy;
% else
%     param.nDummy = info.dummy;
% end
% param.nDummyRemoved = param.nDummy;

if ~all(diff([volTs.tr])<0.01); dbstack; error('runs have different tr'); end
tr = volTs(1).tr/1000;


param.trDecon = tr;
if tr ~= dsgn.dt
    % if isfield(dsgn,'onsetList')
        % diff(dsgn.onsetList)
        % dsgn.onsetList / tr

        warning(strjoin({''...
        ['volume TR   =  ' sprintf('%7.6f ',tr) 'sec']...
        ['stim dt     =  ' sprintf('%7.6f ',dsgn.dt) 'sec']...
        ['stim onsets = [' sprintf('%7.3f ',dsgn.onsetList) ']sec']...
        ['            = [' sprintf('%7.3f ',(dsgn.onsetList / tr)) ']vol']...
        ['Defaulting to deconvolution TR = ' num2str(param.trDecon,'%7.6f') 'sec']},newline))
    % else
    %     diff(volTs.dsgn.onsetList)
    %     volTs.dsgn.onsetList / (volTs.tr/1000)
    %     volTs.dsgn.onsetList / volTs.dsgn.dt
    %     % volTs.dsgn.onsetList / volTs.dsgn.dt - round(volTs.dsgn.onsetList / volTs.dsgn.dt)
    % 
    %     tError = (volTs.dsgn.onsetList(end) / (volTs.tr/1000) - volTs.dsgn.onsetList(end) / (volTs.dsgn.dt)) * 1000;
    % 
    %     warning(strjoin({''...
    %     ['volume TR   =  ' sprintf('%7.3f ',volTs.tr) 'ms']...
    %     ['stim dt     =  ' sprintf('%7.3f ',volTs.dsgn.dt*1000) 'ms']...
    %     ['stim onsets = [' sprintf('%7.0f ',volTs.dsgn.onsetList*1000) ']ms']...
    %     ['            = [' sprintf('%7.3f ',(volTs.dsgn.onsetList / (volTs.tr/1000))) ']vol']...
    %     ''...
    %     [num2str(tError) 'ms error by the last onset']},newline))
    % 
    %     if abs(tError)<10
    %         param.trDecon = volTs.dsgn.dt;
    %     else
    %         dbstack; error('timing problem here')
    %     end
    % 

        
    % end
end



forceThis = force;
verboseThis = 1;
switch info.dataSetLabel
    case 'vsmDriven'
        %% Response from each runs
        param.skipCat = ~info.doCat;
        [files,fRun,fSes,~,param_getResp2] = getResp2(volTs,dsgn,volAnat.mask.head.f,param,forceThis,verboseThis);
        
        % % % % % % % % [files,fRun,fSes,fSes_echoCat,param] =
        % % % % % % % % getRespTmp(volTs,dsgn,volAnat.mask.head.f,param,forceThis,verboseThis);
        % % % % % % % % % getRespTmp.m implements multiple FIR (more than one response
        % % % % % % % % time course, one per stimulus condition, e.g. here one for
        % % % % % % % % regular trials and one for the catch trial). The difficulty here
        % % % % % % % % is that it is hard to get a F test for individual stimulus
        % % % % % % % % condition--by default there is only the omnibus F test. The
        % % % % % % % % solution might be to use glt but I'm not well verse on that.
        % % % % % % % 
        % % % % % % % % %% Response from each runs and concatenated runs
        % % % % % % % % for rr = 1:length(volTs)
        % % % % % % % %     [files,fRun,fSes,fSes_echoCat,param] = getResp2(volTs(rr),dsgn,volAnat.mask.head.f,param,forceThis);
        % % % % % % % % end
        % % % % % % % % %% Response from concatenated runs
        % % % % % % % % [files,fRun,fSes,fSes_echoCat,param] = getResp2(volTs,dsgn,volAnat.mask.head.f,param,forceThis,verboseThis);





        %% SPM double gamma fit (gamma variate + d/dt derivative)
        param.model = 'SPMG2';
        [filesAct,fRunAct,fSesAct,~,param_getAct] = getAct(volTs,dsgn,volAnat.mask.head.f,param,forceThis,verboseThis);


    otherwise
        dbstack; error('double-check that');
        [files,fRun,fSes,fSes_echoCat,param] = getResp(volTs,volAnat,param,forceThis);
end


for rr = 1:size(files.resp.f,1)
    volRespRun(rr,1).ts = MRIread(files.resp.f{rr,1},~readFlag);
    volRespRun(rr,1).tsOnBase = MRIread(files.respOnBase.f{rr,1},~readFlag);
    volRespRun(rr,1).base = MRIread(files.base.f{rr,1});
    volRespRun(rr,1).F  = MRIread(files.respF.f{rr,1});
    volRespRun(rr,1).Fq = MRIread(files.respF_fdr.f{rr,1});
    volRespRun(rr,1).vid = files.respOnBaseMovieHighBit.f{rr,1};
end
for rr = 1:size(filesAct.coef.f,1)
    volRespRun(rr,1).(param_getAct.model).coef    = MRIread(filesAct.coef.f{rr,1},~readFlag);
    volRespRun(rr,1).(param_getAct.model).coefPol = MRIread(filesAct.coefPol.f{rr,1},~readFlag);
    volRespRun(rr,1).(param_getAct.model).F       = MRIread(filesAct.F.f{rr,1});
    volRespRun(rr,1).(param_getAct.model).Fq      = MRIread(filesAct.F_fdr.f{rr,1});
end
if ~isempty(fSes)
    volRespSes.ts = MRIread(char(files.resp.fSes),~readFlag);
    volRespSes.tsOnBase = MRIread(char(files.respOnBase.fSes),~readFlag);
    volRespSes.base = MRIread(char(files.base.fSes));
    volRespSes.baseCat = MRIread(char(files.base.fCat));
    volRespSes.F  = MRIread(char(files.respF.fSes));
    volRespSes.Fcat  = MRIread(char(files.respF.fCat));
    volRespSes.Fq = MRIread(char(files.respF_fdr.fSes));
    volRespSes.FqCat = MRIread(char(files.respF_fdr.fCat));
    volRespSes.vid = char(files.respOnBaseMovieHighBit.fSes);
else
    volRespSes = [];
end


end




%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

disp(repmat('-',1,length(stepLabel)+6)); disp([upper(stepLabel) ': DONE']); toc; disp(repmat('-',1,length(stepLabel)+6)); disp(' '); disp(' ');
eval(['out1 = ' outVar1 '; clear ' outVar1]);
if info.doCat
    eval(['out2 = ' outVar2 '; clear ' outVar2]);
else
    eval(['out2 = [];']);
end