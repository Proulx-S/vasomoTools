function [out, info] = volTsGetResp(do,info,volTs,volAnat,force)
if isempty(do)
    do.loadIt = 0;
    do.doIt   = 1;
    do.saveIt = 0;
end

if ~exist('force','var');     force = []; end
if ~exist('volAnat','var'); volAnat = []; end
if isempty(force);            force = 0; end

% if ~isfield(info,'K');                   info.K = []        ; end
% if ~isfield(info,'win');               info.win = zeros(0,2); end
% if ~isfield(info,'skipSvd');       info.skipSvd = 0         ; end
% if ~isfield(info,'dtrndOrder'); info.dtrndOrder = []        ; end
% if ~isfield(info,'onsets');   info.onsets = []        ; end
% if isempty(info.onsets);   info.onsets = []        ; end
% if ~isfield(info,'ondurList');   info.ondurList = []        ; end


if ~isfield(volTs,'dsgn');   volTs.dsgn = []        ; end
if ~isempty(volTs.dsgn)
    if isfield(volTs.dsgn,'onsets')
        onsets = volTs.dsgn.onsets;
    else
        onsets = volTs.dsgn.onsetList;
    end
end
if ~isempty(volTs.dsgn)
    if isfield(volTs.dsgn,'ondurs')
        ondurs = volTs.dsgn.ondurs;
    else
        ondurs = volTs.dsgn.ondurList;
    end
end




%% User variables
outVar = 'volResp';
stepLabel = 'event-related response processing';




%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if isfield(info,'outDir'); outDir = info.outDir; else, outDir = info.preprocDir; end; if ~exist(outDir,'dir'); mkdir(outDir); end; stepFile = fullfile(outDir,[strjoin({['sub-' info.sub] ['ses-' info.ses] mfilename},'_')]);

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
    if length(volAnat)==1 && ...
            isfield(volAnat,'mask') && isfield(volAnat.mask,'crop') && isfield(volAnat.mask.crop,'mri') && isfield(volAnat.mask.crop.mri,'vol') && ~isempty(volAnat.mask.crop.mri.vol)
        %%%crop
        mask = volAnat.mask.crop.mri.vol;
        %%%head
        mask = mask & any(volAnat.mask.head.mri.vol,4);
        % %%%brain
        % mask = mask & any(volAnat.mask.brain.mri.vol,4);
        %%%apply
        volTs = applyMask(volTs,mask);
    else

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

if length(volTs)==1 && isfield(volTs,'nDummy') && ~isempty(volTs.nDummy)
    param.nDummy = volTs.nDummy;
else
    param.nDummy = info.dummy;
end
param.nDummyRemoved = param.nDummy;
if volTs.tr/1000 ~= volTs.dsgn.dt
    volTs.tr/1000 
    if isfield(volTs.dsgn,'onsets')
        diff(volTs.dsgn.onsets)
        volTs.dsgn.onsets / (volTs.tr/1000)
        volTs.dsgn.onsets / 4

        warning(strjoin({''...
        ['volume TR   =  ' sprintf('%7.3f ',volTs.tr/1000) 'sec']...
        ['stim dt     =  ' sprintf('%7.3f ',volTs.dsgn.dt) 'sec']...
        ['stim onsets = [' sprintf('%7.3f ',volTs.dsgn.onsets) ']sec']...
        ['            = [' sprintf('%7.3f ',(volTs.dsgn.onsets / (volTs.tr/1000))) ']vol']...
        'Defaulting to deconvolution TR = 1sec'},newline))
        param.trDecon = 1;
    else
        diff(volTs.dsgn.onsetList)
        volTs.dsgn.onsetList / (volTs.tr/1000)
        volTs.dsgn.onsetList / volTs.dsgn.dt
        % volTs.dsgn.onsetList / volTs.dsgn.dt - round(volTs.dsgn.onsetList / volTs.dsgn.dt)

        tError = (volTs.dsgn.onsetList(end) / (volTs.tr/1000) - volTs.dsgn.onsetList(end) / (volTs.dsgn.dt)) * 1000;

        warning(strjoin({''...
        ['volume TR   =  ' sprintf('%7.3f ',volTs.tr) 'ms']...
        ['stim dt     =  ' sprintf('%7.3f ',volTs.dsgn.dt*1000) 'ms']...
        ['stim onsets = [' sprintf('%7.0f ',volTs.dsgn.onsetList*1000) ']ms']...
        ['            = [' sprintf('%7.3f ',(volTs.dsgn.onsetList / (volTs.tr/1000))) ']vol']...
        ''...
        [num2str(tError) 'ms error by the last onset']},newline))

        if abs(tError)<10
            param.trDecon = volTs.dsgn.dt;
        else
            dbstack; error('timing problem here')
        end
        

        
    end
    
else
    param.trDecon = volTs.tr/1000; %volTs.tr/1000/4; %1
end
param.verbose = 1;
forceThis = force;
switch info.dataSetLabel
    case 'vsmDriven'
        [files,fRun,fSes,fSes_echoCat,param] = getResp(volTs,[],param,forceThis);
    otherwise
        dbstack; error('double-check that');
        [files,fRun,fSes,fSes_echoCat,param] = getResp(volTs,volAnat,param,forceThis);
end
volResp.ts = MRIread(files.resp.f{1});
volResp.tsOnBase = MRIread(files.respOnBase.f{1});
volResp.base = MRIread(files.base.f{1});
volResp.F  = MRIread(files.respF.f{1});
volResp.Fq = MRIread(files.respF_fdr.f{1});
volResp.vid = files.respOnBaseMovieHighBit.f{1};


end




%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

disp(repmat('-',1,length(stepLabel)+6)); disp([upper(stepLabel) ': DONE']); toc; disp(repmat('-',1,length(stepLabel)+6)); disp(' '); disp(' ');
eval(['out = ' outVar '; clear ' outVar]);