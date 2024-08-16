function [out, info] = volPsdFullMt3(do,info,volTs,volAnat,volPsd)
if isempty(do)
    do.loadIt = 0;
    do.doIt   = 1;
    do.saveIt = 0;
end

if ~exist('volAnat','var'); volAnat = []; end
if ~exist('volPsd','var');   volPsd = []; end % for purpose of a second-step analysis using a mask derived from the first-step analysis

if ~isfield(info,'K');                   info.K = []        ; end
if ~isfield(info,'win');               info.win = zeros(0,2); end
if ~isfield(info,'skipSvd');       info.skipSvd = 0         ; end
if ~isfield(info,'dtrndOrder'); info.dtrndOrder = []        ; end
if ~isfield(info,'onsetList');   info.onsetList = []        ; end
if ~isfield(info,'ondurList');   info.ondurList = []        ; end
if ~isfield(info,'perm');             info.perm = []         ; end

if isempty(info.perm); info.perm = 0; end


%% User variables
outVar = 'volPsd';
stepLabel = 'full multitaper processing';




%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if isfield(info,'outDir'); outDir = info.outDir; else, outDir = info.preprocDir; end; if ~exist(outDir,'dir'); mkdir(outDir); end; stepFile = fullfile(outDir,[strjoin({['sub-' info.sub] ['ses-' info.ses] mfilename},'_')]);

if exist('do','var') && ~isempty(do)
    if isfield(do,'loadIt')  && ~isempty(do.loadIt);   loadIt = do.loadIt;  else,  loadIt=0; end
    if isfield(do,'doIt')    && ~isempty(do.doIt);     doIt   = do.doIt;    else,    doIt=0; end
    if isfield(do,'saveIt')  && ~isempty(do.saveIt);   saveIt = do.saveIt;  else,  saveIt=0; end
    if isfield(do,'writeIt') && ~isempty(do.writeIt); writeIt = do.writeIt; else, writeIt=0; end
else
    loadIt = 0; doIt = 1; saveIt = 0; writeIt = 0;
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
    if ischar(volAnat)
        mask = MRIread(volAnat); mask = mask.vol;
    elseif isstruct(volAnat)
        mask = MRIread(volAnat.f); mask = mask.vol;
    else
        dbstack; error('double-check that');
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
            volTs = vec2vol(vol2vec(volTs));
            for I = 1:length(volAnat)
                if isfield(volAnat(I),'fun') && isfield(volAnat(I).fun,'mask')

                    mask = volAnat(I).fun.mask.crop.vol;

                    if isfield(volAnat(I).fun.mask,'head')
                        mask = mask & volAnat(I).fun.mask.head.mri.vol;
                    elseif isfield(volAnat(I).fun.mask,'brain')
                        mask = mask & volAnat(I).fun.mask.brain.mri.vol;
                    end

                    if ~isempty(volPsd)
                        threshPrctl = 99;
                        Mi = 1;
                        [~,Fi] = min(abs(volPsd.svd.f - 1/mean(diff(volPsd.dsgn.onsets))));
                        threshMask = false(size(volPsd.vol2vec));
                        threshMask(volPsd.vol2vec) = abs(volPsd.svd.spSV(:,:,:,:,Fi,:,:,Mi))>prctile(abs(volPsd.svd.spSV(:,:,:,:,Fi,:,:,Mi)),threshPrctl);
                        % figure('WindowStyle','docked');
                        % imagesc(mask>prctile(abs(volPsd.svd.spSV(:,:,:,:,Fi,:,:,Mi)),threshPrctl))
                        % histogram(abs(volPsd.svd.spSV(:,:,:,:,Fi,:,:,Mi))); xline(prctile(abs(volPsd.svd.spSV(:,:,:,:,Fi,:,:,Mi)),threshPrctl))
                        mask = mask & threshMask;
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
    volTs = applyMask(volTs,mask);
else

end

% %% Add time
% volTs = addTime(volTs);

%% Normalize to thermal noise
% thermalNoiseRange = [0.5 inf];
% modeToRemove = 1:5;
% funTs = normPSD3(funTs,thermalNoiseRange,modeToRemove);

%% Detrend run-by-run
[volTs,~,info.dtrndOrder] = dtrnd2(volTs,[],[],info.dtrndOrder);
% volTsTmp = vec2vol(volTs);
% volTsTmp.vol = volTsTmp.vol - volTsTmp.imMean;


% volTs = vol2vec(volTs);

%% Compute full multitaper spectral analysis
if isempty(info.K);         info.K = 1      ; end
if isempty(info.win);     info.win = [inf 0]; end

K = info.K;
W = [];
winSz   = info.win(1);
if length(info.win)==1 || info.win(1)==inf; info.win(2) = 0; end
winStep = info.win(2);
if winSz~=inf && winStep==0; winStep = 1; end

phaseRand = 0;
taperPerm = info.perm;
% if isempty(info.onsetList)
%     info.onsetList = volTs.dsgn.onsetList';
%     info.ondurList = volTs.dsgn.ondurList';
% end
% volPsd = runFullMT2(volTs,W,K,[winSz winStep],info.onsetList,info.ondurList,[],1,info.skipSvd,[],[],taperPerm,phaseRand);
volPsd = runFullMT3(volTs,W,K,[winSz winStep],[],[],[],1,info.skipSvd,[],[],taperPerm,phaseRand);
volPsd.volAnat = volAnat;

% phaseRand = 2^2;
% taperPerm = 2^2;
% volPsd = runFullMT(volTs,W,K,[winSz winSz/winOver],[],1,[],[],[],taperPerm,phaseRand);


end

if writeIt
%% %%%%%%%%%%%%%
% Write to nii %
%%%%%%%%%%%%% %%
disp('writing to nii')

%%% PSD maps
volPsd.vec = permute(log(volPsd.psd.PSD),[5 6 8 1 2 3 4 7]);
tmpMean = mean(volPsd.vec,2);
volPsd.nframes = size(volPsd.vec,1);
volPsd.tr = mean(diff(volPsd.psd.f))*1000;
tmp = vec2vol(volPsd);
tmp.vol(1:10,1:10,1,:) = repmat(permute(tmpMean,[2 3 4 1]),[10 10 1 1]);
volPsd.psd.fspec = replace(volTs.fspec,'_volTs.nii.gz','_volPsd');
[a,b,~] = fileparts(volPsd.psd.fspec); b = replace(b,'preproc',strjoin({['mask-' volPsd.volAnat.label] ['K-' num2str(volPsd.psd.K)]},'_'));
volPsd.psd.fspec = fullfile(a,[b '.nii.gz']);
% volPsd.psd.fspec = [fullfile(info.preprocDir,strjoin({['sub-' info.sub] ['ses-' info.ses] ['logPsd']},'_')) '.nii.gz'];
MRIwrite(tmp,volPsd.psd.fspec);
disp(volPsd.psd.fspec)

%%% Coherence (first singular vector)
if isfield(volPsd.svd,'spSV') && ~isempty(volPsd.svd.spSV)
    volPsd.vec = permute(abs(volPsd.svd.spSV(:,:,:,:,:,:,:,1)),[5 6 8 1 2 3 4 7]);
    volPsd.nframes = size(volPsd.vec,1);
    volPsd.tr = mean(diff(volPsd.svd.f))*1000;
    tmp = vec2vol(volPsd);
    tmp.vol = abs(tmp.vol);
    tmp.vol(1:10,1:10,1,:)           = repmat(permute(volPsd.svd.COH(:,:,:,:,:,:,:,1)         ,[1 2 3 5 4 6 7 8]),[10 10 1 1]);
    if info.perm
        tmp.vol(end-9:end,end-9:end,1,:) = repmat(permute(volPsd.svd.COH_permMean(:,:,:,:,:,:,:,1),[1 2 3 5 4 6 7 8]),[10 10 1 1]);
    end
    volPsd.svd.fspec.spSVmag = replace(volTs.fspec,'_volTs.nii.gz','_volCohMag');
    [a,b,~] = fileparts(volPsd.svd.fspec.spSVmag); b = replace(b,'preproc',strjoin({['mask-' volPsd.volAnat.label] ['K-' num2str(volPsd.psd.K)]},'_'));
    volPsd.svd.fspec.spSVmag = fullfile(a,[b '.nii.gz']);
    MRIwrite(tmp,volPsd.svd.fspec.spSVmag);
    disp(volPsd.svd.fspec.spSVmag)
    
    tmp = vec2vol(volPsd);
    tmp.vol = angle(tmp.vol);
    tmp.vol(1:10,1:10,1,:)           = repmat(permute(volPsd.svd.COH(:,:,:,:,:,:,:,1)         ,[1 2 3 5 4 6 7 8]),[10 10 1 1]);
    if info.perm
        tmp.vol(end-9:end,end-9:end,1,:) = repmat(permute(volPsd.svd.COH_permMean(:,:,:,:,:,:,:,1),[1 2 3 5 4 6 7 8]),[10 10 1 1]);
    end
    volPsd.svd.fspec.spSVphase = replace(volTs.fspec,'_volTs.nii.gz','_volCohPhase');
    [a,b,~] = fileparts(volPsd.svd.fspec.spSVphase); b = replace(b,'preproc',strjoin({['mask-' volPsd.volAnat.label] ['K-' num2str(volPsd.psd.K)]},'_'));
    volPsd.svd.fspec.spSVphase = fullfile(a,[b '.nii.gz']);
    MRIwrite(tmp,volPsd.svd.fspec.spSVphase);
    disp(volPsd.svd.fspec.spSVphase)
end
if isfield(volPsd.svd,'spSV_pVal')
    dbstack; error('double-check that')
    volPsd.vec = permute(volPsd.svd.spSV_pVal(:,:,:,:,:,:,:,1),[5 6 8 1 2 3 4 7]);
    volPsd.nframes = size(volPsd.vec,1);
    volPsd.tr = mean(diff(volPsd.svd.f))*1000;
    tmp = vec2vol(volPsd);
    tmp.vol(1:10,1:10,1,:)           = 0;
    tmp.vol(end-9:end,end-9:end,1,:) = 0;
    volPsd.svd.fspec.spSVmag_pVal = [fullfile(info.preprocDir,strjoin({['sub-' info.sub] ['ses-' info.ses] ['part-mag'] ['spSVpVal']},'_')) '.nii.gz'];
    MRIwrite(tmp,volPsd.svd.fspec.spSVmag_pVal);
    disp(volPsd.svd.fspec.spSVmag_pVal)
end
if isfield(volPsd.svd,'spSV_fdr')
    dbstack; error('double-check that')
    volPsd.vec = permute(volPsd.svd.spSV_fdr(:,:,:,:,:,:,:,1),[5 6 8 1 2 3 4 7]);
    volPsd.nframes = size(volPsd.vec,1);
    volPsd.tr = mean(diff(volPsd.svd.f))*1000;
    tmp = vec2vol(volPsd);
    tmp.vol(1:10,1:10,1,:)           = 0;
    tmp.vol(end-9:end,end-9:end,1,:) = 0;
    volPsd.svd.fspec.spSVmag_fdr = [fullfile(info.preprocDir,strjoin({['sub-' info.sub] ['ses-' info.ses] ['part-mag'] ['spSVfdr']},'_')) '.nii.gz'];
    MRIwrite(tmp,volPsd.svd.fspec.spSVmag_fdr);
    disp(volPsd.svd.fspec.spSVmag_fdr)
end
if isfield(volPsd.svd,'spSV') && isfield(volPsd.svd,'spSV_fdr')
    dbstack; error('double-check that')
    volPsd.vec = permute(volPsd.svd.spSV(:,:,:,:,:,:,:,1),[5 6 8 1 2 3 4 7]);
    thresh = permute(abs(volPsd.svd.spSV_pVal(:,:,:,:,:,:,:,1)),[5 6 8 1 2 3 4 7]);
    volPsd.vec(thresh>0.05) = 0;
    volPsd.nframes = size(volPsd.vec,1);
    volPsd.tr = mean(diff(volPsd.svd.f))*1000;
    tmp = vec2vol(volPsd);
    tmp.vol = abs(tmp.vol);
    tmp.vol(1:10,1:10,1,:)           = repmat(permute(volPsd.svd.COH(:,:,:,:,:,:,:,1)         ,[1 2 3 5 4 6 7 8]),[10 10 1 1]);
    tmp.vol(end-9:end,end-9:end,1,:) = repmat(permute(volPsd.svd.COH_permMean(:,:,:,:,:,:,:,1),[1 2 3 5 4 6 7 8]),[10 10 1 1]);
    volPsd.svd.fspec.spSVmag_pValThresh = [fullfile(info.preprocDir,strjoin({['sub-' info.sub] ['ses-' info.ses] ['part-mag'] ['spSVpValThresh']},'_')) '.nii.gz'];
    MRIwrite(tmp,volPsd.svd.fspec.spSVmag_pValThresh);
    disp(volPsd.svd.fspec.spSVmag_pValThresh)
    tmp = vec2vol(volPsd);
    tmp.vol = angle(tmp.vol);
    tmp.vol(1:10,1:10,1,:) = repmat(permute(volPsd.svd.COH(:,:,:,:,:,:,:,1),[1 2 3 5 4 6 7 8]),[10 10 1 1]);
    volPsd.svd.fspec.spSVphase_pValThresh = [fullfile(info.preprocDir,strjoin({['sub-' info.sub] ['ses-' info.ses] ['part-phase'] ['spSVpValThresh']},'_')) '.nii.gz'];
    MRIwrite(tmp,volPsd.svd.fspec.spSVphase_pValThresh);
    disp(volPsd.svd.fspec.spSVphase_pValThresh)
end
if isfield(volPsd.svd,'spSV') && isfield(volPsd.svd,'spSV_fdr')
    dbstack; error('double-check that')
    volPsd.vec = permute(volPsd.svd.spSV(:,:,:,:,:,:,:,1),[5 6 8 1 2 3 4 7]);
    thresh = permute(abs(volPsd.svd.spSV_fdr(:,:,:,:,:,:,:,1)),[5 6 8 1 2 3 4 7]);
    volPsd.vec(thresh>0.05) = 0;
    volPsd.nframes = size(volPsd.vec,1);
    volPsd.tr = mean(diff(volPsd.svd.f))*1000;
    tmp = vec2vol(volPsd);
    tmp.vol = abs(tmp.vol);
    tmp.vol(1:10,1:10,1,:)           = repmat(permute(volPsd.svd.COH(:,:,:,:,:,:,:,1)         ,[1 2 3 5 4 6 7 8]),[10 10 1 1]);
    tmp.vol(end-9:end,end-9:end,1,:) = repmat(permute(volPsd.svd.COH_permMean(:,:,:,:,:,:,:,1),[1 2 3 5 4 6 7 8]),[10 10 1 1]);
    volPsd.svd.fspec.spSVmag_fdrThresh = [fullfile(info.preprocDir,strjoin({['sub-' info.sub] ['ses-' info.ses] ['part-mag'] ['spSVfdrThresh']},'_')) '.nii.gz'];
    MRIwrite(tmp,volPsd.svd.fspec.spSVmag_fdrThresh);
    disp(volPsd.svd.fspec.spSVmag_fdrThresh)
    tmp = vec2vol(volPsd);
    tmp.vol = angle(tmp.vol);
    tmp.vol(1:10,1:10,1,:) = repmat(permute(volPsd.svd.COH(:,:,:,:,:,:,:,1),[1 2 3 5 4 6 7 8]),[10 10 1 1]);
    volPsd.svd.fspec.spSVphase_fdrThresh = [fullfile(info.preprocDir,strjoin({['sub-' info.sub] ['ses-' info.ses] ['part-phase'] ['spSVfdrThresh']},'_')) '.nii.gz'];
    MRIwrite(tmp,volPsd.svd.fspec.spSVphase_fdrThresh);
    disp(volPsd.svd.fspec.spSVphase_fdrThresh)
end
volPsd.vec = [];

end




%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

disp(repmat('-',1,length(stepLabel)+6)); disp([upper(stepLabel) ': DONE']); toc; disp(repmat('-',1,length(stepLabel)+6)); disp(' '); disp(' ');
eval(['out = ' outVar '; clear ' outVar]);