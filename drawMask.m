function masks = drawMask(fspec,useSynth,force)
global srcFs srcAfni
if ~exist('force','var');     force = []; end
% if ~exist('verbose','var'); verbose = []; end
if isempty(force);            force = 0; end
% if isempty(verbose);        verbose = 0; end
if ~exist('useSynth','var'); useSynth = []; end
if isempty(useSynth);        useSynth = 0; end




if useSynth
    dbstack; error('double-check that')
    %%% Start with mask from synthstrip
    fspec2 = fspec;
    fSynth = replace(fspec,'_volTs.nii.gz','_volSynthMask.nii.gz');
    [a,b] = fileparts(replace(fspec2,'.nii.gz',''));
    if startsWith(b,'cat_av_') && exist(fullfile(a,['av_' b '.nii.gz']),'file')
        fspec2 = fullfile(a,['av_' b '.nii.gz']);
        fSynth = replace(fspec2,'_volTs.nii.gz','_volSynthMask.nii.gz');
    end
    fMask = replace(fSynth,'_volSynthMask.nii.gz','_volBrainMask.nii.gz');

    if force || ~exist(fSynth,'file') || ~exist(fMask,'file')
        cmd = {srcFs}; cmd{end+1} = srcAfni;
        cmd{end+1} = 'echo mri_synthstrip';
        cmd{end+1} = 'mri_synthstrip \';
        cmd{end+1} = ['-i ' fspec2 ' \'];
        cmd{end+1} = ['-m ' fSynth];
        cmd{end+1} = ['cp ' fSynth ' ' fMask];
        cmd{end+1} = '3dmask_tool -overwrite \';
        cmd{end+1} = ['-input ' fMask ' \'];
        cmd{end+1} = ['-prefix ' fMask ' \'];
        cmd{end+1} = ['-dilate_input 1'];
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end

    %%% Manual drawing
    % fMaskBids = fMask;
    if force || ~exist(fMask,'file')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        disp('Edit brain mask for preprocessing,')
        disp('save it to default name and close window')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        cmd = {srcFs};
        cmd{end+1} = ['fslview -m single ' fspec ' ' fMask];
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
else

    %%% Manual drawing
    fIn = replace(fspec,'_volTs.nii.gz','_volTs-mask.nii.gz');
    fMask = replace(fspec,'_volTs.nii.gz','_volBrainMask.nii.gz');
    
    if force || ~exist(fMask,'file')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        disp('Draw brain mask for preprocessing,')
        disp('save it to default name and close window')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        cmd = {srcFs};
        cmd{end+1} = ['fslview -m single ' fspec ' -b 0,800'];
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
        movefile(fIn,fMask)
    end
end

% %%% Manual drawing
% cmd{end+1} = ['fslview -m single ' fspec ' -b 0,800'];
% 
% fIn = replace(fspec,'_volTs.nii.gz','_volTs-mask.nii.gz');
% fMask = replace(fspec,'_volTs.nii.gz','_volBrainMask.nii.gz');
% fMaskInv = replace(fMask,'_volBrainMask.nii.gz','_volBrainMaskInv.nii.gz');
% 
% if force || ~exist(fMask,'file')
%     disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
%     disp('Draw brain mask for preprocessing,')
%     disp('save it to default name and close window')
%     disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
%     [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
%     movefile(fIn,fMask)
% end

%%% Invert mask
fMaskInv = replace(fMask,'_volBrainMask.nii.gz','_volBrainMaskInv.nii.gz');
if force || ~exist(fMaskInv,'file')
    cmd = {srcAfni};
    cmd{end+1} = '3dcalc -overwrite \';
    cmd{end+1} = ['-prefix ' fMaskInv ' \'];
    cmd{end+1} = ['-a ' fMask ' \'];
    cmd{end+1} = '-expr ''-(a-1)''';
    [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
end

masks.fMask    = fMask;
masks.fMaskInv = fMaskInv;
masks.fUlay    = fspec;

return

%% Image or images to use for creating mask
if isfield(runSet,'finalFiles')
    fIm = fullfile(runSet.bidsDerivDir,'cat_av_preproc_volTs.nii.gz');
    if exist(fIm,'file'); disp('using preprocessed average'); end
    skipHeadFlag = 0;
    % elseif isfield(runSet,'brMocoFiles')
    % elseif isfield(runSet,'wrMocoFiles')
elseif isfield(runSet,'initFiles')
    fIm = fullfile(runSet.initFiles.wd,'cat_av_setPlumb_volTs.nii.gz');
    if exist(fIm,'file'); disp('using unprocessed averages'); end
    skipHeadFlag = 1;
else
    dbstack; error('cannot figure out what to do')
end
if ~exist(fIm,'file'); dbstack; error('image for drawing mask not found'); end



%% Simple crop
volAnat.func.mask.crop.mri = MRIread(fIm,1);
%%% Aim to crop out 'crop'mm in each direction
crop = 5; %mm
trim = ceil(crop./volAnat.func.mask.crop.mri.volres); %voxel
trimActual = trim.*volAnat.func.mask.crop.mri.volres; %mm
volAnat.func.mask.crop.mri.vol = true(volAnat.func.mask.crop.mri.volsize);
volAnat.func.mask.crop.mri.vol([1:trim(1) end-trim(1)+1:end],:                            ,:) = false;
volAnat.func.mask.crop.mri.vol(:                            ,[1:trim(2) end-trim(2)+1:end],:) = false;

trimZ = trim(3);
while volAnat.func.mask.crop.mri.volsize(3) - trimZ*2 < 5
    trimZ = trimZ-1;
end
volAnat.func.mask.crop.mri.vol(:,:,[1:trimZ end-trimZ+1:end]) = false;



%% Brain mask
cmd = {srcFs}; cmd{end+1} = srcAfni;
if info.useSynth
    cmd{end+1} = 'echo mri_synthstrip';
    cmd{end+1} = 'mri_synthstrip \';
    cmd{end+1} = ['-i ' fIm ' \'];
    fOut = replace(fIm,'_volTs.nii.gz','_volSynthMask.nii.gz');
    cmd{end+1} = ['-m ' fOut];

    fIn = fOut;
    fMask = replace(fIn,'_volSynthMask.nii.gz','_volBrainMask.nii.gz');
    fMaskInv = replace(fIn,'_volSynthMask.nii.gz','_volBrainMaskInv.nii.gz');
    cmd{end+1} = ['cp ' fIn ' ' fMask];

    cmd{end+1} = '3dmask_tool -overwrite \';
    cmd{end+1} = ['-input ' fMask ' \'];
    cmd{end+1} = ['-prefix ' fMask ' \'];
    cmd{end+1} = ['-dilate_input 1'];

    cmd{end+1} = 'echo edit brain mask if needed';
    cmd{end+1} = ['fslview ' fIm ' ' fMask];

    if force || ~exist(fMask,'file') || ~exist(fMaskInv,'file')
        % disp('edit brain mask if needed')
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end

        % zero first and last slice
        mask = MRIread(fMask);
        mask.vol(:,:,[1 end]) = 0;
        MRIwrite(mask,fMask);

        % invert mask
        mask.vol = -(mask.vol-1);
        MRIwrite(mask,fMaskInv);
    end

else

    cmd{end+1} = ['fslview -m single ' fIm];

    fIn = replace(fIm,'_volTs.nii.gz','_volTs-mask.nii.gz');
    fMask = replace(fIm,'_volTs.nii.gz','_volBrainMask.nii.gz');
    fMaskInv = replace(fMask,'_volBrainMask.nii.gz','_volBrainMaskInv.nii.gz');

    if force || ~exist(fMask,'file') || ~exist(fMaskInv,'file')
        disp('draw brain mask')
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
        movefile(fIn,fMask)

        cmd = {srcAfni};
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-prefix ' fMaskInv ' \'];
        cmd{end+1} = ['-a ' fMask ' \'];
        cmd{end+1} = '-expr ''-(a-1)''';
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
    end
end


% if ~exist(fOut,'file')
%     if info.useSynth
%         disp('edit brain mask if needed')
%     else
%         disp('draw brain mask')
%     end
%     [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status; dbstack; error(cmdout); error('x'); end
%     if ~info.useSynth
%         movefile(fIn,fOut)
%     end
% end

volAnat.func.mask.brain = volAnat.func.mask.crop;
volAnat.func.mask.brain.mri = MRIread(fMask);
volAnat.func.mask.brainInv = volAnat.func.mask.crop;
volAnat.func.mask.brainInv.mri = MRIread(fMaskInv);

%% Head mask
if ~skipHeadFlag
    cmd = {srcFs};
    cmd{end+1} = ['fslview -m single ' fIm];
    cmd = strjoin(cmd,newline); % disp(cmd)

    fIn = replace(fIm,'_volTs.nii.gz','_volTs-mask.nii.gz');
    fOut = replace(fIm,'_volTs.nii.gz','_volHeadMask.nii.gz');

    if force || ~exist(fOut,'file')
        disp('draw head mask')
        [status,cmdout] = system(cmd,'-echo'); if status; dbstack; error(cmdout); error('x'); end
        movefile(fIn,fOut)
    end

    volAnat.func.mask.head     = volAnat.func.mask.crop;
    volAnat.func.mask.head.mri = MRIread(fOut);
end







%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

disp(repmat('-',1,length(stepLabel)+6)); disp([upper(stepLabel) ': DONE']); toc; disp(repmat('-',1,length(stepLabel)+6)); disp(' '); disp(' ');
eval(['out = ' outVar '; clear ' outVar]);