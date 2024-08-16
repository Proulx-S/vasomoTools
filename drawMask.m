function out = drawMask(ulay,useSynth)
global srcFs
% if ~exist('force','var');     force = []; end
% if ~exist('verbose','var'); verbose = []; end
% if isempty(force);            force = 0; end
% if isempty(verbose);        verbose = 0; end
if ~exist('useSynth','var'); useSynth = []; end
if isempty(useSynth);        useSynth = 0; end


ulay

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






end






%%%%%%%%%%%%%%%%
%% House keeping
%%%%%%%%%%%%%%%%
if saveIt; disp(strjoin({[upper(stepLabel) ': saving to '] [stepFile '.mat']},newline)); tmp = whos(outVar); if tmp.bytes/1e9<2; save(stepFile,outVar); else, save(stepFile,outVar,'-v7.3'); end; disp([upper(stepLabel) ': saved']); end

disp(repmat('-',1,length(stepLabel)+6)); disp([upper(stepLabel) ': DONE']); toc; disp(repmat('-',1,length(stepLabel)+6)); disp(' '); disp(' ');
eval(['out = ' outVar '; clear ' outVar]);