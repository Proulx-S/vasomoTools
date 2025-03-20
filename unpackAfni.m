function fRes = unpackAfni(fRes,fMask,force,verbose)
global src
if ~exist('fMask','var');     fMask = []; end
if ~exist('force','var');     force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force);     force = 0; end
if isempty(verbose); verbose = 0; end        
if numel(fRes)>1
    for i = 1:numel(fRes)
        fRes(i) = unpackAfni(fRes(i),fMask,force,verbose);
    end
    return;
end
if isempty(fMask); fMask = fRes.fMask; end

stats.fStat = fRes.fStat;
if isempty(fMask)
    stats.fMask = fRes.fMask;
else
    stats.fMask = fMask;
end

stats.model    = fRes.param.model;
stats.task     = fRes.param.dsgn.task;
stats.condList = fRes.param.dsgn.condLabel;

cmd = {src.afni};

%% Extract baseline -- fitted
%%% from each run
fIn = stats.fStat;
for r = 1:size(fRes.fIn,1)
    dOut = strsplit(fRes.fIn{r},filesep); dOut = strjoin(dOut(1:end-1),filesep);
    fOut = strsplit(fIn        ,filesep); fOut = fOut{end};
    fOut = replace(fOut,'_stats.nii.gz','_poly0base.nii.gz');
    fOut = fullfile(dOut,fOut);
    stats.fPoly0Base{r,1} = fOut;
    if force || ~exist(fOut,'file')    
        buck = ['Run#' num2str(r) 'Pol#0_Coef'];
        cmd{end+1} = '3dbucket -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = [fIn '[' buck ']'];
    end
end
%%% concatenate across runs
if fRes.r==0 % only for analysis on catenated runs
    fIn = stats.fPoly0Base;
    fOut = strsplit(fIn{1}     ,'_'); fOut{contains(fOut,'run-')} = 'run-cat'; fOut = strjoin(fOut,'_');
    stats.fPoly0Base_cat = fOut;
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTcat -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = strjoin(fIn,' ');
    end
end
%%% average across runs
if fRes.r==0 % only for analysis on catenated runs
    fIn = stats.fPoly0Base_cat;
    fOut = replace(fIn,'_poly0base.nii.gz','_poly0baseAv.nii.gz');
    stats.fPoly0Base_catAv = fOut;
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTstat -overwrite -mean \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = fIn;
    end
end

%% Extract baseline -- temporal average
%%% from each run
for r = 1:size(fRes.fIn,1)
    fIn  = char(fRes.fIn(r));
    fOut = replace(fIn,'preproc_volTs.nii.gz','preproc_volTsAv.nii.gz');
    stats.fTsAvBase{r,1} = fOut;
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTstat -overwrite -mean \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = fIn;
    end
end
%%% concatenate across runs
if fRes.r==0 % only for analysis on catenated runs
    fIn = stats.fTsAvBase;
    fOut = strsplit(fIn{1},'_'); fOut{contains(fOut,'run-')} = 'run-cat'; fOut = strjoin(fOut,'_');
    stats.fTsAvBase_cat = fOut;
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTcat -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = strjoin(fIn,' ');
    end
end
if fRes.r==0 % only for analysis on catenated runs
    %%% average across runs
    fIn = stats.fTsAvBase_cat;
    fOut = replace(fIn,'preproc_volTsAv.nii.gz','av_preproc_volTsAv.nii.gz');
    stats.fTsAvBase_catAv = fOut;
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dTstat -overwrite -mean \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = fIn;
    end
end


for k = 0:fRes.param.dsgn.condK % 0 for the full model; >=1 for each event conditions
    if k>0
        %% Extract response ts and add baselines to it
        switch fRes.param.model
            case {'SPMG2'}
                stats.fResp{1,k}   = '';
                stats.fRespSd{1,k} = '';
            case {'TENT' 'TENTzero'}
                stats.fResp{1,k}   = fRes.fResp{1,k};
                stats.fRespSd{1,k} = fRes.fRespStd{1,k};
                fIn   = stats.fResp{1,k};
                % baseline from fit
                fOut  = replace(fIn,'_respAv.nii.gz','_respAvOnPoly0Base.nii.gz');
                if fRes.r==0
                    fBase = stats.fPoly0Base_catAv;
                else
                    fBase = char(stats.fPoly0Base);
                end
                stats.fRespOnPoly0Base{1,k} = fOut;
                if force || ~exist(fOut,'file')
                    cmd{end+1} = '3dcalc -overwrite \';
                    cmd{end+1} = ['-prefix ' fOut ' \'];
                    cmd{end+1} = ['-a ' fIn   ' \'];
                    cmd{end+1} = ['-b ' fBase ' \'];
                    cmd{end+1} = '-expr ''a+b''';
                end
                % baseline from temporal average
                fOut  = replace(fIn,'_respAv.nii.gz','_respAvOnTsAvBase.nii.gz');
                if fRes.r==0
                    fBase = stats.fTsAvBase_catAv;
                else
                    fBase = char(stats.fTsAvBase);
                end
                stats.fRespOnTsAvBase{1,k} = fOut;
                if force || ~exist(fOut,'file')
                    cmd{end+1} = '3dcalc -overwrite \';
                    cmd{end+1} = ['-prefix ' fOut ' \'];
                    cmd{end+1} = ['-a ' fIn   ' \'];
                    cmd{end+1} = ['-b ' fBase ' \'];
                    cmd{end+1} = '-expr ''a+b''';
                end
        otherwise
            dbstack; error('code that');
        end
    end

    %% Extract F-value
    fIn = fRes.fStat;
    if k==0 % full-model
        fOut = replace(fIn,'_stats.nii.gz','_fVal.nii.gz');
        stats.fFullF = fOut;
    else    % individual conditions of the model
        fOut = replace(fIn,'_stats.nii.gz','_fVal.nii.gz');
        fOut = replace(fOut,'cond-FULL',['cond-' stats.condList{k}]);
        stats.fCondF{k,1} = fOut;
    end
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dbucket -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        if k==0 % full-model
            cmd{end+1} = [fIn '[Full_Fstat]'];
        else    % individual conditions of the model
            cmd{end+1} = [fIn '[' stats.task '_' stats.condList{k} '_Fstat]'];
        end
    end

    %% Compute p-value
    if k==0 % full-model
        fIn  = stats.fFullF;
        fOut = replace(fIn,'_fVal.nii.gz','_fValP.nii.gz');
        stats.fFullF_pVal = fOut;
    else    % individual conditions of the model
        fIn  = stats.fCondF{k,1};
        fOut = replace(fIn,'_fVal.nii.gz','_fValP.nii.gz');
        stats.fCondF_pVal{k,1} = fOut;
    end
    if force || ~exist(fOut,'file')
        cmd{end+1} = ['df=$(3dAttribute BRICK_STATAUX ' fIn ')'];
        cmd{end+1} = 'df1=$(echo $df | awk ''{print $(NF-1)}'')';
        cmd{end+1} = 'df2=$(echo $df | awk ''{print $NF}'')';
        cmd{end+1} = '3dcalc -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = ['-a ' fIn ' \'];
        cmd{end+1} = '-expr "1-stat2cdf(a,4,$df1,$df2,0)" 2> /dev/null';
    end

    %% Compute q-value (fdr)
    if k==0 % full-model
        fIn  = stats.fFullF;
        fOut = replace(fIn,'_fVal.nii.gz','_fValQ.nii.gz');
        stats.fFullF_qVal      = fOut;
    else    % individual conditions of the model
        fIn  = stats.fCondF{k,1};
        fOut = replace(fIn,'_fVal.nii.gz','_fValQ.nii.gz');
        stats.fCondF_qVal{k,1} = fOut;
    end
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dFDR -overwrite -qval \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = ['-input ' fIn ' \'];
        cmd{end+1} = ['-mask ' fMask];
    end

    %% Transform coefficients
    switch fRes.param.model
        case {'SPMG2'}
            % dbstack; error('code that')
            % fIn = fRes.fStat;
            % fOut = replace(fIn,'_stats.nii.gz','_coef.nii.gz');
            % stats.fCoef = fOut;
            % if force || ~exist(fOut,'file')
            %     cmd{end+1} = '3dbucket -overwrite \';
            %     cmd{end+1} = ['-prefix ' fOut ' \'];
            %     cmd{end+1} = [fIn '[' param.funDsgn.label '#0_Coef,' param.funDsgn.label '#1_Coef]'];
            % end
        case {'TENT' 'TENTzero'}
        otherwise
            dbstack; error('code that');
    end
end



%% Run system commands
if length(cmd)>1
    if verbose
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end
    disp(' done')
else
    disp(' already done, skipping')
end
