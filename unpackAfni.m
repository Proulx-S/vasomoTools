function fRes = unpackAfni(fRes,fMask,force,verbose)
global src
if ~exist('fMask','var');     fMask = []; end
if ~exist('force','var');     force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force);     force = 0; end
if isempty(verbose); verbose = 0; end        
if numel(fRes)>1
    for i = 1:numel(fRes)
        fRes(i).stats = unpackAfni(fRes(i),fMask,force,verbose);
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
    if fRes.param.PCflag
        dOut = strsplit(fRes.fIn{r,1,1},filesep); dOut = strjoin(dOut(1:end-1),filesep);
        dOut = replace(dOut,'part-real','part-realImag');
    else
        dOut = strsplit(fRes.fIn{r},filesep); dOut = strjoin(dOut(1:end-1),filesep);
    end
    if fRes.r==0
        rStr = strsplit(dOut,'_'); rStr = rStr{contains(rStr,'run-')};
        dOut = strsplit(dOut,'_'); dOut{contains(dOut,'run-')} = 'run-cat'; dOut = strjoin(dOut,'_');
        dOut = strsplit(dOut,'_'); dOut(contains(dOut,'chunk-')) = [];      dOut = strjoin(dOut,'_');
    end
    if ~exist(dOut,'dir'); mkdir(dOut); end
    fOut = strsplit(char(fIn)        ,filesep); fOut = fOut{end};
    if fRes.r==0
        fOut = replace(fOut,'_stats',[ '_' rStr '_poly0base.nii.gz']);
    else
        fOut = replace(fOut,'_stats','_poly0base.nii.gz');
    end
    fOut = fullfile(dOut,fOut);
    if fRes.param.PCflag
        % stats.fPoly0Base{r,1} = replace(fOut,'.nii.gz'      ,'Mag.nii.gz'  );
        % stats.fPoly0Base{r,2} = replace(fOut,'.nii.gz'      ,'Phase.nii.gz');
        stats.fPoly0Base{r,1} = replace(fOut,'part-realImag','part-real'   );
        stats.fPoly0Base{r,2} = replace(fOut,'part-realImag','part-imag'   );
        if force || ...
            ~exist(stats.fPoly0Base{r,1},'file') || ...
            ~exist(stats.fPoly0Base{r,2},'file')
            % real
            buck = ['Run#' num2str(r) 'Pol#0_Coef'];
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' stats.fPoly0Base{r,1} ' \'];
            cmd{end+1} = [char(fIn) '+orig[' buck ']'];
            % imaginary
            buck = ['Run#' num2str(r+size(fRes.fIn,1)) 'Pol#0_Coef'];
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' stats.fPoly0Base{r,2} ' \'];
            cmd{end+1} = [char(fIn) '+orig[' buck ']'];
        end
    else
        stats.fPoly0Base{r,1} = fOut;
        if force || ~exist(fOut,'file')
            buck = ['Run#' num2str(r) 'Pol#0_Coef'];
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [char(fIn) '+orig[' buck ']'];
        end
    end
end
%%% concatenate across runs
if fRes.r==0 % only for analysis on catenated runs
    for v = 1:size(stats.fPoly0Base,2)
        fIn = stats.fPoly0Base(:,v);
        fOut = strsplit(fileparts(fIn{1})     ,'_'); fOut{contains(fOut,'run-')} = 'run-cat'; fOut = strjoin(fOut,'_');
        tmp = strsplit(fIn{1},filesep); tmp = tmp{end};
        tmp = strsplit(tmp,'_'); tmp{contains(tmp,'run-')} = 'run-cat'; tmp = strjoin(tmp,'_');
        fOut = fullfile(fOut,tmp);
        stats.fPoly0Base_cat(:,v) = cellstr(fOut);
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dTcat -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = strjoin(fIn,' ');
        end
    end
end
%%% average across runs
if fRes.r==0 % only for analysis on catenated runs
    for v = 1:size(stats.fPoly0Base_cat,2)
        fIn = stats.fPoly0Base_cat(:,v);
        fOut = strsplit(char(fIn),filesep); fOut{end} = replace(fOut{end},'run-cat','run-catAv'); fOut = strjoin(fOut,filesep);
        stats.fPoly0Base_catAv(:,v) = cellstr(fOut);
        if force || ~exist(char(fOut),'file')
            cmd{end+1} = '3dTstat -overwrite -mean \';
            cmd{end+1} = ['-prefix ' char(fOut) ' \'];
            cmd{end+1} = char(fIn);
        end
    end
end



if fRes.param.PCflag
    %%% transform individual runs
    for r = 1:size(stats.fPoly0Base,1)
        % phase
        stats.fPoly0Base{r,3} = replace(replace(stats.fPoly0Base{r,1},'part-real','part-realImag'),'_poly0base.nii.gz','_poly0basePhase.nii.gz');
        if force || ~exist(stats.fPoly0Base{r,3},'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' stats.fPoly0Base{r,3} ' \'];
            cmd{end+1} = ['-a '      stats.fPoly0Base{r,1} ' \'];
            cmd{end+1} = ['-b '      stats.fPoly0Base{r,2} ' \'];
            cmd{end+1} = '-expr ''atan2(b,a)''';
        end
        % abs
        stats.fPoly0Base{r,4} = replace(replace(stats.fPoly0Base{r,2},'part-imag','part-realImag'),'_poly0base.nii.gz','_poly0baseMag.nii.gz');
        if force || ~exist(stats.fPoly0Base{r,4},'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' stats.fPoly0Base{1,4} ' \'];
            cmd{end+1} = ['-a '      stats.fPoly0Base{r,1} ' \'];
            cmd{end+1} = ['-b '      stats.fPoly0Base{r,2} ' \'];
            cmd{end+1} = '-expr ''sqrt(a*a+b*b)''';
        end
    end
    if fRes.r==0
        %%% transform catenated runs
        % phase
        stats.fPoly0Base_cat{1,3} = replace(replace(stats.fPoly0Base_cat{1,1},'part-real','part-realImag'),'_poly0base.nii.gz','_poly0basePhase.nii.gz');
        if force || ~exist(stats.fPoly0Base_cat{1,3},'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' stats.fPoly0Base_cat{1,3} ' \'];
            cmd{end+1} = ['-a '      stats.fPoly0Base_cat{1,1} ' \'];
            cmd{end+1} = ['-b '      stats.fPoly0Base_cat{1,2} ' \'];
            cmd{end+1} = '-expr ''atan2(b,a)''';
        end
        % abs
        stats.fPoly0Base_cat{1,4} = replace(replace(stats.fPoly0Base_cat{1,2},'part-imag','part-realImag'),'_poly0base.nii.gz','_poly0baseMag.nii.gz');
        if force || ~exist(stats.fPoly0Base_cat{1,4},'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' stats.fPoly0Base_cat{1,4} ' \'];
            cmd{end+1} = ['-a '      stats.fPoly0Base_cat{1,1} ' \'];
            cmd{end+1} = ['-b '      stats.fPoly0Base_cat{1,2} ' \'];
            cmd{end+1} = '-expr ''sqrt(a*a+b*b)''';
        end
        %%% transform averaged runs
        % phase
        stats.fPoly0Base_catAv{1,3} = replace(replace(stats.fPoly0Base_catAv{1,1},'part-real','part-realImag'),'_poly0base.nii.gz','_poly0basePhase.nii.gz');
        if force || ~exist(stats.fPoly0Base_catAv{1,3},'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' stats.fPoly0Base_catAv{1,3} ' \'];
            cmd{end+1} = ['-a '      stats.fPoly0Base_catAv{1,1} ' \'];
            cmd{end+1} = ['-b '      stats.fPoly0Base_catAv{1,2} ' \'];
            cmd{end+1} = '-expr ''atan2(b,a)''';
        end
        % abs
        stats.fPoly0Base_catAv{1,4} = replace(replace(stats.fPoly0Base_catAv{1,2},'part-imag','part-realImag'),'_poly0base.nii.gz','_poly0baseMag.nii.gz');
        if force || ~exist(stats.fPoly0Base_catAv{1,4},'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' stats.fPoly0Base_catAv{1,4} ' \'];
            cmd{end+1} = ['-a '      stats.fPoly0Base_catAv{1,1} ' \'];
            cmd{end+1} = ['-b '      stats.fPoly0Base_catAv{1,2} ' \'];
            cmd{end+1} = '-expr ''sqrt(a*a+b*b)''';
        end
    end
end


tmp = strsplit(fRes.fIn{1,1},'_');
if ~any(contains(tmp,'part-')) || ~contains(tmp{contains(tmp,'part-')},'Mag1')
    %% Extract baseline -- temporal average
    %%% from each run
    for r = 1:size(fRes.fIn,1)
        fIn  = char(fRes.fIn(r));
        if fRes.param.PCflag
            fIn = replace(fIn,'part-real','part-mag');   fIn = strsplit(fIn,'_');
            fIn{contains(fIn,'rec-venc')} = 'rec-venc0'; fIn = strjoin(fIn,'_');
        end
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
end







for k = 0:fRes.param.dsgn.condK % 0 for the full model; >=1 for each event conditions
    if k>0
        %% Extract response ts and add baselines to it
        switch fRes.param.model
            case {'SPMG2'}
                stats.fResp{1,k}   = '';
                stats.fRespSd{1,k} = '';
            case {'TENT' 'TENTzero'}
                if fRes.param.PCflag
                    stats.fResp(1,k,:)   = fRes.fResp(1,k,:);
                    stats.fResp(1,k,3)   = replace(replace(stats.fResp(1,k,1),'part-real','part-realImag'),'_respAv.nii.gz','_respAvPhase.nii.gz');
                    stats.fResp(1,k,4)   = replace(replace(stats.fResp(1,k,1),'part-real','part-realImag'),'_respAv.nii.gz','_respAvMag.nii.gz');
                    
                    respR = stats.fResp{1,k,1};
                    respI = stats.fResp{1,k,2};
                    if fRes.r==0
                        baseR = stats.fPoly0Base_catAv{1,1};
                        baseI = stats.fPoly0Base_catAv{1,2};
                    else
                        baseR = stats.fPoly0Base{1,1};
                        baseI = stats.fPoly0Base{1,2};
                    end
                    
                    % get phase
                    if force || ~exist(stats.fResp{1,k,3},'file')
                        cmd{end+1} = '3dcalc -overwrite \';
                        cmd{end+1} = ['-prefix ' stats.fResp{1,k,3} ' \'];
                        cmd{end+1} = ['-a '      respR              ' \'];
                        cmd{end+1} = ['-b '      respI              ' \'];
                        cmd{end+1} = ['-c '      baseR              ' \'];
                        cmd{end+1} = ['-d '      baseI              ' \'];
                        cmd{end+1} = '-expr ''atan2(b+d,a+c)''';
                    end
                    % get mag
                    if force || ~exist(stats.fResp{1,k,4},'file')
                        cmd{end+1} = '3dcalc -overwrite \';
                        cmd{end+1} = ['-prefix ' stats.fResp{1,k,4} ' \'];
                        cmd{end+1} = ['-a '      respR              ' \'];
                        cmd{end+1} = ['-b '      respI              ' \'];
                        cmd{end+1} = ['-c '      baseR              ' \'];
                        cmd{end+1} = ['-d '      baseI              ' \'];
                        cmd{end+1} = '-expr ''sqrt((a+c)*(a+c)+(b+d)*(b+d))''';
                    end

                    % stats.fRespSd(1,k,:) = fRes.fRespStd(1,k,:);
                    % stats.fRespSd(1,k,3) = replace(replace(stats.fRespSd(1,k,1),'part-real','part-realImag'),'_respAv.nii.gz','_respAvPhase.nii.gz');
                    % stats.fRespSd(1,k,4) = replace(replace(stats.fRespSd(1,k,1),'part-real','part-realImag'),'_respAv.nii.gz','_respAvMag.nii.gz');
                    % % get phase
                    % if force || ~exist(stats.fRespSd{1,k,3},'file')
                    %     cmd{end+1} = '3dcalc -overwrite \';
                    %     cmd{end+1} = ['-prefix ' stats.fRespSd{1,k,3} ' \'];
                    %     cmd{end+1} = ['-a '      stats.fRespSd{1,k,1} ' \'];
                    %     cmd{end+1} = ['-b '      stats.fRespSd{1,k,2} ' \'];
                    %     cmd{end+1} = '-expr ''atan2(b,a)''';
                    % end
                    % % get mag
                    % if force || ~exist(stats.fRespSd{1,k,4},'file')
                    %     cmd{end+1} = '3dcalc -overwrite \';
                    %     cmd{end+1} = ['-prefix ' stats.fRespSd{1,k,4} ' \'];
                    %     cmd{end+1} = ['-a '      stats.fRespSd{1,k,1} ' \'];
                    %     cmd{end+1} = ['-b '      stats.fRespSd{1,k,2} ' \'];
                    %     cmd{end+1} = '-expr ''sqrt(a*a+b*b)''';
                    % end
                    % 
                    % % too lazy to add baseline to mag
                    % stats.fRespOnPoly0Base = {};
                    % stats.fRespOnTsAvBase = {};
                else
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
                        cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                        cmd{end+1} = ['-a ' char(fIn)   ' \'];
                        cmd{end+1} = ['-b ' char(fBase) ' \'];
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
                        cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                        cmd{end+1} = ['-a ' char(fIn)   ' \'];
                        cmd{end+1} = ['-b ' char(fBase) ' \'];
                        cmd{end+1} = '-expr ''a+b''';
                    end
                end
        otherwise
            dbstack; error('code that');
        end
        if fRes.param.PCflag
            % Need a way to get stats for the combined the real and imaginay regressors
            % 3dDeconvolve glt?
            % Right now we just rely on the full model, which works as long as there is only 2 regressors, one real and one imaginary
            continue
        end
    end

    %% Extract F-value
    fIn = char(stats.fStat);
    if k==0 % full-model
        fOut = replace(fIn,'_stats','_fVal.nii.gz');
        stats.fFullF = fOut;
    else    % individual conditions of the model
        fOut = replace(fIn,'_stats','_fVal.nii.gz');
        fOut = replace(fOut,'cond-FULL',['cond-' stats.condList{k}]);
        stats.fCondF{1,k} = char(fOut);
    end
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dbucket -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        if k==0 % full-model
            cmd{end+1} = [char(fIn) '+orig[Full_Fstat]'];
        else    % individual conditions of the model
            cmd{end+1} = [char(fIn) '+orig[' stats.task '_' stats.condList{k} '_Fstat]'];
        end
    end

    %% Compute p-value
    if k==0 % full-model
        fIn  = char(stats.fFullF);
        fOut = replace(fIn,'_fVal.nii.gz','_fValP.nii.gz');
        stats.fFullF_pVal = fOut;
    else    % individual conditions of the model
        fIn  = stats.fCondF{1,k};
        fOut = replace(fIn,'_fVal.nii.gz','_fValP.nii.gz');
        stats.fCondF_pVal{1,k} = fOut;
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
        fIn  = char(stats.fFullF);
        fOut = replace(fIn,'_fVal.nii.gz','_fValQ.nii.gz');
        stats.fFullF_qVal      = fOut;
    else    % individual conditions of the model
        fIn  = stats.fCondF{1,k};
        fOut = replace(fIn,'_fVal.nii.gz','_fValQ.nii.gz');
        stats.fCondF_qVal{1,k} = fOut;
    end
    if force || ~exist(fOut,'file')
        cmd{end+1} = '3dFDR -overwrite -qval \';
        cmd{end+1} = ['-prefix ' char(fOut) ' \'];
        cmd{end+1} = ['-input ' char(fIn) ' \'];
        cmd{end+1} = ['-mask ' char(fMask)];
    end
end



%% Run system commands
if length(cmd)>1
    if verbose
        [status,cmdout] = system(strjoin(cmd,newline),'-echo');
    else
        [status,cmdout] = system(strjoin(cmd,newline));
    end
    if status || isempty(cmdout) || contains(cmdout,'ERROR','IgnoreCase',false); dbstack; error(cmdout); end
    disp(' done')
else
    disp(' already done, skipping')
end



%% Transform coefficients -- from the double gamma fit, find the main vector (principal response delay across voxels) and adjust accordingly
for k = 1:fRes.param.dsgn.condK
    if (fRes.r==0 || fRes.R==1) && strcmp(fRes.param.model,'SPMG2') % only for analysis on catenated runs for sufficient precision in delay estimation
        
        % Extract coefficients
        fIn = char(stats.fStat);
        fOut = replace(fIn,'cond-FULL',['cond-' stats.condList{k}]);
        fOut = char(replace(fOut,'_stats','_coefs.nii.gz'));
        stats.fCondCoef{1,k} = fOut;
        if force || ~exist(fOut,'file')
            cmd = {src.afni};
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '+orig[' stats.task '_' stats.condList{k} '#0_Coef,' stats.task '_' stats.condList{k} '#1_Coef]'];
            [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        end
        
        % Convert to complex values
        fIn = fOut;
        fFig = replace(fIn,'.nii.gz','MainVector.fig');        
        stats.fCondCoef_mainVector{1,k} = fFig;
        fOutPolar = replace(fIn,'_coefs.nii.gz','_polar.nii.gz');
        stats.fCondPolar{1,k} = fOutPolar;
        fOutCoef = replace(fIn,'_coefs.nii.gz','_coefsAdj.nii.gz');
        stats.fCondCoef_adj{1,k} = fOutCoef;
        fOutCoefFlag = replace(fOutCoef,'.nii.gz','.flag');
        stats.fCondCoef_adjFlag{1,k} = fOutCoefFlag;
        if force || ~exist(fOutPolar,'file') || ~exist(fOutCoef,'file') || ~exist(fOutCoefFlag,'file') || ~exist(fFig,'file')
            mriCoef = MRIread(fIn);
            mriCoef.vol = complex(mriCoef.vol(:,:,:,1), mriCoef.vol(:,:,:,2));
            
            % Get slope (main vector) of the data from significant voxels
            mriQ = MRIread(stats.fCondF_qVal{1,k});
            mask = mriQ.vol<0.05;
            if nnz(mask(:))>15
                coefAdjFlag = 1;
            else
                coefAdjFlag = 0;
                mriP = MRIread(stats.fCondF_pVal{1,k});
                mask = mriP.vol<0.05;
            end
            fid = fopen(fOutCoefFlag, 'w'); fprintf(fid, '%d', coefAdjFlag); fclose(fid);
            
            slp = real(mriCoef.vol(mask))\imag(mriCoef.vol(mask));
            v = complex(1,slp); v = v./abs(v);

            % Visualize the main response vector with individual voxels
            hFig = figure('Visible','off');
            scatter(real(mriCoef.vol(mask)),imag(mriCoef.vol(mask)),'.k'); hold on;
            lim = [-1 1].*max(abs([real(mriCoef.vol(mask)); imag(mriCoef.vol(mask))]));
            line(lim,lim.*slp,'Color','r');
            xlim(lim); ylim(lim); grid on; ax = gca; ax.DataAspectRatio = [1 1 1];
            xlabel('Real'); ylabel('Imaginary');
            legend('Voxels','Main Vector');
            [a,b,~] = fileparts(replace(fIn,'.nii.gz',''));
            [~,a,~] = fileparts(a);
            if ~coefAdjFlag
                title([a newline b newline '!!!WARNING!!! showing p<0.05 voxels, delay not corrected (less than 15 q<0.05 voxels)'],'Interpreter','none');
            else
                title([a newline b],'Interpreter','none');
            end
            % and save
            set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
            savefig(hFig, fFig, 'compact');
            if verbose > 0
                hFig.Visible = 'on';
                hFig.WindowStyle = 'docked';
                drawnow
            else
                close(hFig);
            end

            % Adjusting relative to principal vector
            rho   = abs(mriCoef.vol);                         % Magnitude
            if ~coefAdjFlag
                theta = angle(mriCoef.vol);
            else
                theta = wrapToPi(angle(mriCoef.vol) - angle(v));  % Phase
            end
            mriPol = mriCoef; mriPol.fspec = fOutPolar;
            mriPol.vol(:,:,:,1) = rho;                         % Magnitude
            mriPol.vol(:,:,:,2) = theta;                       % Phase
            MRIwrite(mriPol,fOutPolar);

            mriCoef.fspec = fOutCoef;
            [mriCoef.vol(:,:,:,1),mriCoef.vol(:,:,:,2)] = pol2cart(theta,rho);
            MRIwrite(mriCoef,fOutCoef);
        else
            coefAdjFlag = readmatrix(fOutCoefFlag,'FileType','text');
        end
        if ~coefAdjFlag
            % Erase fCondCoef_adj because fCondCoef was not adjusted
            stats.fCondCoef_adj{1,k} = '';
        end
    else
        stats.fCondCoef{1,k}            = '';
        stats.fCondCoef_mainVector{1,k} = '';
        stats.fCondPolar{1,k}           = '';
        stats.fCondCoef_adj{1,k}        = '';
    end
end

%% Output
fRes.stats = stats;