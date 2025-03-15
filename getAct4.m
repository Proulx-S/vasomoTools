function [fRun,fSes,fSes_echoCat,param] = getAct4(volTs,dsgn,fMask,param,force,verbose)
% see /autofs/space/takoyaki_001/users/proulxs/tools/vasomoTools/getResp2.m
fSes_echoCat = [];
% global srcAfni srcFs
if ~exist('force','var');     force = []; end
if ~exist('verbose','var'); verbose = []; end
if ~exist('dsgn','var');       dsgn = []; end
if ~exist('fMask','var');   fMask = []; end
if isempty(force);     force = 0; end
if isempty(verbose); verbose = 0; end
if isempty(dsgn)
    if isfield(volTs,'dsgn') && isequal(dsgn)
        dsgn = dsgn;
    else
        error('badly specified dsgn')
    end
end
if ~isfield(param,'skipMov'); param.skipMov = []; end
if ~isfield(param,'skipCat'); param.skipCat = []; end
if ~isfield(param,'skipRun'); param.skipRun = []; end
if ~isfield(param,'dryRun'); param.dryRun = []; end
if isempty(param.skipMov); param.skipMov = 0; end
if isempty(param.skipCat); param.skipCat = 0; end
if isempty(param.skipRun); param.skipRun = 0; end
if isempty(param.dryRun); param.dryRun = 0; end

%% Data files
if ~isempty(volTs)
    fVolTs = {volTs.fspec}';
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
else
    fMask;
end
if ~all(diff([volTs.tr])<0.01); dbstack; error('runs have different tr'); end

%% Functional design
if ~isa(dsgn,'runDsgn')
    dbstack; error('old version, reconciliate')
end
if ~isempty(dsgn.cond)
    k = sort(unique(dsgn.cond)); if any(diff(k)-1); dbstack; error('cond indices are not monotonically increasing'); end
    dsgn.condLabel = repmat({'stim'},size(k));
    dsgn.condLabel{k==0} = 'catch';
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
param.tr = [volTs.tr]./1000;
param.dsgn = dsgn;


% param.tr = volTs(1).tr/1000;
% param.funDsgn.trStim    = dsgn.dt;
% param.funDsgn.k         = 1;
% param.funDsgn.startSeq  = dsgn.onsetList;
% param.funDsgn.durSeq    = dsgn.ondurList;
% if isa(dsgn,'runDsgn') && ~isempty(dsgn.cond)
%     k = unique(dsgn.cond); if any(diff(k)-1); dbstack; error('cond indices are not monotonically increasing'); end
%     dsgn.condLabel = repmat({'stim'},size(k));
%     dsgn.condLabel{k==0} = 'catch';
%     dsgn.condK = length(k);
% else
%     param.funDsgn.condSeq   = ones(size(param.funDsgn.startSeq));
%     param.funDsgn.condLabel = {'stim'};
%     if isfield(dsgn,'nullTrial') && ~isempty(dsgn.nullTrial)
%         param.funDsgn.condSeq(dsgn.nullTrial) = 2;
%         param.funDsgn.condLabel{end+1} = 'catch';
%         param.funDsgn.k = 2;
%     end
% end
% param.funDsgn.label    = dsgn.task;
% param.funDsgn.trDecon  = dsgn.trDecon;


% %% Run afni's 3dDeconvolve for response timecourse estimation
% [fRun,fSes,param] = runAfni(fVolTs,param,fMask,force,verbose); % analysis performed on each echoe within that function
% nEcho = size(fRun,2);
% nRun = size(fRun,1);

%% Run afni's 3dDeconvolve for double-gamma response amplitude (and delay) fit
%  'SPMG1'       = 1 parameter SPM gamma variate basis function
%      exp(-t)*(A1*t^P1-A2*t^P2) where
%    A1 = 0.0083333333  P1 = 5  (main positive lobe)
%    A2 = 1.274527e-13  P2 = 15 (undershoot part)
%    This function is NOT normalized to have peak=1!
% 'SPMG2'       = 2 parameter SPM: gamma variate + d/dt derivative
%    [For backward compatibility: 'SPMG' == 'SPMG2']
%  'SPMG3'       = 3 parameter SPM basis function set
%            ==> ** The SPMGx functions now can take an optional
%                   (duration) argument, specifying that the primal
%                    SPM basis functions should be convolved with
%                     a square wave 'duration' seconds long and then
%                  be normalized to have peak absolute value = 1;
%                  e.g., 'SPMG3(20)' for a 20 second duration with
%                  three basis function.  [28 Apr 2009]
%               ** Note that 'SPMG1(0)' will produce the usual
%                     'SPMG1' wavefunction shape, but normalized to
%                   have peak value = 1 (for example).
if ~isfield(param,'model') || isempty(param.model); param.model = 'SPMG2'; end
[fRun,fSes,param] = runAfni(fVolTs,param,fMask,force,verbose); % analysis performed on each echoe within that function

%% Get model
verboseThis = verbose;
if ~isempty(fRun)
    for R = 1:size(fRun,1)
        paramCur = param; paramCur.nFrame = paramCur.nFrame(R); paramCur.tr = paramCur.tr(R);
        fRun(R,1) = plotDsgnMat(fRun(R,1),paramCur,fVolTs(R,1),volTs(R,1),verboseThis);
    end
end
if ~isempty(fSes)
    dbstack; error('double-check that')
    fSes = plotDsgnMat(fSes,param,fVolTs,volTs,verboseThis);
end

%% Refactor
if ~isempty(fRun)
    for R = 1:size(fRun,1)
        tmp(R,1).afni = fRun(R,1);
        if isfield(fRun,'fResp')
            % tmp(R,1).afni = rmfield(fRun(R,1),'fResp');
            tmp(R,1).fs.fRespTs = fRun(R,1).fResp;
        end
        if isfield(fRun,'fRespStd')
            tmp(R,1).fs.fRespTsTrialSd = fRun(R,1).fRespStd;
        end
        tmp(R,1).fs.fMask = fRun(R,1).fMask;
    end
    fRun = tmp; clear tmp
end
if ~isempty(fSes)
    dbstack; error('double-check that')
    if isfield(fSes,'fResp')
        tmp.afni = rmfield(fSes,'fResp');
        tmp.fs.fRespTs = fSes.fResp;
    else
        tmp.afni = fSes;
    end
    tmp.fs.fMask   = fMask;
    fSes = tmp; clear tmp
end
nEcho = size(fVolTs,2);
nRun = size(fVolTs,1);


if ~isempty(fRun)
    f = fRun;
else
    f = [];
end
if ~isempty(fSes)
    dbstack; error('double-check that')
    f = [f; fSes];
end
% if ~isempty(fRun) && ~isempty(fSes)
%     f = [fRun; fSes];
% else
%     f = [];
%     dbstack; error('X');
% end


% 
% 
% 
% if ~isempty(fSes)
%     fSes = plotDsgnMat(fSes,param,fVolTs,volTs,verboseThis);
%     if isfield(fSes,'fResp')
%         tmp.afni = rmfield(fSes,'fResp');
%         if ~isempty(fSes.fResp)
%             tmp.fs.fRespTs = fSes.fResp;
%         end
%         fSes = tmp; clear tmp
%     else
%         tmp = fSes; clear fSes
%         fSes.afni = tmp;
%     end
%     fSes.fs.fMask = fMask;
% end
% nEcho = size(fVolTs,2);
% nRun = size(fVolTs,1);


%% Unpack outputs
disp('Unpacking outputs')
forceThis = force;
if ~isempty(f)
    cmd = {};
    if exist('srcAfni','var') && ~isempty(srcAfni); cmd{end+1} = srcAfni; end
    for i = 1:numel(f)

        %%% Extract baseline -- fitted
        fIn = f(i).afni.fStat;
        fOut = replace(fIn,'_stats.nii.gz','_fBasePoly0.nii.gz');
        f(i).fs.fBasePoly0 = fOut;
        if forceThis || ~exist(fOut,'file')
            buck = num2str(1:size(f(i).afni.fIn,1),'Run#%iPol#0_Coef,'); buck(end) = [];
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[' buck ']'];
        end

        %%% Extract baseline -- temporal average
        fIn  = char(f(i).afni.fIn);
        fOut = replace(fIn,'preproc_volTs.nii.gz','av_preproc_volTs.nii.gz');
        f(i).fs.fBaseTsAv = fOut;
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dTstat -overwrite -mean \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = fIn;
        end


        for k = 0:param.dsgn.condK % 0 for the full model; >=1 for each event conditions

            %%% Add baselines to response ts
            if k>0
                fIn   = f(i).fs.fRespTs{k};
                %%%% baseline from fit
                fOut  = replace(fIn,'_respAv.nii.gz','_respAvOnBasePoly0.nii.gz');
                fBase = f(i).fs.fBasePoly0;
                f(i).fs.fRespTsOnBasePoly0{k,1} = fOut;
                if force || ~exist(fOut,'file')
                    cmd{end+1} = '3dcalc -overwrite \';
                    cmd{end+1} = ['-prefix ' fOut ' \'];
                    cmd{end+1} = ['-a ' fIn   ' \'];
                    cmd{end+1} = ['-b ' fBase ' \'];
                    cmd{end+1} = '-expr ''a+b''';
                end
                %%%% baseline from temporal average
                fOut = replace(fIn,'_respAv.nii.gz','_respAvOnBaseTsAv.nii.gz');
                fBase = f(i).fs.fBaseTsAv;
                f(i).fs.fRespTsOnBaseTsAv{k,1} = fOut;
                if force || ~exist(fOut,'file')
                    cmd{end+1} = '3dcalc -overwrite \';
                    cmd{end+1} = ['-prefix ' fOut ' \'];
                    cmd{end+1} = ['-a ' fIn   ' \'];
                    cmd{end+1} = ['-b ' fBase ' \'];
                    cmd{end+1} = '-expr ''a+b''';
                end
            end

            %%% Extract F-value
            fIn = f(i).afni.fStat;
            if k==0 % full-model
                fOut = replace(fIn,'_stats.nii.gz','_fVal.nii.gz');
                f(i).fs.fFullF = fOut;
            else    % individual conditions of the model
                fOut = replace(fIn,'_stats.nii.gz','_fVal.nii.gz');
                fOut = replace(fOut,'cond-FULL',['cond-' param.dsgn.condLabel{k}]);
                f(i).fs.fCondF{k,1} = fOut;
            end
            if forceThis || ~exist(fOut,'file')
                cmd{end+1} = '3dbucket -overwrite \';
                cmd{end+1} = ['-prefix ' fOut ' \'];
                if k==0 % full-model
                    cmd{end+1} = [fIn '[Full_Fstat]'];
                else    % individual conditions of the model
                    cmd{end+1} = [fIn '[' param.dsgn.task '_' param.dsgn.condLabel{k} '_Fstat]'];
                end
            end

            %%% Compute p-value
            if k==0 % full-model
                fIn  = f(i).fs.fFullF;
                fOut = replace(fIn,'_fVal.nii.gz','_fValP.nii.gz');
                f(i).fs.fFullF_pVal      = fOut;
            else    % individual conditions of the model
                fIn  = f(i).fs.fCondF{k,1};
                fOut = replace(fIn,'_fVal.nii.gz','_fValP.nii.gz');
                f(i).fs.fCondF_pVal{k,1} = fOut;
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

            %%% Compute q-value (fdr)
            if k==0 % full-model
                fIn  = f(i).fs.fFullF;
                fOut = replace(fIn,'_fVal.nii.gz','_fValQ.nii.gz');
                f(i).fs.fFullF_qVal      = fOut;
            else    % individual conditions of the model
                fIn  = f(i).fs.fCondF{k,1};
                fOut = replace(fIn,'_fVal.nii.gz','_fValQ.nii.gz');
                f(i).fs.fCondF_qVal{k,1} = fOut;
            end


            %%% Coef
            switch param.model
                case {'SPMG2'}
                    dbstack; error('code that')
                    fIn = f(i).afni.fStat;
                    fOut = replace(fIn,'_stats.nii.gz','_coef.nii.gz');
                    f(i).fs.fCoef = fOut;
                    if forceThis || ~exist(fOut,'file')
                        cmd{end+1} = '3dbucket -overwrite \';
                        cmd{end+1} = ['-prefix ' fOut ' \'];
                        cmd{end+1} = [fIn '[' param.funDsgn.label '#0_Coef,' param.funDsgn.label '#1_Coef]'];
                    end

                case {'SPMG3'}
                    dbstack; error('code that')
                case {'TENT' 'TENTzero'}
                otherwise
                    dbstack; error('code that');
            end

        end
    end
end


%%% Run system commands
if length(cmd)>1
    if verbose
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
        % [status,cmdout] = system(strjoin(cmd(1:4),newline)); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
    end
    disp(' done')
else
    disp(' already done, skipping')
end





return





% % % %%% Convert SPMG2 cartesian responses coefficient (gamma + first derivative) to polar (amplitude + delay) coefficient
% % % switch param.model
% % %     case {'SPMG2'}
% % %         dbstack; error('double-check that')
% % %         disp('convert hrf+derivative cartesian coefficients to polar coefficients')
% % % 
% % %         for i = 1:size(f,1)
% % %             fIn     = f(i).fs.fCoef;
% % %             fOut    = replace(fIn,'_coef.nii.gz','_coefPol.nii.gz');
% % %             fOutFig = replace(fIn,'_coef.nii.gz','_coefPol.fig');
% % %             fFDR    = f(i).fs.fFullQ;
% % %             f(i).fs.fCoefPol = fOut;
% % % 
% % %             verboseThis = verbose;
% % %             forceThis   = force;
% % %             if forceThis || ~exist(fOut,'file') || ~exist(fOutFig,'file')
% % %                 % hMat = figure('WindowStyle','docked');
% % %                 hMat = figure('Visible','off');
% % % 
% % %                 coef = MRIread(fIn);
% % %                 fdr  = MRIread(fFDR);
% % %                 mask = MRIread(fMask);
% % %                 mask = mask.vol & fdr.vol<0.05;
% % % 
% % %                 coef.vol = complex(coef.vol(:,:,:,1),coef.vol(:,:,:,2));
% % %                 scatter(real(coef.vol(mask)),imag(coef.vol(mask)));
% % %                 ax = gca; ax.DataAspectRatio = [1 1 1];
% % %                 grid on
% % %                 axis([-1 1 -1 1].*max(abs(axis)))
% % %                 xline(0,'k'); yline(0,'k');
% % %                 xlabel('SPM canon (coef)')
% % %                 ylabel('SPM canon derivative (coef)')
% % % 
% % %                 %get principal vector
% % %                 slp = real(coef.vol(mask))\imag(coef.vol(mask));
% % %                 hRef = refline(slp,0); hRef.Color = 'r';
% % %                 v = complex(1,slp); v = v./abs(v);
% % %                 title([num2str(angle(v)/pi*180,'%0.1f°') ' deviation from expected HR delay'])
% % % 
% % %                 % subtract that vector orientation from data
% % %                 coefPol = coef;
% % %                 coefPol.vol(:,:,:,1) = abs(coef.vol);
% % %                 coefPol.vol(:,:,:,2) = wrapToPi( angle(coef.vol) - angle(v) );
% % %                 MRIwrite(coefPol,fOut);
% % % 
% % %                 if verboseThis>1
% % %                     hMat.Visible = 'on';
% % %                     hMat.WindowStyle = 'docked';
% % %                     savefig(hMat,fOutFig,'compact')
% % %                 else
% % %                     set(hMat, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
% % %                     savefig(hMat,fOutFig,'compact')
% % %                     close(hMat)
% % %                 end
% % % 
% % %                 disp(' done')
% % %             else
% % %                 disp(' already done, skipping')
% % %             end
% % %         end
% % % 
% % %     case {'SPMG3'}
% % %         dbstack; error('code that')
% % %     case {'TENT' 'TENTzero'}
% % %     otherwise
% % %         dbstack; error('code that');
% % % end
% % % 
% % % 
% % % 
% % % 
% % % 
% % % if param.skipRun && ~param.skipCat
% % %     fRun = [];
% % %     fSes = f;
% % % elseif ~param.skipRun && ~param.skipCat
% % %     fRun = f(1:end-1);
% % %     fSes = f(end);
% % %     clear f
% % % else
% % %     dbstack; error('fix that mess')
% % %     if ~isempty(fRun) && ~isempty(fSes)
% % %         fRun = f(1:end-1,:);
% % %         fSes = f(end,:);
% % %         f    = [];
% % %     else
% % %         dbstack; error('X');
% % %     end
% % % end
% % % 
% % % %%% Run command
% % % disp('Writing baselines')
% % % if length(cmd)>1
% % %     if verbose
% % %         [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || isempty(cmdout) || contains(cmdout,'error','IgnoreCase',true); dbstack; error(cmdout); error('x'); end
% % %     else
% % %         [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout) || contains(cmdout,'error','IgnoreCase',true); dbstack; error(cmdout); error('x'); end
% % %     end
% % %     disp(' done')
% % % else
% % %     disp(' already done, skipping')
% % % end
% % % %% %%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%
%% Make videos %%
%%%%%%%%%%%%%%%%%
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







function [fRun,fSes,param] = runAfni(fList,param,fMask,force,verbose)
% global srcAfni srcFs
if ~exist('fMask','var'); fMask = []; end
if ~exist('force','var'); force = []; end
if ~exist('verbose','var'); verbose = []; end
if isempty(force); force = 0; end
if isempty(verbose); verbose = 0; end


%% %%%%%%%%%%%%%%%%%%
% Functional design %
%%%%%%%%%%%%%%%%%% %%
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


if size(fList,1)==1
    param.skipCat = 1;
end
if param.skipCat && param.skipRun
    param.skipRun = 0;
end


%% %%%%%
% Mask %
%%%%%%%%
mriMask = MRIread(fMask);
mriMask.vol([1:5 end-4:end],:              ) = 0;
mriMask.vol(:              ,[1:5 end-4:end]) = 0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed-effect model on individual runs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if param.skipRun
    fRun = [];
else
    cmd = {};
    if size(fList,2)>1; dbstack; error('double-check that'); end
    for E = 1:size(fList,2)
        for R = 1:size(fList,1)

            %% Define files
            fIn = fList(R,E);
            fOut = char(fIn);
            if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end

            fStat  = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats.nii.gz']);
            fFit   = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_fit.nii.gz'  ]);
            fResid = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_resid.nii.gz']);
            fMask  = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_mask.nii.gz' ]);
            mriMask.fspec = fMask; MRIwrite(mriMask,fMask);
            switch HRmodel
                case {'TENT' 'TENTzero'}
                    fResp    = cell(size(param.dsgn.condLabel));
                    fRespStd = cell(size(param.dsgn.condLabel));
                    for k = 1:param.dsgn.condK
                        fResp{k}    = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_respAv.nii.gz']);
                        fRespStd{k} = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_respSd.nii.gz']);
                    end
                case {'SPMG2' 'SPMG3'}
                otherwise
                    dbstack; error('figure that out')
            end
            if param.dryRun
                tmpName = tempname;
                fStim   = [tmpName '_startTime.1D'  ];
                fMat    = [tmpName '_stats.xmat.1D' ];
                fMatFig = [tmpName '_stats.xmat.fig'];
            else
                fStim   = cell(size(param.dsgn.condLabel));
                for k = 1:param.dsgn.condK
                    fStim{k} = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-' param.dsgn.condLabel{k} '_model-' HRmodel '_startTime.1D']);
                end
                fMat    = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats.xmat.1D']);
                fMatFig = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_cond-FULL_model-' HRmodel '_stats.xmat.fig']);
            end

            curParam = param;
            curParam.nFrame = param.nFrame(R);
            curParam.tr     = param.tr(R);
            


            fRun(R,E).fIn     = fIn;
            fRun(R,E).fFit    = fFit;
            fRun(R,E).fResid  = fResid;
            fRun(R,E).fStat   = fStat;
            fRun(R,E).fMask   = fMask;
            fRun(R,E).fMat    = fMat;
            fRun(R,E).fMatFig = fMatFig;
            switch HRmodel
                case {'TENT' 'TENTzero'}
                    fRun(R,E).fResp    = fResp;
                    fRun(R,E).fRespStd = fRespStd;
                case {'SPMG2' 'SPMG3'}
                otherwise
                    dbstack; error('figure that out')
            end



            %% Contruct afni command
            cmdTmp = {};
            if size(fList,2)>1
                cmdTmp{end+1} = ['echo ''  ''run' num2str(R) '/' num2str(size(fList,1)) ' -- echo' num2str(E) '/' num2str(size(fList,2))];
            else
                cmdTmp{end+1} = ['echo ''  ''run' num2str(R) '/' num2str(size(fList,1))];
            end

            switch HRmodel
                case {'TENT' 'TENTzero'}
                    [cmdTmpTmp,param.dsgn.nReg] = afniCmd(fIn,fStim,fMask,curParam,fResp,fRespStd,fFit,fResid,fMat,fStat,verbose,param.dryRun);
                    % [cmdTmpTmp,param.funDsgn.nReg] = afniCmd(fIn,fMask,fStim,param.nDummy,param.tr,startSeq,durSeq,condSeq,HRmodel,param.funDsgn.label,[],fResp,fFit,fResid,fMat,fStat,verbose,param.nDummyRemoved,param.trDecon,param.dryRun);
                case {'SPMG2' 'SPMG3'}
                    dbstack; error('double-check')
                    [cmdTmpTmp,param.funDsgn.nReg] = afniCmd(fIn,fMask,fStim,param.nDummy,param.tr,startSeq,durSeq,condSeq,HRmodel,param.funDsgn.label,[],[]   ,fFit,fResid,fMat,fStat,verbose,param.nDummyRemoved,[]           ,param.dryRun);
                otherwise
                    dbstack; error('figure that out')
            end
            if param.dryRun
                system(strjoin([{srcAfni} cmdTmpTmp],newline))
            end

            if ~exist(fStat,'file') || force
                cmdTmp = [cmdTmp cmdTmpTmp];
                if exist('fResp','var') && ~isempty(fResp); cmdTmp{end+1} = ['echo ''   ''' strjoin(cellstr(fResp),' ')]; end
                cmdTmp{end+1} = ['echo ''   ''' fStat];
                cmdTmp{end+1} = ['echo ''   ''' fMat];
            else
                if exist('fResp','var') && ~isempty(fResp); cmdTmp{end+1} = ['echo ''   ''' strjoin(cellstr(fResp),' ')]; end
                cmdTmp{end+1} = ['echo ''   ''' fStat];
                cmdTmp{end+1} = ['echo ''   ''' fMat];
                cmdTmp{end+1} = 'echo ''   ''already done, skipping';
            end

            fRun(R,E).cmd =  strjoin(cmdTmp,newline);

            cmd = [cmd cmdTmp];  
        end
    end

    %% Run afni command
    [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end


    % %% Get per-trial dsgn matrix
    % switch HRmodel
    %     case {'TENT' 'TENTzero'}
    %         condSeqX = 1:length(condSeq);
    %         for R = 1:size(fList,1)
    %             fIn = fList(R,E);
    %             fOut = char(fIn);
    %             fMatPerTrial = fullfile(fileparts(replace(fOut,'.nii.gz','')),['cond-visOn_model-' HRmodel '_stats.xmatPerTrial.1D']);
    %             [cmdXmat,paramX.funDsgn.nReg] = afniCmd(fIn,fMask,fStim,param.nDummy,param.tr,startSeq,durSeq,condSeqX,HRmodel,param.funDsgn.label,param,fResp,fFit,fResid,fMatPerTrial,fStat,verbose,param.nDummyRemoved,param.trDecon,1,param.nFrame(R));
    %             cmdXmat{end} = [cmdXmat{end} ' > /dev/null 2>&1']; cmdXmat = [{srcAfni} cmdXmat];
    %             if force || ~exist(fMatPerTrial,'file')
    %                 system(strjoin(cmdXmat,newline));
    %             end
    %             cmdX = {srcAfni}; cmdX{end+1} = ['1dcat ' fMatPerTrial];
    %             [~,cmdoutX] = system(strjoin(cmdX,newline));
    %             param.perTrialXmat(R).f = fMatPerTrial;
    %             param.perTrialXmat(R).mat  = str2num(cmdoutX);
    %             param.perTrialXmat(R).nReg = paramX.funDsgn.nReg;
    %         end
    %     case {'SPMG2' 'SPMG3'}
    %     otherwise
    %         dbstack; error('figure that out')
    % end


    
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed-effect model on concatenated runs (separate baselines) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if param.skipCat
    fSes = [];
else
    % dbstack; error('double-check that')
    cmd = {srcAfni};
    for E = 1:size(fList,2)

        %% Define files
        fIn = fList(:,E);
        tmp = strsplit(fIn{1},'_'); tmp = char(tmp(contains(tmp,'run-')));
        fOut = replace(fIn{1},tmp,'run-cat'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end

        fStat   = fullfile(fileparts(replace(fOut,'.nii.gz','')),['cond-visOn_model-' HRmodel '_stats.nii.gz']);
        fFit    = fullfile(fileparts(replace(fOut,'.nii.gz','')),['cond-visOn_model-' HRmodel '_fit.nii.gz']);
        fResid  = fullfile(fileparts(replace(fOut,'.nii.gz','')),['cond-visOn_model-' HRmodel '_resid.nii.gz']);
        switch HRmodel
            case {'TENT' 'TENTzero'}
                fResp   = fullfile(fileparts(replace(fOut,'.nii.gz','')),['cond-visOn_model-' HRmodel '_resp.nii.gz']);
            case {'SPMG2' 'SPMG3'}
            otherwise
                dbstack; error('figure that out')
        end
        if param.dryRun
            tmpName = tempname;
            fStim   = [tmpName '_startTime.1D'  ];
            fMat    = [tmpName '_stats.xmat.1D' ];
            fMatFig = [tmpName '_stats.xmat.fig'];
        else
            fStim   = fullfile(fileparts(replace(fOut,'.nii.gz','')),['cond-visOn_model-' HRmodel '_startTime.1D']);
            fMat    = fullfile(fileparts(replace(fOut,'.nii.gz','')),['cond-visOn_model-' HRmodel '_stats.xmat.1D']);
            fMatFig = fullfile(fileparts(replace(fOut,'.nii.gz','')),['cond-visOn_model-' HRmodel '_stats.xmat.fig']);
        end
        
        fSes(1,E).fIn     = fIn;
        fSes(1,E).fFit    = fFit;
        fSes(1,E).fResid  = fResid;
        fSes(1,E).fStat   = fStat;
        fSes(1,E).fMat    = fMat;
        fSes(1,E).fMatFig = fMatFig;
        switch HRmodel
            case {'TENT' 'TENTzero'}
                fSes(1,E).fResp   = fResp;
            case {'SPMG2' 'SPMG3'}
            otherwise
                dbstack; error('figure that out')
        end
        


        %% Contruct afni command
        cmdTmp = {};
        if size(fList,2)>1
            cmdTmp{end+1} = ['echo ''  ''runCat/' num2str(size(fList,1)) ' -- echo' num2str(E) '/' num2str(size(fList,2))];
        else
            cmdTmp{end+1} = ['echo ''  ''runCat/' num2str(size(fList,1))];
        end
        
        switch HRmodel
            case {'TENT' 'TENTzero'}
                [cmdTmpTmp,param.funDsgn.nReg] = afniCmd(fIn,fMask,fStim,param.nDummy,param.tr,startSeq,durSeq,condSeq,HRmodel,param.funDsgn.label,[],fResp,fFit,fResid,fMat,fStat,verbose,param.nDummyRemoved,param.trDecon,param.dryRun);
            case {'SPMG2' 'SPMG3'}
                [cmdTmpTmp,param.funDsgn.nReg] = afniCmd(fIn,fMask,fStim,param.nDummy,param.tr,startSeq,durSeq,condSeq,HRmodel,param.funDsgn.label,[],[]   ,fFit,fResid,fMat,fStat,verbose,param.nDummyRemoved,[]           ,param.dryRun);
            otherwise
                dbstack; error('figure that out')
        end
        if param.dryRun
            system(strjoin([{srcAfni} cmdTmpTmp],newline))
        end
        
        if ~exist(fStat,'file') || force
            cmdTmp = [cmdTmp cmdTmpTmp];
            if exist('fResp','var') && ~isempty(fResp); cmdTmp{end+1} = ['echo ''   ''' fResp]; end
            cmdTmp{end+1} = ['echo ''   ''' fStat];
            cmdTmp{end+1} = ['echo ''   ''' fMat];
        else
            if exist('fResp','var') && ~isempty(fResp); cmdTmp{end+1} = ['echo ''   ''' fResp]; end
            cmdTmp{end+1} = ['echo ''   ''' fStat];
            cmdTmp{end+1} = ['echo ''   ''' fMat];
            cmdTmp{end+1} = 'echo ''   ''already done, skipping';
        end

        fSes(1,E).cmd =  strjoin(cmdTmp,newline);

        cmd = [cmd cmdTmp];
    end

    %% Run afni command
    [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end

end


function [cmd,nReg] = afniCmd(fIn,fStim,fMask,param,fResp,fRespStd,fFit,fResid,fMat,fStat,verbose,dryRun)
% function [cmd,nReg] = afniCmd2(fIn,fMask,fStim,nDummyIgnore,tr,startSeq,durSeq,condSeq,HRmodel,label,param,fResp,fFit,fResid,fMat,fStat,verbose,nDummyRemoved,trDecon,dryRun,nFrame)
% param.nDummyRemoved [int]: number of initial frames that are already removed from the
% timeseries. The stimulus timeseries must therefore be adjusted
% accordingly.
% param.nDummyIgnore [int]: number of initial frames to ignore from the input
% timeseries. Only frames after these will be feed to 3dDeconvolve using
% the [nDummyIgnore..$] notation.
% nDummy [int]: total number of dummy initial frames
tr      = param.tr;
trDecon = param.trDecon;
nDummy = param.nDummyIgnore + param.nDummyRemoved;
fIn = cellstr(fIn);
cmd = {'3dDeconvolve -overwrite \'};
% cmdTmp{end+1} = ['-force_TR ' num2str(trStim) ' \'];
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
cmd{end+1} = '-polort A \';
cmd{end+1} = ['-stim_times_subtract ' num2str(tr*nDummy,'%f') ' \'];

% New way (only implemented for dry runs so far)
dsgn = param.dsgn;
nRegAll = [];
cmd{end+1} = ['-num_stimts ' num2str(dsgn.condK) ' \'];
for k = 1:dsgn.condK
    cmd{end+1} = ['-stim_label ' num2str(k) ' ' [char(dsgn.task) '_' dsgn.condLabel{k}] ' \'];

    % write design to file
    if dryRun
        fStim = [tempname '_startTime.1D' ];
    end
    fido = fopen(fStim{k}, 'w');
    if ~iscell(fIn); dbstack; error('fIn must be type cell'); end
    for i = 1:length(fIn)
        fprintf(fido,'%.3f ',dsgn.onsetList((k-1)==dsgn.cond));
        fprintf(fido,'\n');
    end
    fclose(fido);


    % k = 1;

    switch param.model
        case 'SPMG2'
            dbstack; error('double-check that')
            nReg = 2;
            if max(abs(diff(durSeq)))/max(durSeq) > 0.0001; dbstack; error('stim duration cannot be different across trials'); end
            cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' ''' HRmodel '(' num2str(mean(durSeq),'%0.3f') ')'' \'];
        case 'SPMG3'
            dbstack; error('double-check that')
            nReg = 3;
            if max(abs(diff(durSeq)))/max(durSeq) > 0.0001; dbstack; error('stim duration cannot be different across trials'); end
            cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' ''' HRmodel '(' num2str(mean(durSeq),'%0.3f') ')'' \'];
        case 'TENT'
            dbstack; error('code that')
        case 'TENTzero'
            eTime     = dsgn.onsetList(dsgn.cond==(k-1));
            eTimeNext = find(dsgn.cond==(k-1))+1;
            if eTimeNext(end) > length(dsgn.onsetList)
                eTimeNext(end) = [];
                eTimeNext = dsgn.onsetList(eTimeNext);
                eTimeNext(end+1) = (param.nFrame + param.nDummyRemoved) * tr;
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
cmd{end+1} = ['-TR_times ' num2str(trDecon,'%f') ' \'];

if dryRun
    dbstack; error('code that')
    fMat = replace(fMat,'_stats.xmat.1D','_stats.xmatPerTrial.1D');
    fMat  = [tempname '_stats.xmat.1D'];    
end




if ~dryRun
    cmd{end+1} = ['-fitts ' fFit ' \'];
    cmd{end+1} = ['-errts ' fResid ' \'];
    cmd{end+1} = '-bout -fout -tout \';
end
% cmdTmp{end+1} = ['-TR_times ' num2str(trStim) ' \'];
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
