function [fRun,fSes,fSes_echoCat,param] = getAct3(volTs,dsgn,fMask,param,force,verbose)
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
        fRun(R,1) = plotDsgnMat(fRun(R,1),param,fVolTs(R,1),volTs(R,1),verboseThis);
    end
end
if ~isempty(fSes)
    fSes = plotDsgnMat(fSes,param,fVolTs,volTs,verboseThis);
end

%% Refactor
if ~isempty(fRun)
    for R = 1:size(fRun,1)
        if isfield(fRun,'fResp')
            tmp(R,1).afni = rmfield(fRun(R,1),'fResp');
            tmp(R,1).fs.fRespTs = fRun(R,1).fResp;
        else
            tmp(R,1).afni = fRun(R,1);
        end
        tmp(R,1).fs.fMask   = fMask;
    end
    fRun = tmp; clear tmp
end
if ~isempty(fSes)
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


%% Simplify stat outputs
disp('Simplifying stat outputs')
forceThis = force;
if ~isempty(f)
    cmd = {srcAfni};
    for i = 1:numel(f)
        %%% Fstat
        fIn = f(i).afni.fStat;
        fOut = replace(fIn,'_stats.nii.gz','_fullFval.nii.gz');
        f(i).fs.fFullF = fOut;
        if forceThis || ~exist(fOut,'file')
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fIn '[Full_Fstat]'];
        end

        %%% Pval
        fIn  = f(i).fs.fFullF;
        fOut = replace(fIn,'_fullFval.nii.gz','_fullPval.nii.gz');
        f(i).fs.fFullP = fOut;
        if force || ~exist(fOut,'file')
            cmd{end+1} = ['df=$(3dAttribute BRICK_STATAUX ' fIn ')'];
            cmd{end+1} = 'df1=$(echo $df | awk ''{print $(NF-1)}'')';
            cmd{end+1} = 'df2=$(echo $df | awk ''{print $NF}'')';
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = ['-a ' fIn ' \'];
            cmd{end+1} = '-expr "1-stat2cdf(a,4,$df1,$df2,0)" 2> /dev/null';
        end

        %%% Qval (fdr)
        fIn = f(i).afni.fStat;
        fOut = replace(fIn,'_stats.nii.gz','_fullQval.nii.gz');
        f(i).fs.fFullQ = fOut;
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dFDR -overwrite -qval \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = ['-input ' fIn ' \'];
            cmd{end+1} = ['-mask '  fMask];
            cmd{end+1} = '3dbucket -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = [fOut '[FDRq:Full_Fstat]'];
        end

        %%% Coef
        switch param.model
            case {'SPMG2'}
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
    %%% Baseline (simple temporal average)
    fOut      = cell(size(f));
    fOutExist = false(size(f));
    for i = 1:numel(f)
        % average within run
        fOutR = cell(size(f(i).afni.fIn));
        for R = 1:size(f(i).afni.fIn,1)
            fIn  = f(i).afni.fIn{R};
            fOutR{R} = replace(fIn,'preproc_volTs.nii.gz','av_preproc_volTs.nii.gz');
            if force || ~exist(fOutR{R},'file')
                cmd{end+1} = '3dTstat -overwrite \';
                cmd{end+1} = ['-prefix ' fOutR{R} ' \'];
                cmd{end+1} = fIn;
            end
        end

        if size(f(i).afni.fIn,1)==1
            f(i).fs.fBaseAv = char(fOutR);
        else
            % catenate runs
            fIn = fOutR;
            fOut{i} = fIn{1}; fOut{i} = strsplit(fOut{i},'_run-'); fOut{i}{2} = strsplit(fOut{i}{2},'_'); fOut{i}{2}{1} = 'cat'; fOut{i}{2} = strjoin(fOut{i}{2},'_'); fOut{i} = strjoin(fOut{i},'_run-');
            if force || ~exist(fOut{i},'file')
                cmd{end+1} = '3dTcat -overwrite \';
                cmd{end+1} = ['-prefix ' fOut{i} ' \'];
                cmd{end+1} = strjoin(fIn,' ');
            end
            % average across runs
            fIn = fOut{i};
            fOut{i} = replace(fOut{i},'_run-cat','_run-avCat');
            if force || ~exist(fOut{i},'file')
                cmd{end+1} = '3dTstat -overwrite \';
                cmd{end+1} = ['-prefix ' fOut{i} ' \'];
                cmd{end+1} = fIn;
            end
            f(i).fs.fBaseAv = fOut{i};
        end
    end
    % if any(~fOutExist)
    %     nFrame = [];
    %     tsAv   = [];
    %     for R = 1:size(f,1)
    %         fIn = f(R,1).afni.fIn;
    %         if length(fIn)>1 && R==size(f,1)
    %             mri.vol = tsAv./nFrame;
    %             MRIwrite(mri,fOut{R,1});
    %         elseif length(fIn)==1
    %             mri = MRIload3(char(fIn),[],[],0);
    %             mri.vol = sum(mri.vol,4);
    %             if R==1
    %                 tsAv   = mri.vol;
    %                 nFrame = mri.nframes;
    %             else
    %                 tsAv   = tsAv + mri.vol;
    %                 nFrame = nFrame + mri.nframes;
    %             end
    %             mri.vol = mri.vol./mri.nframes;
    %             MRIwrite(mri,fOut{R,1});
    %         else
    %             dbstack; error('X');
    %         end
    %     end
    % end
end
% if ~isempty(fRun) && ~isempty(fSes)
%     fRun = f(1:end-1,:);
%     fSes = f(end,:);
%     f    = [];
% else
%     dbstack; error('X');
% end





% if nEcho>1
%     dbstack; error('code that');
%     switch param.model
%         case {'SPMG2' 'SPMG3'}
%             dbstack; error('code that');
%         case {'TENT' 'TENTzero'}
%             for i = 1:numel(fSes_echoRms)
%                 fIn = fSes_echoRms(i).fStat;
%                 fOut = replace(fIn,'_stats.nii.gz','_respFval.nii.gz');
%                 fSes_echoRms(i).fRespStat = fOut;
%                 if forceThis || ~exist(fOut,'file')
%                     cmd{end+1} = '3dbucket -overwrite \';
%                     cmd{end+1} = ['-prefix ' fOut ' \'];
%                     cmd{end+1} = [fIn '[Full_Fstat]'];
%                 end
%             end
%             for i = 1:numel(fRun_echoRms)
%                 fIn = fRun_echoRms(i).fStat;
%                 fOut = replace(fIn,'_stats.nii.gz','_respFval.nii.gz');
%                 fRun_echoRms(i).fRespStat = fOut;
%                 if forceThis || ~exist(fOut,'file')
%                     cmd{end+1} = '3dbucket -overwrite \';
%                     cmd{end+1} = ['-prefix ' fOut ' \'];
%                     cmd{end+1} = [fIn '[Full_Fstat]'];
%                 end
%             end
%         otherwise
%             dbstack; error('code that');
%     end
% end


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



%%% Convert SPMG2 cartesian responses coefficient (gamma + first derivative) to polar (amplitude + delay) coefficient
switch param.model
    case {'SPMG2'}
        % dbstack; error('double-check that')
        disp('convert hrf+derivative cartesian coefficients to polar coefficients')

        for i = 1:size(f,1)
            fIn     = f(i).fs.fCoef;
            fOut    = replace(fIn,'_coef.nii.gz','_coefPol.nii.gz');
            fOutFig = replace(fIn,'_coef.nii.gz','_coefPol.fig');
            fFDR    = f(i).fs.fFullQ;
            f(i).fs.fCoefPol = fOut;

            verboseThis = verbose;
            forceThis   = force;
            if forceThis || ~exist(fOut,'file') || ~exist(fOutFig,'file')
                % hMat = figure('WindowStyle','docked');
                hMat = figure('Visible','off');

                coef = MRIread(fIn);
                fdr  = MRIread(fFDR);
                mask = MRIread(fMask);
                mask = mask.vol & fdr.vol<0.05;

                coef.vol = complex(coef.vol(:,:,:,1),coef.vol(:,:,:,2));
                scatter(real(coef.vol(mask)),imag(coef.vol(mask)));
                ax = gca; ax.DataAspectRatio = [1 1 1];
                grid on
                axis([-1 1 -1 1].*max(abs(axis)))
                xline(0,'k'); yline(0,'k');
                xlabel('SPM canon (coef)')
                ylabel('SPM canon derivative (coef)')

                %get principal vector
                slp = real(coef.vol(mask))\imag(coef.vol(mask));
                hRef = refline(slp,0); hRef.Color = 'r';
                v = complex(1,slp); v = v./abs(v);
                title([num2str(angle(v)/pi*180,'%0.1f°') ' deviation from expected HR delay'])

                % subtract that vector orientation from data
                coefPol = coef;
                coefPol.vol(:,:,:,1) = abs(coef.vol);
                coefPol.vol(:,:,:,2) = wrapToPi( angle(coef.vol) - angle(v) );
                MRIwrite(coefPol,fOut);

                if verboseThis>1
                    hMat.Visible = 'on';
                    hMat.WindowStyle = 'docked';
                    savefig(hMat,fOutFig,'compact')
                else
                    set(hMat, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
                    savefig(hMat,fOutFig,'compact')
                    close(hMat)
                end

                disp(' done')
            else
                disp(' already done, skipping')
            end
        end



        % fIn     = fSes(i).fs.fCoef;
        % fOut    = replace(fIn,'_coef.nii.gz','_coefPol.nii.gz');
        % fOutFig = replace(fIn,'_coef.nii.gz','_coefPol.fig');
        % fSes(i).fs.fCoefPol = fOut;
        % 
        % 
        % verboseThis = verbose;
        % forceThis   = 1;
        % if forceThis || ~exist(fOut,'file') || ~exist(fOutFig,'file')
        %     % hMat = figure('WindowStyle','docked');
        %     hMat = figure('Visible','off');
        % 
        % 
        %     coef = MRIread(fIn);
        %     fdr  = MRIread(fSes(i).fs.fFullQ);
        %     mask = MRIread(fMask);
        %     mask = mask.vol & fdr.vol<0.05;
        % 
        %     coef.vol = complex(coef.vol(:,:,:,1),coef.vol(:,:,:,2));
        %     scatter(real(coef.vol(mask)),imag(coef.vol(mask)));
        %     ax = gca; ax.DataAspectRatio = [1 1 1];
        %     grid on
        %     axis([-1 1 -1 1].*max(abs(axis)))
        %     xline(0,'k'); yline(0,'k');
        %     xlabel('SPM canon (coef)')
        %     ylabel('SPM canon derivative (coef)')
        % 
        %     %get principal vector
        %     slp = real(coef.vol(mask))\imag(coef.vol(mask));
        %     hRef = refline(slp,0); hRef.Color = 'r';
        %     v = complex(1,slp); v = v./abs(v);
        %     title([num2str(angle(v)/pi*180,'%0.1f°') ' deviation from expected HR delay'])
        % 
        %     % subtract that vector orientation from data
        %     coefPol = coef;
        %     coefPol.vol(:,:,:,1) = abs(coef.vol);
        %     coefPol.vol(:,:,:,2) = wrapToPi( angle(coef.vol) - angle(v) );
        %     MRIwrite(coefPol,fOut);
        % 
        %     if verboseThis>1
        %         hMat.Visible = 'on';
        %         hMat.WindowStyle = 'docked';
        %     else
        %         close(hMat)
        %     end
        % 
        %     disp(' done')
        % else
        %     disp(' already done, skipping')
        % end
    case {'SPMG3'}
        dbstack; error('code that')
    case {'TENT' 'TENTzero'}
    otherwise
        dbstack; error('code that');
end






%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manipulate baselines %
%%%%%%%%%%%%%%%%%%%%%%%%%

forceThis = force;
cmd = {srcAfni};
if ~isempty(f)
    for i = 1:numel(f)
        % Extract fitted baseline
        fIn  = f(i).afni.fStat;
        fOut = replace(f(i).fs.fFullF,'_fullFval.nii.gz','_fitBase.nii.gz');
        f(i).fs.fBaseFit = fOut;
        if force || ~exist(fOut,'file')
            buck = num2str(1:size(f(i).afni.fIn,1),'Run#%iPol#0_Coef,'); buck(end) = [];
            cmd{end+1} = '3dcalc -overwrite -TR 999 \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = ['-a ' fIn '[' buck '] \'];
            cmd{end+1} = '-expr a';
        end
        % Average fitted baseline across runs
        fIn = fOut;
        fOut = replace(fOut,'run-cat','run-avCat');
        if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
        f(i).fs.fBaseFitAv = fOut;
        cmd{end+1} = '3dTstat -overwrite \';
        cmd{end+1} = ['-prefix ' fOut ' \'];
        cmd{end+1} = fIn;
        
        if ~ismember(param.model,{'TENT' 'TENTzero'}); continue; end
        % Add this fitted baseline to response fits
        fIn   = f(i).fs.fRespTs;
        fOut = replace(fIn,'_resp.nii.gz','_respOnBaseFit.nii.gz');
        fBase   = f(i).fs.fBaseFitAv;
        f(i).fs.fRespTsOnBaseFit = fOut;
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = ['-a ' fIn   ' \'];
            cmd{end+1} = ['-b ' fBase ' \'];
            cmd{end+1} = '-expr ''a+b''';
        end
        % Add ts average baseline to response fits
        fIn   = f(i).fs.fRespTs;
        fOut  = replace(fIn,'_resp.nii.gz','_respOnBaseAv.nii.gz');
        fBase = f(i).fs.fBaseAv;
        f(i).fs.fRespTsOnBaseAv = fOut;
        if force || ~exist(fOut,'file')
            cmd{end+1} = '3dcalc -overwrite \';
            cmd{end+1} = ['-prefix ' fOut ' \'];
            cmd{end+1} = ['-a ' fIn   ' \'];
            cmd{end+1} = ['-b ' fBase ' \'];
            cmd{end+1} = '-expr ''a+b''';
        end
    end
end



if param.skipRun && ~param.skipCat
    fRun = [];
    fSes = f;
elseif ~param.skipRun && ~param.skipCat
    fRun = f(1:end-1);
    fSes = f(end);
    clear f
else
    dbstack; error('fix that mess')
    if ~isempty(fRun) && ~isempty(fSes)
        fRun = f(1:end-1,:);
        fSes = f(end,:);
        f    = [];
    else
        dbstack; error('X');
    end
end

% 
% % Fitted baseline
% forceThis = force;
% cmd = {srcAfni};
% %%% Run by run
% if ~isempty(fRun)
%     for i = 1:numel(fRun)
%         fIn  = fRun(i).afni.fStat;
%         fOut = replace(fRun(i).fs.fFullF,'_fullFval.nii.gz','_fitBase.nii.gz');
%         fRun(i).fs.fBaseFit = fOut;
%         if force || ~exist(fOut,'file')
%             buck = num2str(1,'Run#%iPol#0_Coef,'); buck(end) = [];
% 
%             cmd{end+1} = '3dcalc -overwrite \';
%             cmd{end+1} = ['-prefix ' fOut ' \'];
%             cmd{end+1} = ['-a ' fIn '[' buck '] \'];
%             cmd{end+1} = '-expr a';
%         end
%     end
% end
% %%% Whole session
% if ~isempty(fSes)
%     for i = 1:numel(fSes)
%         % Fitted baseline
%         fIn  = fSes(i).afni.fStat;
%         fOut = replace(fSes(i).fs.fFullF,'_fullFval.nii.gz','_fitBase.nii.gz');
%         fSes(i).fs.fBaseFit = fOut;
%         if force || ~exist(fOut,'file')
%             buck = num2str(1:nRun,'Run#%iPol#0_Coef,'); buck(end) = [];
% 
%             cmd{end+1} = '3dcalc -overwrite \';
%             cmd{end+1} = ['-prefix ' fOut ' \'];
%             cmd{end+1} = ['-a ' fIn '[' buck '] \'];
%             cmd{end+1} = '-expr a';
% 
%             % %add the mean cross-run mean as the first frame
%             % fTmp = replace(fOut,'.nii.gz','TMP.nii.gz');
%             % cmd{end+1} = '3dTstat -overwrite \';
%             % cmd{end+1} = '-mean \';
%             % cmd{end+1} = ['-prefix ' fTmp ' \'];
%             % cmd{end+1} = [fIn '[' buck ']'];
%             % cmd{end+1} = '3dTcat -overwrite \';
%             % cmd{end+1} = ['-prefix ' fOut ' \'];
%             % cmd{end+1} = [fTmp ' ' fOut];
%         end
%     end
% 
%     % Add baseline to response
%     switch param.model
%         case {'SPMG2' 'SPMG3'}
%         case {'TENT' 'TENTzero'}
%             fIn   = fSes(i).fs.fRespTs;
% 
%             fOut  = replace(fIn,'_resp.nii.gz','_respOnBaseAv.nii.gz');
%             fBase = fSes(i).fs.fBaseAv;
%             fSes(i).fs.fRespTsOnBaseAv = fOut;
%             if force || ~exist(fOut,'file')
%                 cmd{end+1} = '3dcalc -overwrite \';
%                 cmd{end+1} = ['-prefix ' fOut ' \'];
%                 cmd{end+1} = ['-a ' fIn   ' \'];
%                 cmd{end+1} = ['-b ' fBase '[0] \'];
%                 cmd{end+1} = '-expr ''a+b''';
%             end
%             fOut = replace(fIn,'_resp.nii.gz','_respOnBaseFit.nii.gz');
%             fBase   = fSes(i).fs.fBaseFit;
%             fSes(i).fs.fRespTsOnBaseFit = fOut;
%             if force || ~exist(fOut,'file')
%                 cmd{end+1} = '3dcalc -overwrite \';
%                 cmd{end+1} = ['-prefix ' fOut ' \'];
%                 cmd{end+1} = ['-a ' fIn   ' \'];
%                 cmd{end+1} = ['-b ' fBase '[0] \'];
%                 cmd{end+1} = '-expr ''a+b''';
%             end
%         otherwise
%             dbstack; error('X');
%     end
% 
% end
% 
% 
% 
% 
% %%% Each run
% if ~isempty(fRun)
%     dbstack; error('code that')
%     for R = 1:nRun
%         for E = 1:nEcho+1
%             %%%% Define file names
%             if E>nEcho
%                 if nEcho==1; continue; end
%                 fResp = fRun_echoRms(R,1).fResp;
%                 fBase = replace(fResp,'_resp.nii.gz','_base.nii.gz');
%                 fRespOnBase = replace(fResp,'_resp.nii.gz','_respOnBase.nii.gz');
%                 fRun_echoRms(R,1).fBase = fBase;
%                 fRun_echoRms(R,1).fRespOnBase = fRespOnBase;
%                 % if useFittedBaseline
%                 % Use fitted baseline
%                 fStat = fRun_echoRms(R,1).fStat;
%                 % else
%                 %     % Use plain average as baseline
%                 %     copyfile(fFunc.fAvCatAvEchoRms,fBase);
%                 % end
%             else
%                 fResp = fRun(R,E).fResp;
%                 fBase = replace(fResp,'_resp.nii.gz','_base.nii.gz');
%                 fRespOnBase = replace(fResp,'_resp.nii.gz','_respOnBase.nii.gz');
%                 fRun(R,E).fBase = fBase;
%                 fRun(R,E).fRespOnBase = fRespOnBase;
%                 % if useFittedBaseline
%                 % Use fitted baseline
%                 fStat = fRun(R,E).fStat;
%                 % else
%                 %     % Use plain average as baseline
%                 %     copyfile(fFunc.fAvCatAv{E},fBase);
%                 % end
%             end
% 
%             % if useFittedBaseline
%             %%% Average fitted baselines across runs
%             cmd{end+1} = '3dTstat -overwrite \';
%             cmd{end+1} = '-mean \';
%             cmd{end+1} = ['-prefix ' fBase ' \'];
%             % buck = num2str(1:size(fRun,1),'Run#%iPol#0_Coef,'); buck(end) = [];
%             buck = 'Run#1Pol#0_Coef';
%             cmd{end+1} = [fStat '[' buck ']'];
%             % end
% 
%             %%% Add baseline to response
%             cmd{end+1} = '3dcalc -overwrite \';
%             cmd{end+1} = ['-prefix ' fRespOnBase ' \'];
%             cmd{end+1} = ['-a ' fBase ' \'];
%             cmd{end+1} = ['-b ' fResp ' \'];
%             cmd{end+1} = '-expr ''a+b''';
%         end
%     end
% end
% 
% %%% Whole session
% if ~isempty(fSes)
%     for i = 1:numel(fSes)
%         % Simple temporal average as a baseline
%         fOut = replace(fSes(i).fs.fFullF,'_fullFval.nii.gz','_tsAv.nii.gz');
%         fSes(i).fs.fBaseAv = fOut;
%         if force || ~exist(fOut,'file')
%             tsAv = volTs(1);
%             tsAv.vol = zeros([tsAv.volsize size(volTs,1)]);
%             nframes  = zeros([1 1 1        size(volTs)  ]);
%             for r = 1:size(volTs,1)
%                 if ~isempty(volTs(r).vol)
%                     tmp = volTs(r);
%                 else
%                     tmp = MRIload3(volTs(r),[],[],0);
%                 end
%                 tsAv.vol(:,:,:,r) = mean(tmp.vol,4);
%                 nframes( 1,1,1,r) = tmp.nframes;
%             end
%             % add the mean cross-run mean as the first frame
%             tsAv.vol = cat(4,sum(tsAv.vol.*nframes,4)./sum(nframes),tsAv.vol);
%             MRIwrite(tsAv,fOut);
%         end
% 
%         % Fitted baseline
%         fIn  = fSes(i).afni.fStat;
%         fOut = replace(fSes(i).fs.fFullF,'_fullFval.nii.gz','_fitBase.nii.gz');
%         fSes(i).fs.fBaseFit = fOut;
%         if force || ~exist(fOut,'file')
%             buck = num2str(1:nRun,'Run#%iPol#0_Coef,'); buck(end) = [];
% 
%             cmd{end+1} = '3dcalc -overwrite \';
%             cmd{end+1} = ['-prefix ' fOut ' \'];
%             cmd{end+1} = ['-a ' fIn '[' buck '] \'];
%             cmd{end+1} = '-expr a';
% 
%             %add the mean cross-run mean as the first frame
%             fTmp = replace(fOut,'.nii.gz','TMP.nii.gz');
%             cmd{end+1} = '3dTstat -overwrite \';
%             cmd{end+1} = '-mean \';
%             cmd{end+1} = ['-prefix ' fTmp ' \'];
%             cmd{end+1} = [fIn '[' buck ']'];
%             cmd{end+1} = '3dTcat -overwrite \';
%             cmd{end+1} = ['-prefix ' fOut ' \'];
%             cmd{end+1} = [fTmp ' ' fOut];
%         end
%     end
% 
%     % Add baseline to response
%     switch param.model
%         case {'SPMG2' 'SPMG3'}
%         case {'TENT' 'TENTzero'}
%             fIn   = fSes(i).fs.fRespTs;
% 
%             fOut  = replace(fIn,'_resp.nii.gz','_respOnBaseAv.nii.gz');
%             fBase = fSes(i).fs.fBaseAv;
%             fSes(i).fs.fRespTsOnBaseAv = fOut;
%             if force || ~exist(fOut,'file')
%                 cmd{end+1} = '3dcalc -overwrite \';
%                 cmd{end+1} = ['-prefix ' fOut ' \'];
%                 cmd{end+1} = ['-a ' fIn   ' \'];
%                 cmd{end+1} = ['-b ' fBase '[0] \'];
%                 cmd{end+1} = '-expr ''a+b''';
%             end
%             fOut = replace(fIn,'_resp.nii.gz','_respOnBaseFit.nii.gz');
%             fBase   = fSes(i).fs.fBaseFit;
%             fSes(i).fs.fRespTsOnBaseFit = fOut;
%             if force || ~exist(fOut,'file')
%                 cmd{end+1} = '3dcalc -overwrite \';
%                 cmd{end+1} = ['-prefix ' fOut ' \'];
%                 cmd{end+1} = ['-a ' fIn   ' \'];
%                 cmd{end+1} = ['-b ' fBase '[0] \'];
%                 cmd{end+1} = '-expr ''a+b''';
%             end
%         otherwise
%             dbstack; error('X');
%     end
% 
% end

%%% Run command
disp('Writing baselines')
if length(cmd)>1
    if verbose
        [status,cmdout] = system(strjoin(cmd,newline),'-echo'); if status || isempty(cmdout) || contains(cmdout,'error','IgnoreCase',true); dbstack; error(cmdout); error('x'); end
    else
        [status,cmdout] = system(strjoin(cmd,newline)); if status || isempty(cmdout) || contains(cmdout,'error','IgnoreCase',true); dbstack; error(cmdout); error('x'); end
    end
    disp(' done')
else
    disp(' already done, skipping')
end


%% %%%%%%%%%%%%%%%%%%%%%%


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







if 0

    %% Reconcile output with convention
    forceThis = force;
    disp('Catenating functional analyses')
    cmd = {srcAfni};
    fieldListIn =  {'fBase' 'fResp' 'fRespOnBase' 'fRespOnBaseMovie' 'fRespOnBaseMovieHighBit' 'fRespStat' 'fFit' 'fResid' 'fStat' 'fMat'    'fF' 'fCoef'};
    fieldListOut = {'base'  'resp'  'respOnBase'  'respOnBaseMovie'  'respOnBaseMovieHighBit'  'respF'     'fit'  'resid'  'stat'  'dsgnMat' 'F'  'coef' };
    for i = 1:length(fieldListIn)
        fieldIn = fieldListIn{i};
        fieldOut = fieldListOut{i};
        % if strcmp(fieldOut,'respF'); keyboard; end
        if isfield(fRun,fieldIn)
            files.(fieldOut).f = reshape({fRun.(fieldIn)},size(fRun));
        end
        if nEcho>1 && isfield(fRun_echoRms,fieldIn)
            files.(fieldOut).fEchoRms = reshape({fRun_echoRms.(fieldIn)},size(fRun_echoRms));
        end
        if nRun>1 && ~isempty(fSes) && isfield(fSes,fieldIn) && ~param.skipCat
            files.(fieldOut).fSes = reshape({fSes.(fieldIn)},size(fSes));
        end
        if nRun>1 && nEcho>1 && isfield(fSes_echoRms,fieldIn) && ~isempty(fSes) && ~param.skipCat
            files.(fieldOut).fSesEchoRms = reshape({fSes_echoRms.(fieldIn)},size(fSes_echoRms));
        end

        %%% catenate echo
        if nEcho>1
            if any(ismember(fieldOut,{'base' 'respF'}))
                if isfield(fRun,fieldIn)
                    fIn = reshape({fRun.(fieldIn)},size(fRun));
                    files.(fieldOut).fEchoCat = cell([size(fIn,1) 1]);
                    for R = 1:size(fIn,1)
                        fOut = fIn{R,1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-cat'; fOut = {strjoin(fOut,'_')};
                        files.(fieldOut).fEchoCat(R,1) = fOut; if ~exist(fileparts(char(fOut)),'dir'); mkdir(fileparts(char(fOut))); end
                        if forceThis || ~exist(char(fOut),'file')
                            if any(ismember(fieldOut,{'respF'}))
                                cmd{end+1} = '3dbucket -overwrite \';
                            else
                                cmd{end+1} = '3dTcat -overwrite \';
                            end
                            cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                            cmd{end+1} = strjoin(fIn(R,:));
                        end
                    end
                end
                if ~isempty(fSes) && isfield(fSes,fieldIn) && ~param.skipCat
                    fIn = {fSes.(fieldIn)};
                    fOut = fIn{1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'echo-')} = 'echo-cat'; fOut = {strjoin(fOut,'_')};
                    files.(fieldOut).fSesEchoCat = fOut; if ~exist(fileparts(char(fOut)),'dir'); mkdir(fileparts(char(fOut))); end
                    if forceThis || ~exist(char(fOut),'file')
                        if any(ismember(fieldOut,{'respF'}))
                            cmd{end+1} = '3dbucket -overwrite \';
                        else
                            cmd{end+1} = '3dTcat -overwrite \';
                        end
                        cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                        cmd{end+1} = strjoin(fIn);
                    end
                end
            end
        end

        %%% catenate runs
        if nRun>1
            if any(ismember(fieldOut,{'base' 'respF'}))
                if isfield(fRun,fieldIn)
                    fIn = reshape({fRun.(fieldIn)},size(fRun));
                    files.(fieldOut).fCat = cell([1 size(fIn,2)]);
                    for E = 1:size(fIn,2)
                        fOut = fIn{1,E}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'run-')} = 'run-cat'; fOut = {strjoin(fOut,'_')};
                        fOut = strsplit(char(fOut),filesep); fOut{end} = ['run-cat_' fOut{end}]; fOut = {strjoin(fOut,filesep)};
                        files.(fieldOut).fCat(1,E) = fOut; if ~exist(fileparts(char(fOut)),'dir'); mkdir(fileparts(char(fOut))); end
                        if forceThis || ~exist(char(fOut),'file')
                            if any(ismember(fieldOut,{'respF'}))
                                cmd{end+1} = '3dbucket -overwrite \';
                            else
                                cmd{end+1} = '3dTcat -overwrite \';
                            end
                            cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                            cmd{end+1} = strjoin(fIn(:,E));
                        end
                    end
                    if nEcho>1
                        fIn = reshape({fRun_echoRms.(fieldIn)},size(fRun_echoRms));
                        fOut = fIn{1}; fOut = strsplit(fOut,'_'); fOut{contains(fOut,'run-')} = 'run-cat'; fOut = {strjoin(fOut,'_')};
                        fOut = strsplit(char(fOut),filesep); fOut{end} = ['run-cat_' fOut{end}]; fOut = {strjoin(fOut,filesep)};
                        files.(fieldOut).fCatEchoRms = fOut; if ~exist(fileparts(char(fOut)),'dir'); mkdir(fileparts(char(fOut))); end
                        if forceThis || ~exist(char(fOut),'file')
                            if any(ismember(fieldOut,{'respF'}))
                                cmd{end+1} = '3dbucket -overwrite \';
                            else
                                cmd{end+1} = '3dTcat -overwrite \';
                            end
                            cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                            cmd{end+1} = strjoin(fIn(:));
                        end
                    end
                end
            end
        end
    end


    %% Produce fdr maps
    %clean up the old catenation approach
    switch param.model
        case {'TENT' 'TENTzero'}
            if param.skipCat
                if isfield(files.respF,'fCat')
                    if exist(char(files.respF.fCat),'file')
                        delete(char(files.respF.fCat))
                    end
                    files.respF = rmfield(files.respF,'fCat');
                end
            end

            fieldList = fields(files.F);
            for i = 1:length(fieldList)
                f = files.respF.(fieldList{i});
                fOut = replace(f,'_respFval.nii.gz','_respFfdr.nii.gz');
                for E = 1:size(f,2)
                    for R = 1:size(f,1)
                        if force || ~exist(fOut{R,E},'file')
                            cmd{end+1} = ['df=$(3dAttribute BRICK_STATAUX ' f{R,E} ')'];
                            cmd{end+1} = 'df1=$(echo $df | awk ''{print $(NF-1)}'')';
                            cmd{end+1} = 'df2=$(echo $df | awk ''{print $NF}'')';
                            cmd{end+1} = '3dcalc -overwrite \';
                            cmd{end+1} = ['-prefix ' fOut{R,E} ' \'];
                            cmd{end+1} = ['-a ' f{R,E} ' \'];
                            cmd{end+1} = '-expr "1-stat2cdf(a,4,$df1,$df2,0)" 2> /dev/null';
                        end
                    end
                end
                files.respF_fdr.(fieldList{i}) = fOut;
            end

        case {'SPMG2' 'SPMG3'}
            f = files.F.f;
            for E = 1:size(f,2)
                for R = 1:size(f,1)
                    fOut = replace(f{R,E},'_fullF.nii.gz','_fullFfdr.nii.gz');
                    if force || ~exist(fOut,'file')
                        cmd{end+1} = ['df=$(3dAttribute BRICK_STATAUX ' char(f{R,E}) ')'];
                        cmd{end+1} = 'df1=$(echo $df | awk ''{print $(NF-1)}'')';
                        cmd{end+1} = 'df2=$(echo $df | awk ''{print $NF}'')';
                        cmd{end+1} = '3dcalc -overwrite \';
                        cmd{end+1} = ['-prefix ' char(fOut) ' \'];
                        cmd{end+1} = ['-a ' char(f{R,E}) ' \'];
                        cmd{end+1} = '-expr "1-stat2cdf(a,4,$df1,$df2,0)" 2> /dev/null';
                    end
                    files.F_fdr.f{R,E} = fOut;
                end
            end

        otherwise
            dbstack; error('figure that out')
    end



    % % if nEcho>1
    %     %%% for each echo
    %     allCandidate = {'f' 'fCat' 'fSes'};
    %     for i = 1:length(allCandidate)
    %         candidate = allCandidate{i};
    %         if isfield(files.respF,candidate)
    %             f = files.respF.(candidate);
    %             fOut = cell(size(f));
    %             for E = 1:nEcho
    %                 fOut{E} = replace(f{E},'_respFval.nii.gz','_respFfdr.nii.gz');
    %                 if force || ~exist(fOut{E},'file')
    %                     cmd{end+1} = ['df=$(3dAttribute BRICK_STATAUX ' f{E} ')'];
    %                     cmd{end+1} = 'df1=$(echo $df | awk ''{print $(NF-1)}'')';
    %                     cmd{end+1} = 'df2=$(echo $df | awk ''{print $NF}'')';
    %                     cmd{end+1} = '3dcalc -overwrite \';
    %                     cmd{end+1} = ['-prefix ' fOut{E} ' \'];
    %                     cmd{end+1} = ['-a ' f{E} ' \'];
    %                     cmd{end+1} = '-expr "1-stat2cdf(a,4,$df1,$df2,0)" 2> /dev/null';
    %                 end
    %             end
    %             files.respF_fdr.(candidate) = fOut;
    %         end
    %     end
    %
    %     %%% for echo rms
    %     candidate = {'fSesEchoRms' 'fEchoRms'};
    %     candidate = candidate(ismember(candidate,fields(files.respF)));
    %     if ~isempty(candidate)
    %         candidate = candidate{1};
    %         f = files.respF.(candidate);
    %         fOut = cell(size(f));
    %         E = 1;
    %         fOut{E} = replace(f{E},'_respFval.nii.gz','_respFfdr.nii.gz');
    %         if force || ~exist(fOut{E},'file')
    %             cmd{end+1} = ['df=$(3dAttribute BRICK_STATAUX ' f{E} ')'];
    %             cmd{end+1} = 'df1=$(echo $df | awk ''{print $(NF-1)}'')';
    %             cmd{end+1} = 'df2=$(echo $df | awk ''{print $NF}'')';
    %             cmd{end+1} = '3dcalc -overwrite \';
    %             cmd{end+1} = ['-prefix ' fOut{E} ' \'];
    %             cmd{end+1} = ['-a ' f{E} ' \'];
    %             cmd{end+1} = '-expr "1-stat2cdf(a,4,$df1,$df2,0)" 2> /dev/null';
    %         end
    %         files.respF_fdr.(candidate) = fOut;
    %     end
    % % else
    % %     % dbstack; error('code this');
    % % end

    files.cmd = strjoin(cmd,newline);
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


end








% %% Generate visualization commands
% if ~isempty(fAnat) && ~isempty(fAnat{1})
%     fT1w = fAnat{contains(fAnat,'_T1w.nii.gz')};
% 
%     ind = contains(fAnat,'echo-rms_satIndex.nii.gz');
%     if any(ind); fSatinIndex_echoRms = fAnat{ind}; else; fSatinIndex_echoRms = []; end
%     ind = contains(fAnat,'echo-rms_satIndexNumOnBack.nii.gz');
%     if any(ind); fSatinIndexPlus_echoRms = fAnat{ind}; else; fSatinIndexPlus_echoRms = []; end
%     ind = contains(fAnat,'echo-rms_satIndexNumAbs.nii.gz');
%     if any(ind); fSatinIndexAbs_echoRms = fAnat{ind}; else; fSatinIndexAbs_echoRms = []; end
%     ind = contains(fAnat,'echo-cat_satIndex.nii.gz');
%     if any(ind); fSatinIndex_echoCat = fAnat{ind}; else; fSatinIndex_echoCat = []; end
%     ind = contains(fAnat,'echo-cat_satIndexNumOnBack.nii.gz');
%     if any(ind); fSatinIndexPlus_echoCat = fAnat{ind}; else; fSatinIndexPlus_echoCat = []; end
%     ind = contains(fAnat,'echo-cat_satIndexNumAbs.nii.gz');
%     if any(ind); fSatinIndexAbs_echoCat = fAnat{ind}; else; fSatinIndexAbs_echoCat = []; end
% else
%     fT1w = [];
%     fSatinIndex_echoRms = [];
%     fSatinIndexPlus_echoRms = [];
%     fSatinIndexAbs_echoRms = [];
%     fSatinIndex_echoCat = [];
%     fSatinIndexPlus_echoCat = [];
%     fSatinIndexAbs_echoCat = [];
% end
% if ~isempty(fFmap) && ~isempty(fFmap{1})
%     ind = contains(fFmap,'rec-FA_TB1SRGE.nii.gz');
%     if any(ind); fB1 = fFmap{ind}; else; fB1 = []; end
%     % fB1 = fFmap{contains(fFmap,'rec-FA_TB1SRGE.nii.gz')};
% else
%     fB1 = [];
% end
% 
% %%% individual runs
% for I = 1:size(fRun,1)
%     E = 1; % multi-echo visualization not implemented
%     fRun(I,E).fUnder = fFunc.fAv{I,E};
% 
%     cmd = {srcAfni};
%     if isfield(param,'layout') && ~isempty(param.layout) && exist(param.layout,'file')
%         cmd{end+1} = ['afni -layout ' param.layout ' \'];
%     else
%         cmd{end+1} = 'afni \';
%     end
%     cmd{end+1} = ['-tbar run' num2str(I) 'echo ' num2str(E) ' \'];
%     cmd{end+1} = [fRun(I,E).fUnder ' \'];
%     cmd{end+1} = [fRun(I,E).fStat ' \'];
%     cmd{end+1} = [fRun(I,E).fResp ' &'];
%     cmd = strjoin(cmd,newline); % disp(cmd)
%     fRun(I,E).cmdVisAfni = cmd;
% 
%     if verbose; disp(cmd); end
%     if verbose>1; [status,cmdout] = system(cmd); end
% 
%     % https://afni.nimh.nih.gov/pub/dist/doc/program_help/README.driver.html
%     % https://afni.nimh.nih.gov/pub/dist/doc/program_help/plugout_drive.html
%     % afni -yesplugouts
%     % plugout_drive  -com 'SWITCH_SESSION A.afni'                       \
%     % -com 'OPEN_WINDOW A.axialimage geom=600x600+416+44 \
%     % ifrac=0.8 opacity=9'                         \
%     % -com 'OPEN_WINDOW A.sagittalimage geom=+45+430     \
%     % ifrac=0.8 opacity=9'                         \
%     % -com 'SWITCH_UNDERLAY anat'                        \
%     % -com 'SWITCH_OVERLAY strip'                        \
%     % -com 'SEE_OVERLAY +'                               \
%     % -com 'SET_DICOM_XYZ 7 12 2'                        \
%     % -com 'OPEN_WINDOW A.axialimage keypress=v'         \
%     % -quit
%     %
%     % SET_SUBBRICKS
% 
% 
% 
%     % fT1w = dir(fullfile(bidsDir,'anat','*proc-RMS_T1w.nii.gz')); fT1w = fullfile(fT1w.folder,fT1w.name);
%     % fB1 = dir(fullfile(bidsDir,'fmap','*rec-FA_TB1SRGE.nii.gz')); fB1 = fullfile(fB1.folder,fB1.name);
%     cmd = {srcFs};
%     cmd{end+1} = 'freeview \';
%     cmd{end+1} = [fRun(I).fUnder ' \'];
%     if ~isempty(fT1w)
%         cmd{end+1} = [fT1w ':resample=cubic \'];
%     end
%     if ~isempty(fB1)
%         cmd{end+1} = [fB1 ':resample=cubic:visible=0 \'];
%     end
%     cmd{end+1} = [fRun(I).fResp ':resample=cubic:visible=0 \'];
%     if ~isempty(fSatinIndexPlus_echoRms)
%         cmd{end+1} = [fSatinIndexPlus_echoRms ':resample=cubic:visible=0 \'];
%     end
%     thresh = abs(norminv(0.95));
%     if ~isempty(fSatinIndex_echoRms)
%         cmd{end+1} = [fSatinIndex_echoRms ':colormap=heat:resample=cubic:visible=0 \'];
%     end
%     cmd{end+1} = [fRun(I).fStat ':colormap=heat:resample=cubic:visible=0 &'];
%     cmd = strjoin(cmd,newline);
% 
%     fRun(I).cmdVisFs = cmd;
% end
% 
% 
% %%% full session
% if isempty(fSes)
%     fSes = fRun;
% end
% if isempty(fSes_echoRms)
%     fSes_echoRms = fRun_echoRms;
% end
% 
% if nEcho>1
%     if isempty(fSes_echoRms)
%         fSes_echoRms = fRun_echoRms;
%     end
%     %%%% Multi echo
%     %%%%% catenate echo
%     fSes_echoCat = fSes(:,1);
%     fSes_echoCat.fResp = {fSes.fResp};
%     fSes_echoCat.fBase = {fSes.fBase};
%     fSes_echoCat.fRespOnBase = {fSes.fRespOnBase};
%     fSes_echoCat.fFit = {fSes.fFit};
%     fSes_echoCat.fResid = {fSes.fResid};
% 
%     fOut = strsplit(fSes(1).fStat,'_'); fOut{contains(fOut,'echo-')} = 'echo-cat';
%     fOut = replace(strjoin(fOut,'_'),'_stats.nii.gz','_fullF.nii.gz'); if ~exist(fileparts(fOut),'dir'); mkdir(fileparts(fOut)); end
% 
%     cmd = {srcAfni};
%     cmd{end+1} = '3dbucket -overwrite \';
%     cmd{end+1} = ['-prefix ' fOut ' \'];
%     cmd{end+1} = [strjoin({fSes.fStat},'[0] ') '[0]'];
%     % cmd = {srcAfni};
%     % cmd{end+1} = '3dTcat -overwrite \';
%     % cmd{end+1} = ['-prefix ' fOut ' \'];
%     % cmd{end+1} = [strjoin({fSes.fStat},'[0] ') '[0]'];
%     cmd = strjoin(cmd,newline);
%     [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end
%     fSes_echoCat.cmd = cmd;
%     fSes_echoCat.fStat = fOut;
% 
%     cadidateField = {'fAvCatAvEchoCat' 'fAvEchoCat' 'fAvCatAv' 'fAv'};
%     cadidateField = cadidateField(ismember(cadidateField,fields(fFunc)));
%     fSes_echoCat.fUnder = fFunc.(cadidateField{1});
% 
%     %%%%% catenate cross-echo rms
%     fSes_echoCat.rms.fResp = fSes_echoRms.fResp;
%     fSes_echoCat.rms.fBase = fSes_echoRms.fBase;
%     fSes_echoCat.rms.fRespOnBase = fSes_echoRms.fRespOnBase;
%     fSes_echoCat.rms.fRespOnBaseMovie = fSes_echoRms.fRespOnBaseMovie;
%     fSes_echoCat.rms.fRespOnBaseMovieHighBit = fSes_echoRms.fRespOnBaseMovieHighBit;
%     fSes_echoCat.rms.fFit = fSes_echoRms.fFit;
%     fSes_echoCat.rms.fResid = fSes_echoRms.fResid;
%     fSes_echoCat.rms.fStat = fSes_echoRms.fStat;
% 
%     cadidateField = {'fAvCatAvEchoRms' 'fAvEchoRms' 'fAvCatAv' 'fAv'};
%     cadidateField = cadidateField(ismember(cadidateField,fields(fFunc)));
%     fSes_echoCat.rms.fUnder = fFunc.(cadidateField{1});
% 
% 
% 
% 
% 
%     % %%%%% visulaize with afni
%     % cmd = {srcAfni};
%     % if isfield(param,'layout') && ~isempty(param.layout) && exist(param.layout,'file')
%     %     cmd{end+1} = ['afni -layout ' param.layout ' \'];
%     % else
%     %     cmd{end+1} = 'afni \';
%     % end
%     % cmd{end+1} = ['-tbar ' num2str(size(fRun,1)) 'runs \'];
%     % cmd{end+1} = [fSes_echoCat.fUnder ' \'];
%     % cmd{end+1} = [fSes_echoCat.rms.fUnder ' \'];
%     % cmd{end+1} = [fSes_echoCat.fStat ' \'];
%     % cmd{end+1} = [fSes_echoCat.rms.fStat ' \'];
%     % cmd{end+1} = [strjoin([fSes_echoCat.fResp fSes_echoCat.rms.fResp],' ') ' &'];
%     % cmd = strjoin(cmd,newline); % disp(cmd)
%     % fSes_echoCat.cmdVisAfni = cmd;
%     %
%     % if verbose; disp(cmd); end
%     % if verbose>1; [status,cmdout] = system(cmd); end
% 
%     %%%%% visulaize with freeview
%     cmd = {srcFs};
%     cmd{end+1} = 'freeview \';
%     cmd{end+1} = '-timecourse \';
%     % cmd{end+1} = '-subtitle \';
%     cmd{end+1} = [char(fSes_echoCat.fUnder) ':name=echo-cat_funcAv:visible=0 \'];
%     cmd{end+1} = [char(fSes_echoCat.rms.fUnder) ':name=echo-rms_funcAv:visible=1 \'];
%     if ~isempty(fT1w)
%         cmd{end+1} = [fT1w ':resample=cubic:name=T1w:visible=0 \'];
%     end
%     if ~isempty(fB1)
%         cmd{end+1} = [fB1 ':colormap=turbo:resample=cubic:name=B1map:colorscale=50,130:visible=0 \'];
%     end
%     for E = 1:length(fSes_echoCat.fResp)
%         tmp = strsplit(fSes_echoCat.fResp{E},'_');
%         cmd{end+1} = [fSes_echoCat.fResp{E} ':name=' tmp{contains(tmp,'echo-')} '_resp:visible=0 \'];
%         cmd{end+1} = [fSes_echoCat.fBase{E} ':name=' tmp{contains(tmp,'echo-')} '_base:visible=0 \'];
%         cmd{end+1} = [fSes_echoCat.fRespOnBase{E} ':name=' tmp{contains(tmp,'echo-')} '_respOnBase:visible=0 \'];
%     end
%     cmd{end+1} = [fSes_echoCat.rms.fResp ':name=echo-rms_resp:visible=0 \'];
%     cmd{end+1} = [fSes_echoCat.rms.fBase ':name=echo-rms_base:visible=0 \'];
%     cmd{end+1} = [fSes_echoCat.rms.fRespOnBase ':name=echo-rms_respOnBase:visible=0 \'];
%     if ~isempty(fSatinIndexPlus_echoRms)
%         cmd{end+1} = [fSatinIndexPlus_echoRms ':resample=cubic:name=echo-rms_satIndexNumPlusBack:visible=0 \'];
%     end
%     if ~isempty(fSatinIndexPlus_echoCat)
%         cmd{end+1} = [fSatinIndexPlus_echoCat ':resample=cubic:name=echo-cat_satIndexNumPlusBack:visible=0 \'];
%     end
%     if ~isempty(fSatinIndexAbs_echoRms)
%         cmd{end+1} = [fSatinIndexAbs_echoRms ':resample=cubic:name=echo-rms_satIndexNumAbs:visible=0 \'];
%     end
%     if ~isempty(fSatinIndexAbs_echoCat)
%         cmd{end+1} = [fSatinIndexAbs_echoCat ':resample=cubic:name=echo-cat_satIndexNumAbs:visible=0 \'];
%     end
%     if ~isempty(fSatinIndex_echoRms)
%         cmd{end+1} = [fSatinIndex_echoRms ':colormap=heat:resample=cubic:name=echo-rms_satIndex:visible=0 \'];
%     end
%     if ~isempty(fSatinIndex_echoCat)
%         cmd{end+1} = [fSatinIndex_echoCat ':colormap=heat:resample=cubic:name=echo-cat_satIndex:visible=0 \'];
%     end
% 
%     for E = 1:size(fSes,2)
%         cmdTmp = {srcAfni};
%         cmdTmp{end+1} = ['fdrval -qinput ' fSes(E).fStat ' 0 0.05'];
%         [~,tmp] = system(strjoin(cmdTmp,newline));
%         tmp = strsplit(tmp,newline);
%         thresh(E) = str2num(tmp{end-1});
% 
%         cmdTmp = {srcAfni};
%         cmdTmp{end+1} = ['3dBrickStat ' fSes(E).fStat];
%         [~,tmp] = system(strjoin(cmdTmp,newline));
%         tmp = strsplit(tmp,newline);
%         maxF(E) = str2num(tmp{end-1});
%         if maxF(E)>20; maxF(E) = 20; end
%         cmd{end+1} = [fSes(E).fStat ':colormap=heat:heatscale=' num2str(thresh(E)) ',' num2str(maxF(E)) ':name=echo-' num2str(E) '_fullF:visible=0 \'];
%     end
%     thresh = mean(thresh);
%     maxF = min(maxF); if maxF>20; maxF = 20; end
%     cmd{end+1} = [fSes_echoCat.fStat ':colormap=heat:heatscale=' num2str(thresh) ',' num2str(maxF) ':name=echo-cat_fullF:visible=0 \'];
% 
%     cmdTmp = {srcAfni};
%     cmdTmp{end+1} = ['fdrval -qinput ' fSes_echoCat.rms.fStat ' 0 0.05'];
%     [~,tmp] = system(strjoin(cmdTmp,newline));
%     tmp = strsplit(tmp,newline);
%     thresh = str2num(tmp{end-1});
% 
%     cmdTmp = {srcAfni};
%     cmdTmp{end+1} = ['3dBrickStat ' fSes_echoCat.rms.fStat];
%     [~,tmp] = system(strjoin(cmdTmp,newline));
%     tmp = strsplit(tmp,newline);
%     maxF = str2num(tmp{end-1});
%     maxF = min(maxF); if maxF>20; maxF = 20; end
%     cmd{end+1} = [fSes_echoCat.rms.fStat ':colormap=heat:heatscale=' num2str(mean(thresh)) ',' num2str(maxF) ':name=echo-rms_fullF &'];
% 
%     cmd = strjoin(cmd,newline);
%     fSes_echoCat.cmdVisFs = cmd;
% else
%     %%%% single-echo
%     fSes_echoCat = [];
%     E = 1;
% 
%     % %%%%% visualize with afni
%     % fSes(1,E).fUnder = fFunc.fAvCatAv{1,E};
%     % cmd = {srcAfni};
%     % if isfield(param,'layout') && ~isempty(param.layout) && exist(param.layout,'file')
%     %     cmd{end+1} = ['afni -layout ' param.layout ' \'];
%     % else
%     %     cmd{end+1} = 'afni \';
%     % end
%     % cmd{end+1} = ['-tbar ' num2str(size(fRun,1)) 'runs \'];
%     % cmd{end+1} = [fSes(1,E).fUnder ' \'];
%     % cmd{end+1} = [fSes(1,E).fStat ' \'];
%     % cmd{end+1} = [fSes(1,E).fResp ' &'];
%     % cmd = strjoin(cmd,newline); %disp(cmd)
%     % fSes(1,E).cmdVisAfni = cmd;
%     %
%     % if verbose; disp(cmd); end
%     % if verbose>1; [status,cmdout] = system(cmd); if status || isempty(cmdout); dbstack; error(cmdout); error('x'); end; end
% 
%     %%%%% visualize with freeview
%     cmd = {srcFs};
%     cmd{end+1} = 'freeview \';
%     cmd{end+1} = [fSes(1,E).fBase ' \'];
%     if ~isempty(fT1w)
%         cmd{end+1} = [fT1w ':resample=cubic \'];
%     end
%     if ~isempty(fB1)
%         cmd{end+1} = [fB1 ':resample=cubic:visible=0 \'];
%     end
%     cmd{end+1} = [fSes(1,E).fResp ':resample=cubic:visible=0 \'];
%     if ~isempty(fSatinIndexPlus_echoRms)
%         cmd{end+1} = [fSatinIndexPlus_echoRms ':resample=cubic:visible=0 \'];
%     end
%     thresh = abs(norminv(0.95));
%     if ~isempty(fSatinIndex_echoRms)
%         cmd{end+1} = [fSatinIndex_echoRms ':colormap=heat:resample=cubic:visible=0 \'];
%     end
%     cmd{end+1} = [fSes(1,E).fStat ':colormap=heat:resample=cubic &'];
%     cmd = strjoin(cmd,newline);
% 
%     fSes(1,E).cmdVisFs = cmd;
% end






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

            fStat  = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_model-' HRmodel '_stats.nii.gz']);
            fFit   = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_model-' HRmodel '_fit.nii.gz'  ]);
            fResid = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_model-' HRmodel '_resid.nii.gz']);
            fMask  = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_model-' HRmodel '_mask.nii.gz' ]);
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
                fMat    = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_model-' HRmodel '_stats.xmat.1D']);
                fMatFig = fullfile(fileparts(replace(fOut,'.nii.gz','')),['task-' param.dsgn.task '_model-' HRmodel '_stats.xmat.fig']);
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
                    [cmdTmpTmp,param.dsgn.nReg] = afniCmd2(fIn,fStim,fMask,curParam,fResp,fRespStd,fFit,fResid,fMat,fStat,verbose,param.dryRun);
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


function [cmd,nReg] = afniCmd2(fIn,fStim,fMask,param,fResp,fRespStd,fFit,fResid,fMat,fStat,verbose,dryRun)
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









function [cmd,nReg] = afniCmd(fIn,fMask,fStim,nDummy,tr,startSeq,durSeq,condSeq,HRmodel,label,param,fResp,fFit,fResid,fMat,fStat,verbose,nDummyRemoved,trDecon,dryRun,nFrame)
if ~exist('nDummyRemoved','var'); nDummyRemoved = []; end
if isempty(nDummyRemoved);        nDummyRemoved = nDummy; warning('param.nDummyRemoved not specified, assuming param.nDummyRemoved = param.nDummy'); end
if nDummyRemoved && nDummyRemoved~=nDummy; dbstack; error('param.nDummyRemoved specified, but does not match param.nDummy'); end
if ~exist('nFrame','var'); nFrame = []; end
if ~exist('param','var');   param = []; end
cmd = {'3dDeconvolve -overwrite \'};
% cmdTmp{end+1} = ['-force_TR ' num2str(trStim) ' \'];
if ~dryRun
    if iscell(fIn)
        if nDummyRemoved
            cmd{end+1} = ['-input ' strjoin(fIn,' ') ' \'];
        else
            cmd{end+1} = ['-input ' strjoin(fIn,' ') '[' num2str(nDummy) '..$] \'];
        end
    else
        if nDummyRemoved
            cmd{end+1} = ['-input ' fIn ' \'];
        else
            cmd{end+1} = ['-input ' fIn '[' num2str(nDummy) '..$] \'];
        end
    end
    if ~isempty(fMask)
        cmd{end+1} = ['-mask ' fMask ' \'];
    end
else
    if isempty(nFrame)
        nFrame = MRIread(fIn{1},1); nFrame = nFrame.nframes;
    end
    cmd{end+1} = ['-nodata ' num2str(nFrame) ' ' num2str(tr,'%0.16f') ' \'];
end
cmd{end+1} = '-polort A \';
% if ~exist('trMri','var') || isempty(trMri) || ~exist('nFrame','var') || isempty(nFrame)
if isempty(tr)
    mri = MRIread(fIn{1},1);
    tr = mri.tr/1000;
end

% nFrame = mri.nframes;
% end
% trMri = 3;
cmd{end+1} = ['-stim_times_subtract ' num2str(tr*nDummy,'%f') ' \'];

if isempty(param)
    % Old way
    cmd{end+1} = '-num_stimts 1 \';
    cmd{end+1} = ['-stim_label 1 ' char(label) ' \'];

    % write design to file
    fido = fopen(fStim, 'w');
    if iscell(fIn)
        for i = 1:length(fIn)
            fprintf(fido,'%.3f ',startSeq(condSeq==1));
            fprintf(fido,'\n');
        end
    else
        fprintf(fido,'%.3f ',startSeq(condSeq==1));
        fprintf(fido,'\n');
    end
    fclose(fido);


    k = 1;

    switch HRmodel
        case 'SPMG2'
            nReg = 2;
            if max(abs(diff(durSeq)))/max(durSeq) > 0.0001; dbstack; error('stim duration cannot be different across trials'); end
            cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' ''' HRmodel '(' num2str(mean(durSeq),'%0.3f') ')'' \'];
        case 'SPMG3'
            nReg = 3;
            if max(abs(diff(durSeq)))/max(durSeq) > 0.0001; dbstack; error('stim duration cannot be different across trials'); end
            cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' ''' HRmodel '(' num2str(mean(durSeq),'%0.3f') ')'' \'];
        case 'TENTzero'
            deconWin = min(diff(startSeq(condSeq==1)));
            if (deconWin/trDecon)/ceil(deconWin/trDecon)>0.9
                deconWin = ceil(deconWin/trDecon)*trDecon;
            else
                deconWin = floor(deconWin/trDecon)*trDecon;
            end
            b = 0;
            % c = deconWin-trDecon;
            c = round((deconWin-trDecon)/trDecon)*trDecon;
            % (c-b)/(n-1)
            nReg = round( (c-b)/trDecon + 1 );
            cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' ''TENTzero(' num2str(b) ',' num2str(c) ',' num2str(nReg) ')'' \'];
            nReg = nReg - 2;
            cmd{end+1} = ['-TR_times ' num2str(trDecon,'%f') ' \'];
            if ~dryRun
                cmd{end+1} = ['-iresp ' num2str(k) ' ' fResp ' \'];
            end
        otherwise
            dbstak; error('X');
    end

else
    % New way (only implemented for dry runs so far)
    stimCondList  = sort(unique(condSeq))';
    stimCondLabel = replace(cellstr(num2str(stimCondList,'Trial%i')),' ','');
    nRegAll = [];
    cmd{end+1} = ['-num_stimts ' num2str(length(stimCondList)) ' \'];
    for k = 1:length(stimCondList)
        cmd{end+1} = ['-stim_label ' num2str(k) ' ' [char(label) stimCondLabel{k}] ' \'];

        % write design to file
        if dryRun
            fStim = [tempname '_startTime.1D' ];
        else
            dbstack; error('code that')
        end
        fido = fopen(fStim, 'w');
        if iscell(fIn)
            for i = 1:length(fIn)
                fprintf(fido,'%.3f ',startSeq(stimCondList(k)==condSeq));
                fprintf(fido,'\n');
            end
        else
            fprintf(fido,'%.3f ',startSeq(stimCondList(k)==condSeq));
            fprintf(fido,'\n');
        end
        fclose(fido);


        % k = 1;

        switch HRmodel
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
                deconWin = min(diff(startSeq));
                % deconWin = min(diff(startSeq(condSeq==1)));
                if (deconWin/trDecon)/ceil(deconWin/trDecon)>0.9
                    deconWin = ceil(deconWin/trDecon)*trDecon;
                else
                    deconWin = floor(deconWin/trDecon)*trDecon;
                end
                b = 0;
                % c = deconWin-trDecon;
                c = round((deconWin-trDecon)/trDecon)*trDecon;
                % (c-b)/(n-1)
                nReg = round( (c-b)/trDecon + 1 );
                cmd{end+1} = ['-stim_times ' num2str(k) ' ' fStim ' ''TENTzero(' num2str(b) ',' num2str(c) ',' num2str(nReg) ')'' \'];
                nReg = nReg - 2;
                if ~dryRun
                    cmd{end+1} = ['-iresp ' num2str(k) ' ' fResp ' \'];
                end
                nRegAll(k) = nReg;
            otherwise
                dbstak; error('X');
        end
    end
    nReg = nRegAll;
    cmd{end+1} = ['-TR_times ' num2str(trDecon,'%f') ' \'];

    if dryRun
        % fMat = replace(fMat,'_stats.xmat.1D','_stats.xmatPerTrial.1D');
        % fMat  = [tempname '_stats.xmat.1D'];
    else
        dbstack; error('code that')
    end
end




if ~dryRun
    cmd{end+1} = ['-fitts ' fFit ' \'];
    cmd{end+1} = ['-errts ' fResid ' \'];
    cmd{end+1} = '-bout \';
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


