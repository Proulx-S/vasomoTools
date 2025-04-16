function [xMat,hMat] = plotDsgnMat(fMat,force,verbose)
global src
if ~exist('verbose','var'); verbose = []; end
if isempty(verbose);        verbose = 0 ; end
if ~exist('force','var');   force = []; end
if isempty(force);          force = 0 ; end
param = fMat.param;
if ~isfield(param,'nDummyRemoved') || isempty(param.nDummyRemoved)
    param.nDummyRemoved = param.nFrameOrig-param.nFrame;
end
if diff(param.nDummyRemoved)>0; dbstack; error('nDummyRemoved cannot be different across runs'); end

%% Extract design matrix
cmdX = {src.afni};
cmdX{end+1} = ['1dcat ' char(fMat(1).fMat)];
[~,cmdout] = system(strjoin(cmdX,newline));
mat = str2num(cmdout);
nReg  = param.dsgn.nReg;
nPoly = size(mat,2) - sum(nReg);
tStim = 0:(size(mat,2)-1);
tStim = tStim-nPoly;
switch param.model
    case 'TENTzero'
        tStim = tStim + 1;
        tStim = tStim .* param.trDecon;
    case 'TENT'
        dbstack; error('double-check that')
        % nPoly = size(mat,2) - param.funDsgn.nReg;
        % tStim = (0:param.funDsgn.nReg-1).*param.trDecon;
        % tStim = [linspace(tStim(1)-(param.trDecon.*nPoly),tStim(1)-param.trDecon,nPoly) tStim];
    case 'SPMG2'
        % tStim = tStim + 1;
    otherwise
        dbstack; error('code that');
end


% yaxis time after run onset
for r = 1:length(param.nFrame)
    iRun(:,r) = [1 param.nFrame(r)-param.nDummyIgnore];
    if r == 1
        iSes(:,r) = iRun(:,r);
    else
        iSes(:,r) = iRun(:,r) + sum(param.nFrame(1:r-1));
    end
    try
        tRun(:,r) = (iRun(:,r) + param.nDummyRemoved(r) + param.nDummyIgnore -1) .* param.tr(r);
    catch
        tRun(:,r) = (iRun(:,r) + param.nDummyRemoved(r) + param.nDummyIgnore -1) .* param.tr;
        warning(['only one tr found in param.tr' newline 'using the same for all runs'])
    end
end

if param.PCflag
    iRun = cat(2,iRun,iRun);
    iSes = cat(2,iSes,iSes + iSes(end));
    tRun = cat(2,tRun,tRun);
end


% xaxis time after stim onset and baselines
switch param.model
    case {'TENTzero' 'TENT'}
        for k = 1:length(param.dsgn.condLabel)
            iStim{k} = 1:param.dsgn.nReg(k);
            if k == 1
                iStims{k} = iStim{k};
            else
                iStims{k} = iStim{k} + iStims{k-1}(end);
            end
        end
        switch param.model
            case 'TENTzero'
                tStim = [iStim{:}].*param.trDecon;
            case 'TENT'
                tStim = ([iStim{:}]-1).*param.trDecon;
        end
        iStim = [iStims{:}]; clear iStims
    case 'SPMG2'
        iStim = tStim(tStim>=0)+1;
    otherwise
        error('code that')
end
pInd = false([1 size(mat,2)]);
pInd(1:nPoly) = true;
iStim = iStim+nnz(pInd);


%% Plot design matrix
if force || ~exist(char(fMat.fMatFig),'file')
    hMat = figure('Visible','off');
    h = imagesc(mat); colormap gray
    h.Parent.YTick = iSes(:);
    h.Parent.YTickLabel = cellstr(num2str(tRun(:),'%0.3f'));
    ylabel('time after run onset (s)');
    h.Parent.XTick = iStim;
    switch param.model
        case {'TENTzero' 'TENT'}
            h.Parent.XTickLabel = cellstr(num2str(tStim','%0.3f'));
            h.Parent.XTickLabelRotation = 45;
            xlabel('time after stim onset (s)');
        case 'SPMG2'
            h.Parent.XTickLabel = {};
            xlabel('regressors')
        otherwise
            error('code that')
    end
    clim([-1 1])

    [~,b,~] = fileparts(fileparts(fMat.fStat));
    title(b,'interpreter','none')
    set(hMat, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    savefig(hMat,char(fMat.fMatFig),'compact')
    if verbose>0
        hMat.Visible = 'on';
        hMat.WindowStyle = 'docked';
    else
        close(hMat)
    end
    drawnow
end

%% Output design matrix
xMat.mat   = mat;
xMat.tRun  = tRun';
xMat.tStim = tStim;
xMat.nReg  = nReg;
xMat.nPoly = nPoly;



