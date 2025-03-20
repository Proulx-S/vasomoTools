function [xMat,hMat] = plotDsgnMat(fMat,verbose,saveFalg)
% Reads


global src
if ~exist('verbose','var'); verbose = []; end
if isempty(verbose);        verbose = 0 ; end
if ~exist('volTs','var');     volTs = []; end
if ~exist('saveFalg','var'); saveFalg = []; end
if isempty(saveFalg);        saveFalg = 1 ; end
param = fMat.param;

% Extract and plot design matrix
hMat = figure('Visible','off');
% hMat = figure('WindowStyle','docked');
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
% tRun = 0:size(mat,1)-1;


h = imagesc(mat); colormap gray

% yaxis time after run onset
for R = 1:length(param.nFrame)
    iRun(:,R) = [1 param.nFrame(R)-param.nDummyIgnore];
    if R == 1
        iSes(:,R) = iRun(:,R);
    else
        iSes(:,R) = iRun(:,R) + sum(param.nFrame(1:R-1));
    end
    tRun(:,R) = (iRun(:,R) + param.nDummyRemoved + param.nDummyIgnore -1) .* param.tr(R);
end
h.Parent.YTick = iSes(:);
h.Parent.YTickLabel = cellstr(num2str(tRun(:),'%0.3f'));
ylabel('time after run onset (s)');


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
if saveFalg
    savefig(hMat,fMat.fMatFig,'compact')
end
if verbose>0
    hMat.Visible = 'on';
    hMat.WindowStyle = 'docked';
else
    close(hMat)
end
drawnow

xMat.mat   = mat;
xMat.tRun  = tRun';
xMat.tStim = tStim;
xMat.nReg  = nReg;
xMat.nPoly = nPoly;



