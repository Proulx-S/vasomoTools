function [fMat,hMat,xMat] = plotDsgnMat(fMat,volTs,verbose,saveFalg)
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
    case {'SPMG2' 'SPMG3'}
        tStim = tStim + 1;
end
% tRun = 0:size(mat,1)-1;

h = imagesc(mat); colormap gray
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



pInd = false([1 size(mat,2)]);
pInd(1:nPoly) = true;
param.dsgn.condLabel
h.Parent.XTick = tStim;
ax.XTickLabel(pInd) = cell(1,nnz(pInd));
switch param.model
    case {'TENTzero' 'TENT'}
        ax.XTickLabel(~pInd) = cellstr(num2str(tStim(~pInd)','%0.3f'))';
        xlabel('time after stim onset (s)')
    case {'SPMG2' 'SPMG3'}
        ax.XTickLabel(~pInd) = cellstr(num2str(tStim(~pInd)','%i'))';
        xlabel('regressors index')
end
clim([-1 1])
[~,b,~] = fileparts(fileparts(fMat.fIn(1:length(fMat))));
title(b,'interpreter','none')
set(hMat, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
if saveFalg
    savefig(hMat,fMat(1).fMatFig,'compact')
end
if verbose>0
    hMat.Visible = 'on';
    hMat.WindowStyle = 'docked';
else
    close(hMat)
end
drawnow


[fMat.fMatFig] = deal(fMat(1).fMatFig);


xMat.mat   = mat;
xMat.tRun  = tRun';
xMat.tStim = tStim;
xMat.nReg  = nReg;
xMat.nPoly = nPoly;



