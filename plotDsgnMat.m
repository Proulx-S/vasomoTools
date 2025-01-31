function [fMat,hMat,xMat] = plotDsgnMat(fMat,param,fVolTs,volTs,verbose,saveFalg)
global srcAfni
if ~exist('verbose','var'); verbose = []; end
if isempty(verbose);        verbose = 0 ; end
if ~exist('volTs','var'); volTs = []; end
if ~exist('saveFalg','var'); saveFalg = []; end
if isempty(saveFalg);        saveFalg = 1 ; end


% Extract and plot design matrix
hMat = figure('Visible','off');
% hMat = figure('WindowStyle','docked');
cmdX = {srcAfni};
cmdX{end+1} = ['1dcat ' fMat(1).fMat];
[~,cmdout] = system(strjoin(cmdX,newline));
mat = str2num(cmdout);
nReg  = param.funDsgn.nReg;
nPoly = size(mat,2) - nReg;
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
tRun = 0:size(mat,1)-1;

if exist('volTs','var') && ~isempty(volTs)
    tRun = tRun.*param.tr;
    tRunStart = 0;
    for r = 1:length(volTs)
        tRunStart(end+1) = tRunStart(end) + volTs(r).nFrame;
    end
    tRunStart(end) = [];
    tRunEnd   = tRunStart + [volTs.nFrame] - 1;
    tRunStart = tRunStart.*param.tr;
    tRunEnd   = tRunEnd  .*param.tr;

    tLabel = [repmat((param.nDummyRemoved - 1)*param.tr,size(tRunStart)) [volTs.nFrame]-1.*param.tr];
    [t,b]  = sort([tRunStart tRunEnd]);
    tLabel = tLabel(b);

    h = imagesc(tStim,tRun,mat); colormap gray
    h.Parent.YTick = t;
    h.Parent.YTickLabel = cellstr(num2str(tLabel','%0.3f'));
    ylabel('time after run onset (s)');
else
    imagesc(tStim,tRun,mat); colormap gray
    ylabel('volume index');
end




pInd = false([1 size(mat,2)]);
pInd(1:nPoly) = true;
ax = gca;
ax.XTick = tStim;
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
if ~isempty(fVolTs)
    [~,b,~] = fileparts(fileparts(fVolTs(1:length(fMat))));
    title(b,'interpreter','none')
end
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



