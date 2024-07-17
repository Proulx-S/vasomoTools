function [ax,F] = plotTs(volTs,H,onsetList,durList)
if ~exist('H','var');                 H = []; end
if isempty(H);                        H = figure('WindowStyle','docked'); end
if ~exist('onsetList','var'); onsetList = []; end
if ~exist('durList','var');     durList = []; end

switch class(H)
    case 'matlab.graphics.layout.TiledChartLayout'
        F = H.Parent;
    case 'matlab.ui.Figure'
        F = H;
    otherwise
end

figure(F);
ax = {};
ax{end+1} = nexttile;

%%% average time series
volTs = vol2vec(volTs);
plot(squeeze(volTs.t),squeeze(mean(volTs.vec,2)),'k')
grid on
axis tight
xlabel('t (sec)')
ylabel('MR signal (a.u.)')

T = volTs.nframes.*volTs.tr/1000;
xlim([0 T])

paramStr = {['T=' num2str(T,'%0.2f') 'sec']};
if ~isempty(durList)
    paramStr{end+1} = ['dur=' num2str(mean(durList)) 'sec'];
end
title(['spatially averaged timeseries (' strjoin(paramStr,'; ') ')'])


if isempty(onsetList) && isfield(volTs,'dsgn') && isfield(volTs.dsgn,'onsetList')
    onsetList = volTs.dsgn.onsetList;
end
if isempty(durList) && isfield(volTs,'dsgn') && isfield(volTs.dsgn,'durList')
    durList = volTs.dsgn.durList;
end

if ~isempty(onsetList) && isempty(durList)
    addOnset([],onsetList)
elseif ~isempty(durList)
    addOndur([],onsetList,durList)
end


ax = [ax{:}];