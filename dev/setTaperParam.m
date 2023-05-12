function [param,w] = setTaperParam(W,TW,funTs,winSz)

tr = funTs.tr/1000;
if ~exist('winSz','var') || isempty(winSz)
    T = funTs.nframes*tr;
else
    T = winSz;
end

if ~isempty(W) && isempty(TW)
    TW = T*W;
    K = nan;
    k = floor(2*TW-1);
    tw = (k+1)/2;
    w = tw/T;
elseif isempty(W) && ~isempty(TW)
else
    error('Define either W or TW and leave the other one empty')
end

if k<2
    minK=2; minTW = (minK+1)/2;
    error(['too few tapers. Try increasing W or setting winSz to at least ' num2str(minTW/W)]);
end

disp(table([W TW K]',[w tw k]','VariableNames',{'requested' 'obtained'},'RowNames',{'W' 'TW','K'}));
param.tapers = [tw k];
param.Fs = 1/tr;



        