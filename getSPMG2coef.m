function [amp,hFig] = getSPMG2coef(f,mask,plotFlag)
if ~exist('plotFlag','var'); plotFlag = []; end
if isempty(plotFlag);        plotFlag = 0 ; end

coef = MRIread(f); coef = coef.vol;
mask = MRIread(mask); mask = logical(mask.vol);
amp = zeros(size(coef,[1 2 3]));

coef = permute(coef,[4 1 2 3]);
coef2 = coef(:,mask);
slp = coef2(1,:)'\coef2(2,:)';
slpPos = coef2(1,coef2(1,:)>0)'\coef2(2,coef2(1,:)>0)';
slpNeg = coef2(1,coef2(1,:)<0)'\coef2(2,coef2(1,:)<0)';
coefMod = [1 slp];

% [~,v] = min(sqrt(sum((coef2 - [27.681;-4.966]).^2,1)));
% coefMod * coef2(:,v)
% scatter(coef2(1,v),coef2(2,v))

amp(:) = coefMod * coef(:,:);

if plotFlag
    hFig = figure('WindowStyle','docked');
    scatter(coef2(1,:),coef2(2,:),'MarkerEdgeColor','k'); hold on
    grid on
    hRef = refline(slp,0); hRef.Color = 'k';
    xLim = xlim;
    yLim = ylim;
    xx = [0 xLim(2)]; yy = xx.*slpPos;
    line(xx,yy,'color','r')
    xx = [xLim(1) 0]; yy = xx.*slpNeg;
    line(xx,yy,'color','b')
else
    hFig = [];
end
