function addWinAndW(H,gram)
if isempty(H)
    H = gca;
else
    axes(H)
end

addWin([],gram)
addW([],gram)



f = permute(gram.f,[1 5 2 3 4 6 7 8]);;
W = permute(gram.W,[1 5 2 3 4 6 7 8]); if length(unique(W))==1; W = unique(W); end
K = permute(gram.K,[1 5 2 3 4 6 7 8]); if length(unique(K))==1; K = unique(K); end
T = permute(gram.T,[1 5 2 3 4 6 7 8]); if length(unique(T))==1; T = unique(T); end
winSz = gram.lWin;
winStep = gram.param.win(2);
if range(diff(f)) / mean(f) > 1e-15; error('variable dx'); end
df = mean(diff(f));

% x = mean(xlim); x = round(x/winStep)*winStep;
% x = x .* [1 1]; xx = x + [-0.5 0.5].*winSz;
% y  = mean(f); y = round(y/df)*df;
% y = y .* [1 1]; yy = y + [-1   1  ].*W;
im = findobj(H.Children,'type','image');
x = im.XData(1) - mean(diff(im.XData))/2;
x = x .* [1 1];
xx = x + [0 1].*winSz;
y  = im.YData(end) + mean(diff(im.YData))/2;
y = y .* [1 1];
yy = y - [0 2].*W;
line(x,yy,'color','m','linewidth',3);
line(xx,y,'color','m','linewidth',3);
