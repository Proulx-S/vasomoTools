function addW(H,gram)
if isempty(H)
    H = gca;
else
    axes(H)
end
hold on

[~,W] = K2W(gram.T,gram.K,0);
% W = unique(gram.W);

if contains(H.Title.String,'spectrum')
    x  = mean(xlim); x = x.*[1 1];
    xx = x + [-1 1].*W;
    y = ylim; y = y(1).*[1 1];
    line(xx,y,'Color','m','linewidth',3)
elseif contains(H.Title.String,'spectrogram') || contains(H.Title.String,'coherogram')
    im = findobj(H.Children,'type','image');
    x = im.XData(1) - mean(diff(im.XData))/2;
    x = x .* [1 1];
    y  = im.YData(end) + mean(diff(im.YData))/2;
    y  = y .* [1 1];
    yy = y - [0 2].*W;
    line(x,yy,'color','m','linewidth',3);
else
    dbstack; error('code that')
end
