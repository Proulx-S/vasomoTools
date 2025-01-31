function h = addW(H,gram)
if isempty(H)
    H = gca;
else
    axes(H)
end
hold on
Hgram = findobj(H.Children,'type','image');

[~,W] = K2W(gram.T,gram.K,0);
% W = unique(gram.W);

if contains(H.Title.String,'spectrum') || isempty(Hgram)
    x  = mean(xlim); x = x.*[1 1];
    xx = x + [-1 1].*W;
    y = ylim; y = y(1).*[1 1];
    h = line(xx,y,'Color','m','linewidth',3);
elseif contains(H.Title.String,'spectrogram') || contains(H.Title.String,'coherogram') || sum(abs(Hgram.YData - round(Hgram.YData)))>0
    im = findobj(H.Children,'type','image');
    x = 0;
    % x = im.XData(1) - mean(diff(im.XData))/2;
    x = x .* [1 1];
    y  = im.YData(end) + mean(diff(im.YData))/2;
    y  = y .* [1 1];
    yy = y - [0 2].*W;
    h = line(x,yy,'color','m','linewidth',3);
else
    dbstack; error('code that')
end
