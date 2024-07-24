function addWin(H,gram,x,y)
if isempty(H)
    H = gca;
else
    axes(H)
end
if ~exist('x','var'); x = []; end
if ~exist('y','var'); y = []; end

hold on

% winSz = gram.lWin;
Fs = gram.param.Fs;
winSz = gram.win(1) + 0.5/Fs;
if contains(H.Title.String,'timeseries')
    if size(gram.t,2)>1
        y = ylim;
        xx = gram.t(:,:,:,:,:,:,1);
        yy = repmat(y(1),[size(xx,1) 1]);
        line(xx,yy,'Color','m','linewidth',3)

        xx = gram.t(:,:,:,:,:,:,end);
        yy = repmat(y(2),[size(xx,1) 1]);
        line(xx,yy,'Color','m','linewidth',3)
    else
        yLim = ylim;
        xLim = xlim;
        if isempty(x)
            x = xLim(1); x = x.*[1 1];
            xx = x + [0 1].*winSz;
        else
            x = x.*[1 1];
            xx = x + [-0.5 0.5].*winSz;
        end
        if isempty(y) || y==-inf
            y = yLim(1);
        elseif y==inf
            y = yLim(2);
        end
        y = y.*[1 1];
        line(xx,y,'Color','m','linewidth',3)
    end

elseif contains(H.Title.String,'spectrogram') || contains(H.Title.String,'coherogram')
    im = findobj(H.Children,'type','image');
    if isempty(x)
        x = im.XData(1) - mean(diff(im.XData))/2;
        x = x .* [1 1];
        xx = x + [0 1].*winSz;
    else
        x = x .* [1 1];
        xx = x + [-0.5 0.5].*winSz;
    end
    if isempty(y)
        y  = im.YData(end) + mean(diff(im.YData))/2;
    else
        dbstack; error('code that');
    end
    y = y .* [1 1];
    line(xx,y,'color','m','linewidth',3);
else
    dbstack; error('code that')
end
