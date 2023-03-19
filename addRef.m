function addRef(w,ax,winSz)

if ~exist('ax','var') || isempty(ax)
    ax = gca;
end
hold on
if ~exist('winSz','var') || isempty(winSz)
    h = findobj(ax.Children,'Tag','shadedErrorBar_mainLine');
    f = h.XData;
    f0w = exp(linspace(log(f(2)),log(f(end)),6)); f0w([1 end]) = [];
    y = ax.YLim(1).*[1 1];
    for ii = 1:length(f0w)
        x = f0w(ii)+[-1 1].*w;
        if x(1)>f(1) && x(2)<f(end)
            plot(x,y,'g','LineWidth',3)
        end
    end
else
    h = findobj(ax.Children,'Type','Image');
    if isempty(h)
        h = findobj(ax.Children,'Type','Surface');
        z = zlim; z = z(2);
    else
        z = 0;
    end
    f = h.XData;
    f0w = exp(linspace(log(f(2)),log(f(end)),6)); f0w([1 end]) = [];
    t = h.YData;
    y = [t(1) t(1)+winSz];
    for ii = 1:length(f0w)
        x = f0w(ii)+[-1 1].*w;
        if x(1)>f(1) && x(2)<f(end)
            hPlot = plot(x,[1 1].*y(1),'g','LineWidth',3);
            hPlot.ZData = [1 1].*z;
        end
    end
    xLim = xlim;
    hPlot = plot([1 1].*xLim(1),y,'g','LineWidth',3);
    hPlot.ZData = [1 1].*z;
end