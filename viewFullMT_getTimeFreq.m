function viewFullMT_getTimeFreq(hF)

figure(hF)

hP = drawpoint('Label','1');

plot(hP.Position,'.')

ax = findobj(hF.Children,'Type','Axes');

hL = findobj(ax.Children,'Type','Line');
for i = 1:2
    hL(i).XData = hL(i).XData - mean(hL(i).XData) + hP.Position(1);
    hL(i).YData = hL(i).YData - mean(hL(i).YData) + hP.Position(2);
end
