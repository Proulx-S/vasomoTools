function addOndur(H,onsets,ondurs)
if isempty(H)
    H = gca;
else
    axes(H)
end
hold on

if contains(H.Title.String,'trial-locked')
    xx = onsets(1).*[1 1];
else
    % x = onsets.*[1 1];
    xx = [onsets onsets+ondurs];
end
yy = ylim;

hPatch = cell(size(xx,1),1);
for i = 1:size(xx,1)
    hPatch{i,1} = patch(xx(i,[1 2 2 1 1]),yy([1 1 2 2 1]),'r','EdgeColor','none','FaceAlpha',0.1);
    % line(xx(i,:),yy,'Color','r','LineStyle','--')
end
uistack([hPatch{:}],'bottom');
