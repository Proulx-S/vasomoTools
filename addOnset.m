function addOnset(H,onsets)
if isempty(H)
    H = gca;
else
    axes(H)
end
hold on

if contains(H.Title.String,'evoked')
    x = onsets(1).*[1 1];
else
    x = onsets.*[1 1];
end
yy = ylim;

for i = 1:size(x,1)
    line(x(i,:),yy,'Color','r','LineStyle','--')
end
