function addFreq(H,onsets,ondurs)
if ~exist('ondurs','var'); ondurs = []; end
if isempty(H)
    H = gca;
else
    axes(H)
end
hold on

fStim = 1/mean(diff(onsets));
if any(contains(H.Title.String,{'spectrogram' 'coherogram'}))
    yline(fStim,'Color','r','linestyle','--')
    xLim = xlim;
    text(xLim(2),fStim,'fStim','HorizontalAlignment','right','VerticalAlignment','baseline','Color','r')
elseif any(contains(H.Title.String,{'spectrum'}))
    xline(fStim,'Color','r','linestyle','--')
    yLim = ylim;
    text(fStim,yLim(2),'fStim','HorizontalAlignment','right','VerticalAlignment','top','Color','r')
end
if ~isempty(ondurs) && length(unique(ondurs))==1
    fOn = 1/ondurs(1);
    if any(contains(H.Title.String,{'spectrogram' 'coherogram'}))
        yline(fOn,'Color','g','linestyle',':')
        text(xLim(2),fOn,'fOn','HorizontalAlignment','right','VerticalAlignment','top','Color','g')
        yline(fOn/2,'Color','g','linestyle',':')
        text(xLim(2),fOn/2,'fOn','HorizontalAlignment','right','VerticalAlignment','top','Color','g')
    elseif any(contains(H.Title.String,{'spectrum'}))
        xline(fOn,'Color','g','linestyle',':')
        text(fOn,yLim(2),'fOn','HorizontalAlignment','left','VerticalAlignment','top','Color','g')
        xline(fOn/2,'Color','g','linestyle',':')
        text(fOn/2,yLim(2),'fOn','HorizontalAlignment','left','VerticalAlignment','top','Color','g')
    end
end
if ~isempty(ondurs)
    offsetList = onsets + ondurs;
    offdurList = onsets(2:end) - offsetList(1:end-1);
    if length(unique(offdurList))==1
        fOff = 1/offdurList(1);
        if any(contains(H.Title.String,{'spectrogram' 'coherogram'}))
            yline(fOff,'Color','b','linestyle',':')
            text(xLim(2),fOff,'fOff','HorizontalAlignment','right','VerticalAlignment','top','Color','b')
        elseif any(contains(H.Title.String,{'spectrum'}))
            xline(fOff,'Color','b','linestyle',':')
            text(fOff,yLim(2),'fOff','HorizontalAlignment','left','VerticalAlignment','top','Color','b')
        end
    end
end

