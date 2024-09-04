function h = addFreq(H,onsets,ondurs,freqFlag)
if ~exist('ondurs','var');     ondurs = []; end
if ~exist('freqFlag','var'); freqFlag = []; end
if isempty(freqFlag);        freqFlag = 0; end
if isempty(H)
    H = gca;
else
    axes(H)
end
hold on

fStim = 1/mean(diff(onsets));
if any(contains(H.Title.String,{'spectrogram' 'coherogram'}))
    h = yline(fStim,'Color','r','linestyle','--');
    xLim = xlim;
    text(xLim(2),fStim,'fStim','HorizontalAlignment','right','VerticalAlignment','top','Color','r')
elseif any(contains(H.Title.String,{'spectrum'}))
    h = xline(fStim,'Color','r','linestyle','--');
    yLim = ylim;
    text(fStim,yLim(2),'fStim','HorizontalAlignment','left','VerticalAlignment','top','Color','r')
end
if freqFlag && ~isempty(ondurs) && length(unique(ondurs))==1
    fOn = 1/ondurs(1);
    if any(contains(H.Title.String,{'spectrogram' 'coherogram'}))
        h = yline(fOn,'Color','g','linestyle',':');
        text(xLim(2),fOn,'fOn','HorizontalAlignment','right','VerticalAlignment','top','Color','g')
        if freqFlag>1
            h = yline(fOn/2,'Color','g','linestyle',':');
            text(xLim(2),fOn/2,'fOn/2','HorizontalAlignment','right','VerticalAlignment','top','Color','g')
        end
    elseif any(contains(H.Title.String,{'spectrum'}))
        h = xline(fOn,'Color','g','linestyle',':');
        text(fOn,yLim(2),'fOn','HorizontalAlignment','left','VerticalAlignment','top','Color','g')
        if freqFlag>1
            h = xline(fOn/2,'Color','g','linestyle',':');
            text(fOn/2,yLim(2),'fOn/2','HorizontalAlignment','left','VerticalAlignment','top','Color','g')
        end
    end
end

if freqFlag>1 && ~isempty(ondurs)
    offsetList = onsets + ondurs;
    offdurList = onsets(2:end) - offsetList(1:end-1);
    fOff = 1/mean(offdurList);
    if any(contains(H.Title.String,{'spectrogram' 'coherogram'}))
        h = yline(fOff,'Color','b','linestyle',':');
        text(xLim(2),fOff,'fOff','HorizontalAlignment','right','VerticalAlignment','top','Color','b')
        h = yline(fOff,'Color','b','linestyle',':');
        text(xLim(2),fOff/2,'fOff/2','HorizontalAlignment','right','VerticalAlignment','top','Color','b')
    elseif any(contains(H.Title.String,{'spectrum'}))
        h = xline(fOff,'Color','b','linestyle',':');
        text(fOff,yLim(2),'fOff','HorizontalAlignment','left','VerticalAlignment','top','Color','b')
        h = xline(fOff/2,'Color','b','linestyle',':');
        text(fOff/2,yLim(2),'fOff/2','HorizontalAlignment','left','VerticalAlignment','top','Color','b')
    end
end

