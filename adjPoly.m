function adjPoly(hIO,label,color,keepFlag)
    if ~exist('label','var');       label = []; end
    if isempty(label);              label = 'original'; end % original dilate2
    if ~exist('color','var');       color = []; end
    if isempty(color);              color = 'w'; end
    if ~exist('keepFlag','var');    keepFlag = []; end
    if isempty(keepFlag);           keepFlag = false; end
    for i = 1:length(hIO)
        doIt(hIO(i),label,color,keepFlag);
        hAO(i) = hIO(i).Parent;
    end
    set(hAO,'Color','none');


    function doIt(hIO,label,color,keepFlag)
        hPoly = findobj(hIO.Parent.Children,'Type','poly');
        switch label
            case {'off' 'on'}
                set(hPoly,'Visible',label);
                return
            case {'original' 'peakVox' 'dilate1' 'dilate1p5' 'dilate2'}
                if keepFlag==-1
                    delete(hPoly);
                end
                if keepFlag
                    hPoly = plot(hIO.Parent,hIO.UserData.roi.poly(ismember(hIO.UserData.roi.polyLabel,label)));
                    hPoly.FaceColor = 'none';
                else
                    hPoly(1).Shape = hIO.UserData.roi.poly(ismember(hIO.UserData.roi.polyLabel,label));
                end
                hPoly(1).LineWidth = 2.5;
                hPoly(1).EdgeColor = color;
            otherwise
                error('Invalid label: %s',label);
        end
