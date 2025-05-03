function adjIm(hIO,label)
    if ~exist('label','var');       label = []; end
    if isempty(label);              label = 'coef_Q_dilate2'; end % 'coef_Q_im' 'coef_Q_dilate2'
    % if ~exist('lowAlpha','var'); lowAlpha = []  ; end
    % if isempty(lowAlpha);        lowAlpha = 0; end

        
    for i = 1:length(hIO)
        doIt(hIO(i),label);
        hAO(i) = hIO(i).Parent;
    end
    set(hAO,'Color','none');


    function doIt(hIO,label)
        if strcmp(label,'off')
            hIO.AlphaData = 0;
            return;
        end

        hIO.CData
        hIO.UserData.roi

        % Get the transparency
        trans = getTrans(hIO.UserData.roi,label,lowAlpha);
        
        % Apply transparency
        hIO.AlphaData = trans.alphaVal;
        hIO.UserData.trans = trans;



    function trans = getTrans(roi,label,lowAlpha)
        trans.label = label;
        trans.lowAlpha  = lowAlpha;
        switch label
            case 'act_imQ'
                im = roi.im.actP.im;
            case 'resp_imQ'
                im = roi.im.respP.im;
        end
        trans.threshVal = zeros(size(im));
        trans.alphaVal  = ones( size(im));
        trans.threshVal(:) = mafdr(im(:),'BHFDR',true);
        trans.alphaVal(trans.threshVal<=0.05) = 1;
        trans.alphaVal(trans.threshVal>0.05)  = lowAlpha;
