function threshOL(hIO,label,lowAlpha)
    if ~exist('label','var');       label = []; end
    if isempty(label);              label = 'actQ_crop'; end % thresholdMetric_fdrRegion: {actQ_crop, actQ_dilate2}
    if ~exist('lowAlpha','var'); lowAlpha = []  ; end
    if isempty(lowAlpha);        lowAlpha = 0; end

    if ismember(label,{'off' 'on'})
        set(hIO,'Visible',label);
        return;
    end
        
    for i = 1:length(hIO)
        doIt(hIO(i),label,lowAlpha);
        hAO(i) = hIO(i).Parent;
    end
    set(hAO,'Color','none');


    function doIt(hIO,label,lowAlpha)
        
        % Get the transparency
        trans = getTrans(hIO.UserData.roi,label,lowAlpha);
        
        % Apply transparency
        hIO.AlphaData = trans.alphaVal;
        hIO.UserData.trans = trans;



    function trans = getTrans(roi,label,lowAlpha)
        trans.label = label;
        trans.lowAlpha  = lowAlpha;
        switch label
            case 'actQ_crop'
                im2vec = true(size(roi.im.actP.im));
                vec = roi.im.actP.im(im2vec);
            case 'actQ_original'
                im2vec = roi.polyMask{ismember(roi.polyLabel,'original')};
                vec = roi.im.actP.im(im2vec);
            case 'actQ_dilate1'
                im2vec = roi.polyMask{ismember(roi.polyLabel,'dilate1')};
                vec = roi.im.actP.im(im2vec);
            case 'actQ_dilate1p5'
                im2vec = roi.polyMask{ismember(roi.polyLabel,'dilate1p5')};
                vec = roi.im.actP.im(im2vec);
            case 'actQ_dilate2'
                im2vec = roi.polyMask{ismember(roi.polyLabel,'dilate2')};
                vec = roi.im.actP.im(im2vec);
            case 'respQ_crop'
                im2vec = true(size(roi.im.respP.im));
                vec = roi.im.respP.im(im2vec);
            case 'respQ_original'
                im2vec = roi.polyMask{ismember(roi.polyLabel,'original')};
                vec = roi.im.respP.im(im2vec);
            case 'respQ_dilate1'
                im2vec = roi.polyMask{ismember(roi.polyLabel,'dilate1')};
                vec = roi.im.respP.im(im2vec);
            case 'respQ_dilate1p5'
                im2vec = roi.polyMask{ismember(roi.polyLabel,'dilate1p5')};
                vec = roi.im.respP.im(im2vec);
            case 'respQ_dilate2'
                im2vec = roi.polyMask{ismember(roi.polyLabel,'dilate2')};
                vec = roi.im.respP.im(im2vec);
        end
        trans.threshVal = zeros(size(im2vec));
        trans.alphaVal  = zeros(size(im2vec));
        trans.alphaVal(im2vec) = 1;
        trans.threshVal(im2vec) = mafdr(vec,'BHFDR',true);
        trans.alphaVal(trans.threshVal<=0.05 & im2vec) = 1;
        trans.alphaVal(trans.threshVal>0.05  & im2vec)  = lowAlpha;
