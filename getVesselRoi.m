function roi = getVesselRoi(mask,imField,im,cropSz)
if ~iscell(imField); imField = {imField}; end
if ~iscell(im); im = {im}; end

if all(mask(:)==0); roi = []; return; end

Pall = getMaskOutline(mask,10);
P = regions(Pall);
roi = repmat(struct,size(P) + [1 0]);
for p = 1:length(P)+1
    if p~=length(P)+1
        roi(p).poly  = P(p);
        roi(p).mask  = poly2mask(roi(p).poly.Vertices(:,1),roi(p).poly.Vertices(:,2),size(mask,1),size(mask,2));
        roi(p).label = ['vessel' num2str(p,'%02i')];
    else
        roi(p).poly  = Pall;
        roi(p).mask  = mask;
        roi(p).label = 'vesselAll';
    end

    %center of mass
    [rows, cols] = ndgrid(1:size(roi(p).mask, 1), 1:size(roi(p).mask, 2));
    baseInd = ismember(imField,'base');
    if any(baseInd)
        %of the masked baseline image
        imCom = im{baseInd};
    else
        %of the mask
        imCom = roi(p).mask;
    end
    roi(p).com(1) = sum( cols(roi(p).mask) .* imCom(roi(p).mask) ) / sum(imCom(roi(p).mask));
    roi(p).com(2) = sum( rows(roi(p).mask) .* imCom(roi(p).mask) ) / sum(imCom(roi(p).mask));

    %crop all images to all rois
    switch roi(p).label            
        case 'vesselAll'
            x = [find(any(roi(p).mask,1),1,'first') find(any(roi(p).mask,1),1,'last')];
            y = [find(any(roi(p).mask,2),1,'first') find(any(roi(p).mask,2),1,'last')];
            x = mean(x) + [-1 1].*max([diff(x) diff(y)])/2;
            y = mean(y) + [-1 1].*max([diff(x) diff(y)])/2;
            for i = 1:length(imField)
                roi(p).im.(imField{i}).x  = round(x + [-1 1].*cropSz/2);
                roi(p).im.(imField{i}).y  = round(y + [-1 1].*cropSz/2);
                % roi(p).im.(imField{i}).x  = round([find(any(roi(p).mask,1),1,'first') find(any(roi(p).mask,1),1,'last')] + [-1 1].*cropSz/2);
                % roi(p).im.(imField{i}).y  = round([find(any(roi(p).mask,2),1,'first') find(any(roi(p).mask,2),1,'last')] + [-1 1].*cropSz/2);
            end
        otherwise
            x = round( roi(p).com(1) + [-1 1].*cropSz/2 );
            y = round( roi(p).com(2) + [-1 1].*cropSz/2 );
            tmp = roi(p).mask; tmp(y(1):y(2),x(1):x(2)) = false;
            if any(tmp(:))
                warning(['cropped image does not include all of ' roi(p).label newline 'consider increasing cropSz'])
            end
            for i = 1:length(imField)
                roi(p).im.(imField{i}).x = x;
                roi(p).im.(imField{i}).y = y;
                % roi(p).im.(imField{i}).im = im{i}(roi(p).im.(imField{i}).y(1):roi(p).im.(imField{i}).y(2),roi(p).im.(imField{i}).x(1):roi(p).im.(imField{i}).x(2),:,:,:);
            end
    end
    for i = 1:length(imField)
        roi(p).im.(imField{i}).im = im{i}(roi(p).im.(imField{i}).y(1):roi(p).im.(imField{i}).y(2),roi(p).im.(imField{i}).x(1):roi(p).im.(imField{i}).x(2),:,:,:);
    end
    roi(p).ts = [];
    roi(p).fs = [];
end

