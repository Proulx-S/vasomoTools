function roi = getVesselRoi2(label,imField,im,cropSz)

    %% Massage input
    if ~iscell(imField); imField = {imField}; end
    if ~iscell(im); im = {im}; end
    % % read mask when specified as a filename
    % if ischar(mask) || iscell(mask); mask = MRIread(char(mask)); end
    % if isstruct(mask); mask = mask.vol; end
    % if nnz(size(mask)>1)==2; mask = squeeze(mask); else error('Your supposed to have a single slice here'); end
    % % exit if empty
    % if all(mask(:)==0); roi = []; return; end


    %% Convert label to mask
    if isstruct(label)
        mriLabel = MRIread(label.f);
        mask = repmat(struct,1,length(label.label));
        for l = 1:length(label.label)
            mask(l).vol      = mriLabel.vol == label.labelVal(l);
            mask(l).label    = label.label{l};
            mask(l).labelVal = label.labelVal(l);
        end
        % merge masks
        mask(end+1).vol    = any(cat(4,mask(ismember({mask.label},{'Left-vessel' 'Right-vessel'})).vol),4);
        mask(end).label    = 'unknown';
        mask(end).labelVal = 30;
        mask(ismember({mask.label},{'Left-vessel' 'Right-vessel'})) = [];
    else
        dbstack; error('label should be a struct');
    end

    %% Read in data to crop
    fIm = cell(size(im));
    for d = 1:length(im)
        if ischar(im{d}) && ~isempty(im{d})
            fIm{d} = im{d};
                im{d} = MRIread(fIm{d});
                im{d} = im{d}.vol;
        elseif isempty(im{d})
            im{d} = [];
        else
            dbstack; error('im should be a char or a struct');
        end
    end

    %% Get individual vessel ROIs
    roi = cell(1,length(mask));
    for l = 1:length(mask)
        roi{l} = doIt(mask(l).vol,mask(l).label,imField,im,fIm,cropSz);
    end
    roi = cat(1,roi{:});
    roi(cellfun('isempty',{roi.label})) = [];





function roi = doIt(mask,label,imField,im,fIm,cropSz)

if ~any(mask(:))
    roi.class     = [];
    roi.id        = [];
    roi.label     = [];
    roi.cropMask  = [];
    roi.cropSz    = [];
    roi.cropXlim  = [];
    roi.cropYlim  = [];
    roi.com       = [];
    roi.poly      = [];
    roi.polyMask  = {};
    roi.polyLabel = {};
    roi.im        = [];
    return;
end


Pall = getMaskOutline(mask,10);
P = regions(Pall);

% roi = repmat(struct,size(P) + [1 0]);
roi = repmat(struct,size(P));
for p = 1:length(P)
    roi(p).class = label;
    roi(p).id    = p;
    roi(p).label = [label num2str(p,'%02i')];
    
    %center of mass
    maskRoi  = poly2mask(P(p).Vertices(:,1),P(p).Vertices(:,2),size(mask,1),size(mask,2));
    [rows, cols] = ndgrid(1:size(maskRoi, 1), 1:size(maskRoi, 2));
    baseInd = ismember(imField,'base');
    if any(baseInd)
        %of the masked baseline image
        imCom = im{baseInd};
    else
        %of the mask
        imCom = maskRoi;
    end
    com(1) = sum( cols(maskRoi) .* imCom(maskRoi) ) / sum(imCom(maskRoi));
    com(2) = sum( rows(maskRoi) .* imCom(maskRoi) ) / sum(imCom(maskRoi));
    x = round( com(1) + [-1 1].*cropSz/2 );
    y = round( com(2) + [-1 1].*cropSz/2 );

    % get cropping mask
    cropMask = false(size(maskRoi));
    cropMask(y(1):y(2),x(1):x(2)) = true;
    
    %check if the cropped image includes all of the roi
    if any(maskRoi(~cropMask))
        warning(['cropped image does not include all of ' roi(p).label newline 'consider increasing cropSz'])
    end

    %store some info
    roi(p).cropMask = cropMask;
    roi(p).cropSz   = cropSz;
    roi(p).cropXlim = x;
    roi(p).cropYlim = y;
    roi(p).com      = com;

    roi(p).poly = P(p);
    roi(p).polyMask    = {false(cropSz+[1 1])};
    roi(p).polyMask{1}(:) = maskRoi(cropMask);
    roi(p).polyLabel = {'original'};
    for i = 1:length(imField)
        roi(p).im.(imField{i}).fName = fIm{i};
        roi(p).im.(imField{i}).x     = x;
        roi(p).im.(imField{i}).y     = y;
        roi(p).im.(imField{i}).mask  = roi(p).polyMask{1};
        roi(p).im.(imField{i}).im = [];
        if ~isempty(im{i})
            roi(p).im.(imField{i}).im = im{i}(roi(p).im.(imField{i}).y(1):roi(p).im.(imField{i}).y(2),roi(p).im.(imField{i}).x(1):roi(p).im.(imField{i}).x(2),:,:,:);
        end
    end    
end

