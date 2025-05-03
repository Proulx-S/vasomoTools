function roi = modifyRoi(roi,mod)


    
for r = 1:length(roi)
    for m = 1:length(mod)
        roi(r) = doIt(roi(r),mod{m},1);
    end
end



function roi = doIt(roi,mod,roiBaseInd)

    if contains(mod,'dilate')
        if length(mod)>6
            n = str2double(replace(mod(7:end),'p','.'));
        else
            n = 1;
        end
        mod = mod(1:6);
    end
    
    
    switch mod
        case 'dilate'
            if n==1.5
                % Dilate the mask by n pixel (including corners)
                roi.polyMask{end+1} = roi.polyMask{roiBaseInd};
                seOct = strel('octagon', 3);
                seDsk = strel('disk'   , 1);
                roi.polyMask{end} = imerode(imerode(imdilate(roi.polyMask{end}, seOct), seDsk), seDsk);
            else
                % Dilate the mask by n pixel
                seDsk = strel('disk', n);
                roi.polyMask{end+1} = imdilate(roi.polyMask{roiBaseInd}, seDsk);
            end
            roi.poly(end+1) = getMaskOutline(roi.polyMask{end},10);
            roi.poly(end).Vertices = roi.poly(end).Vertices - [1 1] + [roi.cropXlim(1) roi.cropYlim(1)];
            % figure('WindowStyle','docked');
            % imagesc(roi.cropMask); hold on
            % plot(roi.poly(end))
            roi.polyLabel{end+1} = [mod replace(num2str(n),'.','p')];

        case 'peakVox'
            roi.polyMask{end+1} = false(size(roi.polyMask{1}));
            [~,b] = max(roi.im.base.im(roi.im.base.mask));
            ind = false(size(roi.im.base.im(roi.im.base.mask))); ind(b) = true;
            roi.polyMask{end}(roi.im.base.mask) = ind;
            roi.poly(end+1) = getMaskOutline(roi.polyMask{end},10);
            roi.polyLabel{end+1} = mod;
        otherwise
            error('Unknown modification: %s',mod);
    end

    % figure('WindowStyle','docked');
    % tiledlayout(1,length(roi.polyMask)-1);
    % for i = 2:length(roi.polyMask)
    %     nexttile;
    %     imagesc(roi.cropXlim,roi.cropYlim,roi.polyMask{i}); axis image; xlim(roi.im.base.x+[-0.5 0.5]); ylim(roi.im.base.y+[-0.5 0.5]);
    %     hold on
    %     hP = plot(roi.poly(i));
    %     hP.FaceColor = 'none';
    %     hP = plot(roi.poly(1));
    %     hP.FaceColor = 'none';
    % end



