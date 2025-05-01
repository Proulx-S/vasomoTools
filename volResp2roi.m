function roi = volResp2roi(vol,roi)

    
    inFields  = fields(vol.stats);
    outFields = {'fResp' 'fCondF' 'fCondF_qVal' 'fCondF_pVal' 'fCondCoef_adj'};
    dataMask = MRIread(char(vol.fMask)); dataMask = dataMask.vol~=0;

    for d = 1:length(outFields)
        if ~ismember(outFields{d},inFields) || isempty(vol.stats.(outFields{d})); continue; end
        switch outFields{d}
            case 'fResp'
                for r = 1:length(roi)
                    im = MRIread(char(vol.stats.fResp)); im = permute(im.vol,[4 1 2 3]);
                    roi(r).vec.resp.vec   = im(:,roi(r).mask(dataMask));
                    roi(r).vec.resp.t     = (0:size(im,1)-1).*vol.param.trDecon;
                    % roi(r).vec.resp.info
                    roi(r).vec.resp.param = vol.param;
                end            
            case 'fCondF'
                for r = 1:length(roi)
                    fName = char(vol.stats.fCondF);
                    im    = MRIread(fName);
                    im    = permute(im.vol,[4 1 2 3]);
                    %im
                    roi(r).im.F       = roi(r).im.base;
                    roi(r).im.F.fName = fName;
                    roi(r).im.F.im(:) = im(:,roi(r).mask);
                    %vec
                    roi(r).vec.resp.F   = im(:,roi(r).mask(dataMask));
                    roi(r).vec.resp.param = vol.param;
                    figure('WindowStyle','docked');
                    imagesc(dataMask)
                end
            case 'fCondF_qVal'
                for r = 1:length(roi)
                    im = MRIread(char(vol.stats.fCondF_qVal)); im = permute(im.vol,[4 1 2 3]);
                    roi(r).vec.resp.Q   = im(:,roi(r).mask(dataMask));
                    roi(r).vec.resp.param = vol.param;
                end
            case 'fCondF_pVal'
                for r = 1:length(roi)
                    im = MRIread(char(vol.stats.fCondF_pVal)); im = permute(im.vol,[4 1 2 3]);
                    roi(r).vec.resp.p   = im(:,roi(r).mask(dataMask));
                    roi(r).vec.resp.param = vol.param;
                end
            case 'fCondCoef_adj'
                for r = 1:length(roi)
                    im = MRIread(char(vol.stats.fCondCoef_adj)); im = permute(im.vol,[4 1 2 3]);
                    roi(r).vec.resp.coef   = im(:,roi(r).mask(dataMask));
                    roi(r).vec.resp.param = vol.param;
                end
            otherwise
                error('Unknown field: %s',outFields{d});
        end         
    end