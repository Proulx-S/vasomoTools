function roi = volPsd2roi(vol,roi)
    
    inFields  = fields(vol);
    outFields = {'psd' 'psdTrialGramMD'};

    if ~islogical(vol.fMask)
        dataMask = MRIread(vol.fMask); dataMask = dataMask.vol~=0;
    else
        dataMask = vol.fMask;
    end

    for d = 1:length(outFields)
        if ~ismember(outFields{d},inFields); continue; end
        switch outFields{d}
            case 'psd'
                for r = 1:length(roi)
                    roi(r).mt.psd.vec   = vol.psd.PSD(:,:,:,:,:,roi(r).cropMask(dataMask),:,:);
                    % roi(r).vec.mt.psd.vecAv = mean(roi(r).vec.mt.psd.vec   ,6);
                    % roi(r).vec.mt.psd.vecEr = std( roi(r).vec.mt.psd.vec,[],6);
                    roi(r).mt.psd.f     = vol.psd.f;
                    roi(r).mt.psd.info  = vol.psd.info;
                    roi(r).mt.psd.param = vol.psd.param;
                    roi(r).mt.psd.K     = vol.psd.K;
                end            
            case 'psdTrialGramMD'
                for r = 1:length(roi)
                    roi(r).mt.psdTrialGram.vec   = vol.psdTrialGramMD.vec.psdPC(:,:,:,:,:,roi(r).cropMask(dataMask),:,:);
                    % roi(r).vec.mt.psdTrialGram.vecAv = mean(roi(r).vec.mt.psdTrialGram.vec   ,6);
                    % roi(r).vec.mt.psdTrialGram.vecEr = std( roi(r).vec.mt.psdTrialGram.vec,[],6);
                    roi(r).mt.psdTrialGram.f         = vol.psdTrialGramMD.f;
                    roi(r).mt.psdTrialGram.t         = vol.psdTrialGramMD.t;
                    roi(r).mt.psdTrialGram.info      = vol.psdTrialGramMD.info;
                    roi(r).mt.psdTrialGram.param     = vol.psdTrialGramMD.param;
                    roi(r).mt.psdTrialGram.K         = vol.psdTrialGramMD.K;
                    roi(r).mt.psdTrialGram.onsetList = vol.psdTrialGramMD.onsetList;
                    roi(r).mt.psdTrialGram.ondurList = vol.psdTrialGramMD.ondurList;
                end
            otherwise
                error('Unknown field: %s',outFields{d});
        end         
    end