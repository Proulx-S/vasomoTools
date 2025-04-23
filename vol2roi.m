function roi = vol2roi(vol,roi)
    
    inFields  = fields(vol);
    outFields = {'psd' 'psdTrialGramMD'};

    dataMask = MRIread(vol.fMask); dataMask = dataMask.vol~=0;

    for d = 1:length(outFields)
        if ~ismember(outFields{d},inFields); continue; end
        switch outFields{d}
            case 'psd'
                for r = 1:length(roi)
                    roi(r).vec.mt.psd.vec   = vol.psd.PSD(:,:,:,:,:,roi(r).mask(dataMask),:,:);
                    % roi(r).vec.mt.psd.vecAv = mean(roi(r).vec.mt.psd.vec   ,6);
                    % roi(r).vec.mt.psd.vecEr = std( roi(r).vec.mt.psd.vec,[],6);
                    roi(r).vec.mt.psd.f     = vol.psd.f;
                    roi(r).vec.mt.psd.info  = vol.psd.info;
                    roi(r).vec.mt.psd.param = vol.psd.param;
                    roi(r).vec.mt.psd.K     = vol.psd.K;
                end            
            case 'psdTrialGramMD'
                for r = 1:length(roi)
                    roi(r).vec.mt.psdTrialGram.vec   = vol.psdTrialGramMD.vec.psdPC(:,:,:,:,:,roi(r).mask(dataMask),:,:);
                    % roi(r).vec.mt.psdTrialGram.vecAv = mean(roi(r).vec.mt.psdTrialGram.vec   ,6);
                    % roi(r).vec.mt.psdTrialGram.vecEr = std( roi(r).vec.mt.psdTrialGram.vec,[],6);
                    roi(r).vec.mt.psdTrialGram.f         = vol.psdTrialGramMD.f;
                    roi(r).vec.mt.psdTrialGram.t         = vol.psdTrialGramMD.t;
                    roi(r).vec.mt.psdTrialGram.info      = vol.psdTrialGramMD.info;
                    roi(r).vec.mt.psdTrialGram.param     = vol.psdTrialGramMD.param;
                    roi(r).vec.mt.psdTrialGram.K         = vol.psdTrialGramMD.K;
                    roi(r).vec.mt.psdTrialGram.onsetList = vol.psdTrialGramMD.onsetList;
                    roi(r).vec.mt.psdTrialGram.ondurList = vol.psdTrialGramMD.ondurList;
                end
            otherwise
                error('Unknown field: %s',outFields{d});
        end         
    end