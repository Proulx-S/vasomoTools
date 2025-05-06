function roi = volPsd2roi(vol,roi)
    
    inFields  = fields(vol);
    outFields = {'psd' 'psdTrialGramMD'};

    
    if ~islogical(vol(1).fMask)
        dataMask = MRIread(vol.fMask); dataMask = dataMask.vol~=0;
    else
        tmp = diff(cat(4,vol.fMask),[],4);
        if any(tmp(:)); dbstack; error('masks differes'); end
        dataMask = vol(1).fMask;
    end

    for d = 1:length(outFields)
        if ~ismember(outFields{d},inFields); continue; end
        switch outFields{d}
            case 'psd'
                for r = 1:length(roi)
                    for R = 1:size(vol,1)
                        roi(r).mt.psd.vec(:,:,R,:,:,:,:,:)   = vol(R).psd.PSD(:,:,:,:,:,roi(r).cropMask(dataMask),:,:);
                        % roi(r).vec.mt.psd.vecAv = mean(roi(r).vec.mt.psd.vec   ,6);
                        % roi(r).vec.mt.psd.vecEr = std( roi(r).vec.mt.psd.vec,[],6);
                        roi(r).mt.psd.f(:,:,R,:,:,:,:,:)     = vol(R).psd.f;
                        roi(r).mt.psd.info{R}  = vol(R).psd.info;
                        roi(r).mt.psd.param{R} = vol(R).psd.param;
                        roi(r).mt.psd.K{R}     = vol(R).psd.K;
                    end

                    % average over runs
                    roi(r).mt.psd.vec = mean(roi(r).mt.psd.vec,3);
                    if max(max(abs(diff(roi(r).mt.psd.f,[],3))))./max(roi(r).mt.psd.f(end,:)) > 1e-5; dbstack; error('freqs are not equal'); end
                    roi(r).mt.psd.f     = mean(roi(r).mt.psd.f,3);
                    roi(r).mt.psd.info  = roi(r).mt.psd.info{1};
                    roi(r).mt.psd.param = roi(r).mt.psd.param{1};
                    if any(diff([roi(r).mt.psd.K{:}])); dbstack; error('Ks are not equal'); end
                    roi(r).mt.psd.K = roi(r).mt.psd.K{1};
                end            
            case 'psdTrialGramMD'
                for r = 1:length(roi)
                    for R = 1:size(vol,1)
                        roi(r).mt.psdTrialGram.vec(:,:,R,:,:,:,:,:)   = vol(R).psdTrialGramMD.vec.psdPC(:,:,:,:,:,roi(r).cropMask(dataMask),:,:);
                    % roi(r).vec.mt.psdTrialGram.vecAv = mean(roi(r).vec.mt.psdTrialGram.vec   ,6);
                    % roi(r).vec.mt.psdTrialGram.vecEr = std( roi(r).vec.mt.psdTrialGram.vec,[],6);
                        roi(r).mt.psdTrialGram.f(:,:,R,:,:,:,:,:)         = vol(R).psdTrialGramMD.f;
                        roi(r).mt.psdTrialGram.t(:,:,R,:,:,:,:,:)         = vol(R).psdTrialGramMD.t;
                        roi(r).mt.psdTrialGram.info{R}      = vol(R).psdTrialGramMD.info;
                        roi(r).mt.psdTrialGram.param{R}     = vol(R).psdTrialGramMD.param;
                        roi(r).mt.psdTrialGram.K{R}         = vol(R).psdTrialGramMD.K;
                        roi(r).mt.psdTrialGram.onsetList{R} = vol(R).psdTrialGramMD.onsetList;
                        roi(r).mt.psdTrialGram.ondurList{R} = vol(R).psdTrialGramMD.ondurList;
                    end

                    % average over runs
                    roi(r).mt.psdTrialGram.vec = mean(roi(r).mt.psdTrialGram.vec,3);
                    if max(max(abs(diff(roi(r).mt.psdTrialGram.f,[],3))))./max(roi(r).mt.psdTrialGram.f(end,:,:,:,:,:,:,end)) > 1e-5; dbstack; error('freqs are not equal'); end
                    if max(max(max(max(abs(diff(roi(r).mt.psdTrialGram.t,[],3))))))./max(roi(r).mt.psdTrialGram.t(end,end,:,:,:,:,:,end)) > 1e-5; dbstack; error('freqs are not equal'); end
                    roi(r).mt.psdTrialGram.f     = mean(roi(r).mt.psdTrialGram.f,3);
                    roi(r).mt.psdTrialGram.t     = mean(roi(r).mt.psdTrialGram.t,3);
                    roi(r).mt.psdTrialGram.info  = roi(r).mt.psdTrialGram.info{1};
                    roi(r).mt.psdTrialGram.param = roi(r).mt.psdTrialGram.param{1};
                    if any(diff([roi(r).mt.psdTrialGram.K{:}])); dbstack; error('Ks are not equal'); end
                    roi(r).mt.psdTrialGram.K = roi(r).mt.psdTrialGram.K{1};
                    roi(r).mt.psdTrialGram.onsetList = roi(r).mt.psdTrialGram.onsetList{1};
                    roi(r).mt.psdTrialGram.ondurList = roi(r).mt.psdTrialGram.ondurList{1};
                end
            otherwise
                error('Unknown field: %s',outFields{d});
        end         
    end