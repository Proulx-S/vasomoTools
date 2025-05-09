function roi = volPsd2roi(vol,roi,outFields,skipSpSV)
    if ~exist('outFields','var'); outFields = []   ; end
    if ~exist('skipSpSV','var');   skipSpSV = false; end
    if isempty(outFields)
        outFields = {'psd' 'psdTrialGramMD'};
    end
    inFields  = fields(vol);
    
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
            case 'svd'
                for r = 1:length(roi)
                    vec   = cell(size(vol));
                    if ~skipSpSV
                        vecSV = cell(size(vol));
                    end
                    f     = cell(size(vol));
                    info  = cell(size(vol));
                    param = cell(size(vol));
                    K     = cell(size(vol));
                    % figure('WindowStyle','docked');
                    for R = 1:size(vol,1)
                        vec{R}   = vol(R).svd.COH(:,:,:,:,:,:,:,1);
                        if ~skipSpSV
                            vecSV{R} = vol(R).svd.spSV(:,:,:,:,:,:,:,1);
                        end
                        % roi(r).vec.mt.psd.vecAv = mean(roi(r).vec.mt.psd.vec   ,6);
                        % roi(r).vec.mt.psd.vecEr = std( roi(r).vec.mt.psd.vec,[],6);
                        f{R}     = vol(R).svd.f;
                        info{R}  = vol(R).svd.info;
                        param{R} = vol(R).svd.param;
                        K{R}     = vol(R).svd.K;
                    end

                    % average over runs
                    vec   = cat(3,vec{:});
                    if ~skipSpSV
                        vecSV = cat(3,vecSV{:});
                    end
                    f   = cat(3,f{:});
                    if max(max(abs(diff(f,[],3))))./max(f(end,:)) > 1e-5; dbstack; error('freqs are not equal'); end
                    roi(r).mt.svd.vec   = mean(vec,3);
                    if ~skipSpSV
                        roi(r).mt.svd.vecSV = mean(vecSV,3);
                    end
                    roi(r).mt.svd.f     = mean(f,3);
                    roi(r).mt.svd.info  = info{1};
                    roi(r).mt.svd.param = param{1};
                    if any(diff([K{:}])); dbstack; error('Ks are not equal'); end
                    roi(r).mt.svd.K = K{1};
                end  
            case 'psd'
                for r = 1:length(roi)
                    vec   = cell(size(vol));
                    f     = cell(size(vol));
                    info  = cell(size(vol));
                    param = cell(size(vol));
                    K     = cell(size(vol));
                    % figure('WindowStyle','docked');
                    for R = 1:size(vol,1)
                        vec{R}   = vol(R).psd.PSD(:,:,:,:,:,roi(r).cropMask(dataMask),:,:);
                        % roi(r).vec.mt.psd.vecAv = mean(roi(r).vec.mt.psd.vec   ,6);
                        % roi(r).vec.mt.psd.vecEr = std( roi(r).vec.mt.psd.vec,[],6);
                        f{R}     = vol(R).psd.f;
                        info{R}  = vol(R).psd.info;
                        param{R} = vol(R).psd.param;
                        K{R}     = vol(R).psd.K;
                    end

                    % average over runs
                    vec = cat(3,vec{:});
                    f   = cat(3,f{:});
                    if max(max(abs(diff(f,[],3))))./max(f(end,:)) > 1e-5; dbstack; error('freqs are not equal'); end
                    roi(r).mt.psd.vec   = mean(vec,3);
                    roi(r).mt.psd.f     = mean(f,3);
                    roi(r).mt.psd.info  = info{1};
                    roi(r).mt.psd.param = param{1};
                    if any(diff([K{:}])); dbstack; error('Ks are not equal'); end
                    roi(r).mt.psd.K = K{1};
                end
            case 'svdTrialGramMD'
                for r = 1:length(roi)
                    vec   = cell(size(vol));
                    % vecSV = cell(size(vol));
                    f     = cell(size(vol));
                    t     = cell(size(vol));
                    info  = cell(size(vol));
                    param = cell(size(vol));
                    K     = cell(size(vol));
                    onsetList = cell(size(vol));
                    ondurList = cell(size(vol));
                    for R = 1:size(vol,1)
                        vec{R}       = vol(R).svdTrialGramMD.vec.cohEPC(:,:,:,:,:,:,:,1);
                        % vecSV{R}     = vol(R).svdTrialGramMD.vec.spSV(:,:,:,:,:,:,:,1);
                        f{R}         = vol(R).svdTrialGramMD.f;
                        t{R}         = vol(R).svdTrialGramMD.t;
                        info{R}      = vol(R).svdTrialGramMD.info;
                        param{R}     = vol(R).svdTrialGramMD.param;
                        K{R}         = vol(R).svdTrialGramMD.K;
                        onsetList{R} = vol(R).svdTrialGramMD.onsetList;
                        ondurList{R} = vol(R).svdTrialGramMD.ondurList;
                    end

                    % average over runs
                    vec   = cat(3,vec{:});
                    % vecSV = cat(3,vecSV{:});
                    f   = cat(3,f{:});
                    t   = cat(3,t{:});
                    if max(max(abs(diff(f,[],3))))./max(f(end,:,:,:,:,:,:,end)) > 1e-5; dbstack; error('freqs are not equal'); end
                    if max(max(max(max(abs(diff(t,[],3))))))./max(t(end,:,:,:,:,:,:,end)) > 1e-5; dbstack; error('freqs are not equal'); end
                    if any(diff([K{:}])); dbstack; error('Ks are not equal'); end
                    roi(r).mt.svdTrialGram.vec       = mean(vec,3);
                    % roi(r).mt.svdTrialGram.vecSV     = mean(vecSV,3);
                    roi(r).mt.svdTrialGram.f         = mean(f,3);
                    roi(r).mt.svdTrialGram.t         = mean(t,3);
                    roi(r).mt.svdTrialGram.info      = info{1};
                    roi(r).mt.svdTrialGram.param     = param{1};
                    roi(r).mt.svdTrialGram.K         = K{1};
                    roi(r).mt.svdTrialGram.onsetList = onsetList{1};
                    roi(r).mt.svdTrialGram.ondurList = ondurList{1};
                end
            case 'psdTrialGramMD'
                for r = 1:length(roi)
                    vec   = cell(size(vol));
                    f     = cell(size(vol));
                    t     = cell(size(vol));
                    info  = cell(size(vol));
                    param = cell(size(vol));
                    K     = cell(size(vol));
                    onsetList = cell(size(vol));
                    ondurList = cell(size(vol));
                    for R = 1:size(vol,1)
                        vec{R}       = vol(R).psdTrialGramMD.vec.psdPC(:,:,:,:,:,roi(r).cropMask(dataMask),:,:);
                        f{R}         = vol(R).psdTrialGramMD.f;
                        t{R}         = vol(R).psdTrialGramMD.t;
                        info{R}      = vol(R).psdTrialGramMD.info;
                        param{R}     = vol(R).psdTrialGramMD.param;
                        K{R}         = vol(R).psdTrialGramMD.K;
                        onsetList{R} = vol(R).psdTrialGramMD.onsetList;
                        ondurList{R} = vol(R).psdTrialGramMD.ondurList;
                    end

                    % average over runs
                    vec = cat(3,vec{:});
                    f   = cat(3,f{:});
                    t   = cat(3,t{:});
                    if max(max(abs(diff(f,[],3))))./max(f(end,:,:,:,:,:,:,end)) > 1e-5; dbstack; error('freqs are not equal'); end
                    if max(max(max(max(abs(diff(t,[],3))))))./max(t(end,:,:,:,:,:,:,end)) > 1e-5; dbstack; error('freqs are not equal'); end
                    if any(diff([K{:}])); dbstack; error('Ks are not equal'); end
                    roi(r).mt.psdTrialGram.vec       = mean(vec,3);
                    roi(r).mt.psdTrialGram.f         = mean(f,3);
                    roi(r).mt.psdTrialGram.t         = mean(t,3);
                    roi(r).mt.psdTrialGram.info      = info{1};
                    roi(r).mt.psdTrialGram.param     = param{1};
                    roi(r).mt.psdTrialGram.K         = K{1};
                    roi(r).mt.psdTrialGram.onsetList = onsetList{1};
                    roi(r).mt.psdTrialGram.ondurList = ondurList{1};
                end
            otherwise
                error('Unknown field: %s',outFields{d});
        end         
    end