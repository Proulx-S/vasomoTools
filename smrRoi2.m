function [roi,hF,hA,rCond] = smrRoi2(rCond,metric,roi,H)
    if ~exist('roi','var');           roi = struct; end    
    if ~exist('H','var');               H = []    ; end
    if ~exist('metric','var');     metric = {}; end
    if isempty(metric);            metric = {'resp'}; end

    redoMT = 0;

    metric = cellstr(metric);

    %% Assert
    if iscell(roi)
        roi = [roi{:}];
    end
    if iscell(H)
        H = [H{:}];
    end
    roi = roi(:);
    if ~isempty(H)
        H   = H(:);
        % if all(size(roi) == flip(size(H))); roi = roi'; end
        if length(roi) ~= length(H); dbstack; error('roi and H must have the same dimensions'); end

        %% Setup figure
        for m = 1:length(metric)
            hF{m} = figure('WindowStyle','docked');
            if isempty(roi)
            else
                hA{m} = cell(size(H));
                for i = 1:length(H)
                    hA{m}{i} = axes(hF{m},'Position',H(i).Position,'Box','on');
                end
                hA{m} = [hA{m}{:}]; hA{m} = hA{m}(:);
            end
        end
    end

    m = ismember(metric,{'coh_dilate1' 'coh_original' 'cohTrialGram_dilate1' 'cohTrialGram_original'});
    %% Recompute mt (relevant for coherence within single vessels)
    if any(m)
        m = find(m,1);
        % if length(metric)>1; dbstack; error('multiple metrics not supported'); end
        
        polyLabel = strsplit(metric{m},'_'); polyLabel = polyLabel{2};
        

        K = [
            rCond.volMt.run(1).psdGram.K
            rCond.volMt.run(1).psdTrialGramMD.K
            rCond.volMt.run(1).psd.K
            ]; % K(end)->full timeseries, K(1)->time-resolved, K(2)->trial-triggered based on missing data
        % K   = [1 3 5]; % K(end)->full timeseries, K(1)->time-resolved, K(2)->trial-triggered based on missing data
        W   = [];
        win = rCond.volMt.run(1).param.psdTrialGram.dsgn.winSec; % in seconds [lenght, step]
        skipSVD = 0;
        skipPSD = 0;

        for r = 1:length(roi)

            disp('--------------------------------');
            disp(['Doing vessel ',num2str(r),' of ',num2str(length(roi))]);
            disp('--------------------------------');

            % create mask for polyRoi within vesselRoi
            mask = roi(r).cropMask;
            polyInd = ismember(roi(r).polyLabel,polyLabel);
            mask(mask) = roi(r).polyMask{polyInd};

            % recompute svd within the mask
            rCondX = runFullMT6(rCond,W,K,win,rCond.dsgn,mask,skipSVD,skipPSD,0);

            % put new results (svd only) into roi
            roi(r) = volPsd2roi(rCondX.volMt.run,roi(r),{'svd' 'svdTrialGramMD'});
        end



        % f   = squeeze(rCondX.volMt.run(3).psdTrialGramMD.f);
        % t   = squeeze(mean(rCondX.volMt.run(3).psdTrialGramMD.t - rCondX.volMt.run(3).psdTrialGramMD.onsetList',2));
        % psd = mean(rCondX.volMt.run(1).psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,:,1),6);
        % psd = cat(3,psd,mean(rCondX.volMt.run(2).psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,:,1),6));
        % psd = cat(3,psd,mean(rCondX.volMt.run(3).psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,:,1),6));
        % figure('WindowStyle','docked');
        % imagesc(mean(t,1),f,squeeze(mean(psd,3)));
        % set(gca,'ColorScale','log');


        % f   = squeeze(rCondX.volMt.run(3).svd.f);
        % coh = squeeze(rCondX.volMt.run(3).svd.COH(:,:,:,:,:,:,:,1));
        % figure('WindowStyle','docked');
        % plot(f,coh);
        
        % f   = squeeze(rCondX.volMt.run(3).svdTrialGramMD.f);
        % t   = squeeze(mean(rCondX.volMt.run(3).svdTrialGramMD.t - rCondX.volMt.run(3).svdTrialGramMD.onsetList',2));
        % winSz = mean(diff(t,[],1),7);
        % coh = squeeze(rCondX.volMt.run(3).svdTrialGramMD.vec.cohEPC(:,:,:,:,:,:,:,1));
        % coh = cat(3,coh,squeeze(rCondX.volMt.run(2).svdTrialGramMD.vec.cohEPC(:,:,:,:,:,:,:,1)))
        % coh = cat(3,coh,squeeze(rCondX.volMt.run(1).svdTrialGramMD.vec.cohEPC(:,:,:,:,:,:,:,1)));
        % figure('WindowStyle','docked');
        % imagesc(mean(t,1),f,mean(coh,3));
        

        % psd =           squeeze(rCondX.volMt.run(3).psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,:,1));
        % psd = cat(3,psd,squeeze(rCondX.volMt.run(2).psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,:,1)));
        % psd = cat(3,psd,squeeze(rCondX.volMt.run(1).psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,:,1)));
        % figure('WindowStyle','docked');
        % imagesc(mean(t,1),f,mean(psd,3));
        % set(gca,'ColorScale','log');
        % colorscale = 'log';
    end


    if any(ismember(metric,{'psd_dilate1_actQ' 'psdTrialGram_dilate1_actQ'}))
        if redoMT
            K   = [1 3 5]; % K(end)->full timeseries, K(1)->time-resolved, K(2)->trial-triggered based on missing data
            W   = [];
            win = [22 0.840]; % in seconds [lenght, step]
            skipSVD = 1;
            skipPSD = 0;
            mask = any(cat(4,roi.cropMask),4);
            if ~all(mask(:)==roi(1).mt.psd.param.fMask(:)); dbstack; error('mask does not match'); end
            rCond = runFullMT6(rCond,W,K,win,rCond.dsgn,mask,skipSVD,skipPSD,0,1);
            roi = volPsd2roi(rCond.volMt.runAv,roi);
            

            winList = squeeze(mean(roi(1).mt.psdTrialGram.t - roi(1).mt.psdTrialGram.onsetList',2));
            [~,b] = min(abs(winList(1,:)))
            roi(1).mt.psdTrialGram.postStimWin0 = b;
            disp('first post-stimulus window:')
            disp(winList(:,b))
            disp('last post-stimulus window:')
            disp(winList(:,end))
            disp('ISI:')
            disp(mean(diff(roi(1).mt.psdTrialGram.onsetList)))
        end
    end





    if isfield(roi(1),'im') && isfield(roi(1).im,'resp')
        roiDataFlag = true;
    else
        roiDataFlag = false;
    end
    

    %% Plot response time course
    if roiDataFlag
        lineStyle = {'-','-'};
        lineColor = {[0 0 0],[0.5 0.5 0.5]};
        hTs = cell(length(roi),length(metric));
        for i = 1:length(roi)
            for m = 1:length(metric)
                if ~isempty(H)
                    hold(hA{m}(i),'on');
                end
                
                polyLabel = strsplit(metric{m},'_');
                if strcmp(metric{m},'respPhs_peakBasePhsInDilate1')
                    polyLabel = 'dilate1';
                    indIm         = roi(i).polyMask{ismember(roi(i).polyLabel,polyLabel)};
                
                else
                    polyLabel = polyLabel{2};
                    indIm         = roi(i).polyMask{ismember(roi(i).polyLabel,polyLabel)};
                    indSig        = false(size(roi(i).im.actP.im));
                    indSig(indIm) = mafdr(roi(i).im.actP.im(indIm),'BHFDR',true)<0.05;
                    indNeg        = roi(i).im.act.im(:,:,:,1)<0;
                    indPos        = roi(i).im.act.im(:,:,:,1)>0;
                end
                nTrial        = [1; roi(i).nTrial; 1];
                
                
                switch metric{m}
                    case {'coh_dilate1' 'coh_original'}
                        % coherence within polyRoi within vessel crop ROI
                        roi(i).smr{m}.vec     = roi(i).mt.svd.vec;
                        roi(i).smr{m}.nVox    = [nnz(indIm&indNeg) nnz(indIm&indPos)];
                        roi(i).smr{m}.nVoxRoi = nnz(indIm);
                        roi(i).smr{m}.neg     = cat(6,true(1,1,1,1,1,roi(i).smr{m}.nVox(1)),false(1,1,1,1,1,roi(i).smr{m}.nVox(2)));
                        roi(i).smr{m}.f       = roi(i).mt.svd.f;
                        roi(i).smr{m}.label   = 'coh';
                        roi(i).smr{m}.metric  = metric{m};
                        roi(i).smr{m}.info    = roi(i).mt.svd.info;
                        roi(i).smr{m}.nRun    = roi(i).R;

                        if ~isempty(H)
                            hTs{i,m} = plot(hA{m}(i),squeeze(roi(i).smr{m}.f),squeeze(roi(i).smr{m}.vec),'k');
                            yLim{m}    = [1/K(3) 1];
                            yScale{m}    = 'linear';
                            yLabel{m}    = 'coherence';
                            cLim{m}      = [];
                            majorGrid{m} = 'on';
                            minorGrid{m} = 'on';
                        end

                    case {'cohTrialGram_dilate1' 'cohTrialGram_original'}
                        % polyLabel = strsplit(metric{m},'_'); polyLabel = polyLabel{2};
                        % indIm         = roi(i).polyMask{ismember(roi(i).polyLabel,polyLabel)};
                        % indSig        = false(size(roi(i).im.actP.im));
                        % indSig(indIm) = mafdr(roi(i).im.actP.im(indIm),'BHFDR',true)<0.05;
                        
                        % PSDtrialGram (time-frequency spectrogram) within the vessel ROI dilated by 1 voxel,
                        % including only active voxels based on SPMG2 activation detection
                        roi(i).smr{m}.vec     = roi(i).mt.svdTrialGram.vec;
                        roi(i).smr{m}.nVox    = [nnz(indIm&indNeg) nnz(indIm&indPos)];
                        roi(i).smr{m}.nVoxRoi = nnz(indIm);
                        roi(i).smr{m}.neg     = cat(6,true(1,1,1,1,1,roi(i).smr{m}.nVox(1)),false(1,1,1,1,1,roi(i).smr{m}.nVox(2)));
                        roi(i).smr{m}.f       = roi(i).mt.svdTrialGram.f;
                        roi(i).smr{m}.t       = mean(roi(i).mt.svdTrialGram.t - roi(i).mt.svdTrialGram.param.dsgn.onsetList,2);
                        roi(i).smr{m}.label   = 'cohTrialGram';
                        roi(i).smr{m}.metric  = metric{m};
                        roi(i).smr{m}.info    = roi(i).mt.svdTrialGram.info;
                        roi(i).smr{m}.nRun    = roi(i).R;

                        if ~isempty(H)
                            imagesc(hA{m}(i),squeeze(mean(roi(i).smr{m}.t,1)),squeeze(roi(i).smr{m}.f),squeeze(roi(i).smr{m}.vec));
                            hA{m}(i).ColorScale = 'linear';
                            yLim{m}      = [];
                            yScale{m}    = 'linear';
                            yLabel{m}    = 'Hz';
                            cLim{m}      = [1/K(2) 1];
                            majorGrid{m} = 'off';
                            minorGrid{m} = 'off';
                        end


                    case {'psd_dilate1_actQ' 'psd_original_actQ'}
                        % polyLabel = strsplit(metric{m},'_'); polyLabel = polyLabel{2};
                        % including only active voxels based on SPMG2 activation detection
                        % indIm         = roi(i).polyMask{ismember(roi(i).polyLabel,polyLabel)};
                        % indSig        = false(size(roi(i).im.actP.im));
                        % indSig(indIm) = mafdr(roi(i).im.actP.im(indIm),'BHFDR',true)<0.05;
                        
                        % PSD within the vessel ROI dilated by 1 voxel,
                        vecNeg = roi(i).mt.psd.vec(:,:,:,:,:,indIm&indSig&indNeg,:,:);
                        vecPos = roi(i).mt.psd.vec(:,:,:,:,:,indIm&indSig&indPos,:,:);
                        roi(i).smr{m}.vec     = cat(6,vecNeg,vecPos);
                        roi(i).smr{m}.nVox    = [nnz(indIm&indSig&indNeg) nnz(indIm&indSig&indPos)];
                        roi(i).smr{m}.nVoxRoi = nnz(indIm);
                        roi(i).smr{m}.neg     = cat(6,true(1,1,1,1,1,roi(i).smr{m}.nVox(1)),false(1,1,1,1,1,roi(i).smr{m}.nVox(2)));
                        roi(i).smr{m}.f       = roi(i).mt.psd.f;
                        roi(i).smr{m}.label   = 'psd';
                        roi(i).smr{m}.metric  = metric{m};
                        roi(i).smr{m}.info    = roi(i).mt.psd.info;
                        roi(i).smr{m}.nRun    = roi(i).R;

                        if ~isempty(H)
                            hTs{i,m} = plot(hA{m}(i),squeeze(roi(i).smr{m}.f),squeeze(mean(roi(i).smr{m}.vec,6)),'k');
                            yLim{m}      = [];
                            yScale{m}    = 'log';
                            yLabel{m}    = 'PSD';
                            cLim{m}      = [];
                            majorGrid{m} = 'on';
                            minorGrid{m} = 'on';
                        end

                    case {'psdTrialGram_dilate1_actQ' 'psdTrialGram_original_actQ'}
                        % polyLabel = strsplit(metric{m},'_'); polyLabel = polyLabel{2};
                        % indIm         = roi(i).polyMask{ismember(roi(i).polyLabel,polyLabel)};
                        % indSig        = false(size(roi(i).im.actP.im));
                        % indSig(indIm) = mafdr(roi(i).im.actP.im(indIm),'BHFDR',true)<0.05;
                        
                        % PSDtrialGram (time-frequency spectrogram) within the vessel ROI dilated by 1 voxel,
                        % including only active voxels based on SPMG2 activation detection
                        vecNeg = roi(i).mt.psdTrialGram.vec(:,:,:,:,:,indIm&indSig&indNeg,:,:);
                        vecPos = roi(i).mt.psdTrialGram.vec(:,:,:,:,:,indIm&indSig&indPos,:,:);
                        roi(i).smr{m}.vec     = cat(6,vecNeg,vecPos);
                        roi(i).smr{m}.nVox    = [nnz(indIm&indSig&indNeg) nnz(indIm&indSig&indPos)];
                        roi(i).smr{m}.nVoxRoi = nnz(indIm);
                        roi(i).smr{m}.neg     = cat(6,true(1,1,1,1,1,roi(i).smr{m}.nVox(1)),false(1,1,1,1,1,roi(i).smr{m}.nVox(2)));
                        roi(i).smr{m}.f       = roi(i).mt.psdTrialGram.f;
                        roi(i).smr{m}.t       = mean(roi(i).mt.psdTrialGram.t - roi(i).mt.psdTrialGram.param.dsgn.onsetList,2);
                        roi(i).smr{m}.label   = 'psdTrialGram';
                        roi(i).smr{m}.metric  = metric{m};
                        roi(i).smr{m}.info    = roi(i).mt.psdTrialGram.info;
                        roi(i).smr{m}.nRun    = roi(i).R;

                        if ~isempty(H)
                            imagesc(hA{m}(i),squeeze(mean(roi(i).smr{m}.t,1)),squeeze(roi(i).smr{m}.f),squeeze(mean(roi(i).smr{m}.vec,6)));
                            hA{m}(i).ColorScale = 'log';
                            yLim{m}      = [];
                            yScale{m}    = 'linear';
                            yLabel{m}    = 'Hz';
                            cLim{m}      = [];
                            majorGrid{m} = 'off';
                            minorGrid{m} = 'off';
                        end
                        
                    case {'respPhs_peakBasePhsInDilate1'}
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%% CONTINUTE HERE %%%%%%
                        im     = permute(roi(i).im.resp.im,[4 1 2 3]);
                        imBase = roi(i).im.basePhase.im;
                        bckgrndMask = getRoiBckgrndMask(roi(i).im.base.im,0);
                        [xBase,yBase] = pol2cart(imBase,1);
                        [bckgrndPhs,~] = cart2pol(mean(xBase(bckgrndMask)),mean(yBase(bckgrndMask)));
                        im = im - bckgrndPhs;
                        [~,ind] = max(abs(im(1,:)));
                        ts   = im(:,ind);
                        t     = ((0:size(ts,1)-1).*roi(i).im.resp.dt)';

                        roi(i).smr{m}.vec     = ts;
                        roi(i).smr{m}.vecEr   = [];
                        roi(i).smr{m}.nVox    = 1;
                        roi(i).smr{m}.nVoxRoi = [nnz(indIm)];
                        roi(i).smr{m}.neg     = [];
                        roi(i).smr{m}.t       = t;
                        roi(i).smr{m}.nTrial  = nTrial;
                        roi(i).smr{m}.label   = {''};
                        roi(i).smr{m}.metric  = metric{m};
                        roi(i).smr{m}.info    = strjoin({'time' '???' 'vox'},' x ');
                        roi(i).smr{m}.bckgrndPhs = bckgrndPhs;

                        if ~isempty(H)
                            % if ~isempty(roi(i).smr{m}.vecEr)
                                t     = roi(i).smr{m}.t;
                                vecAv =      roi(i).smr{m}.vec;
                                % thi is a very conservative approach because squared deviations cannot average out
                                % (for accurate measure of error of the voxel-averaged responses, one needs to work from the residuals or perform another fit on the voxel-averaged full timeseries)
                                
                                hTs{i,m}(1) = plot(hA{m}(i),t,vecAv);
                                
                                ind   = ~roi(i).smr{m}.neg;
                                t     = roi(i).smr{m}.t;
                                vecAv =      mean( roi(i).smr{m}.vec(  :,:,ind)    ,3);
                                vecEr = sqrt(mean( roi(i).smr{m}.vecEr(:,:,ind).^2 ,3));
                                
                                % hTs{i,m}(2) = shplot2(t,vecAv,vecEr,hA{m}(i));
                            % else
                            %     dbstack; error('double check that');
                            %     hTs{i,m} = plot(hA(i),roi(i).ts{m}.t,roi(i).ts{m}.vec);
                            % end
                            yLim{m}{i}      = [-pi pi] - roi(i).smr{m}.bckgrndPhs;
                            yScale{m}    = 'linear';
                            yLabel{m}    = 'MR pc change rel. baseline';
                            cLim{m}      = [];
                            majorGrid{m} = 'on';
                            minorGrid{m} = 'off';
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        
                    case {'resp_dilate1_actQ_actSgn' 'resp_original_actQ_actSgn'}
                        % polyLabel = strsplit(metric{m},'_'); polyLabel = polyLabel{2};
                        % indIm         = roi(i).polyMask{ismember(roi(i).polyLabel,polyLabel)};
                        % indSig        = false(size(roi(i).im.actP.im));
                        % indSig(indIm) = mafdr(roi(i).im.actP.im(indIm),'BHFDR',true)<0.05;
                        % indNeg        = roi(i).im.act.im(:,:,:,1)<0;
                        % indPos        = roi(i).im.act.im(:,:,:,1)>0;
                        % nTrial        = [1; roi(i).nTrial; 1];
                        
                        % response within the vessel ROI dilated by 1 voxel,
                        % including only active voxels based on SPMG2 activation detection,
                        % and segregated by sign of activation
                        im = permute(roi(i).im.resp.im,[4 1 2 3]);
                        tsNeg = im(:,indIm&indSig&indNeg);
                        tsPos = im(:,indIm&indSig&indPos);
                        if isfield(roi(i).im,'respSd') && ~isempty(roi(i).im.respSd)
                            imSd    = permute(roi(i).im.respSd.im,[4 1 2 3]);
                            tsNegEr = imSd(:,indIm&indSig&indNeg);
                            tsPosEr = imSd(:,indIm&indSig&indPos);
                        else
                            tsNegEr = [];
                            tsPosEr = [];
                        end
                        t     = ((0:size(tsNeg,1)-1).*roi(i).im.resp.dt)';
                        
                        roi(i).smr{m}.vec     = cat(3,permute(tsNeg,[1 3 2]),permute(tsPos,[1 3 2]));
                        roi(i).smr{m}.vecEr   = cat(3,permute(tsNegEr,[1 3 2]),permute(tsPosEr,[1 3 2]));
                        roi(i).smr{m}.nVox    = [nnz(indIm&indSig&indNeg) nnz(indIm&indSig&indPos)];
                        roi(i).smr{m}.nVoxRoi = [nnz(indIm)];
                        roi(i).smr{m}.neg     = cat(3,true(1,1,roi(i).smr{m}.nVox(1)),false(1,1,roi(i).smr{m}.nVox(2)));
                        roi(i).smr{m}.t       = t;
                        roi(i).smr{m}.nTrial  = nTrial;
                        roi(i).smr{m}.label   = {'neg','pos'};
                        roi(i).smr{m}.metric  = metric{m};
                        roi(i).smr{m}.info    = strjoin({'time' '???' 'vox'},' x ');

                        if ~isempty(H)
                            % if ~isempty(roi(i).smr{m}.vecEr)
                                ind   =  roi(i).smr{m}.neg;
                                t     = roi(i).smr{m}.t;
                                vecAv =      mean( roi(i).smr{m}.vec(  :,:,ind)    ,3);
                                vecEr = sqrt(mean( roi(i).smr{m}.vecEr(:,:,ind).^2 ,3));  % pool variance across voxels
                                % thi is a very conservative approach because squared deviations cannot average out
                                % (for accurate measure of error of the voxel-averaged responses, one needs to work from the residuals or perform another fit on the voxel-averaged full timeseries)
                                
                                hTs{i,m}(1) = shplot2(t,vecAv,vecEr,hA{m}(i));
                                
                                ind   = ~roi(i).smr{m}.neg;
                                t     = roi(i).smr{m}.t;
                                vecAv =      mean( roi(i).smr{m}.vec(  :,:,ind)    ,3);
                                vecEr = sqrt(mean( roi(i).smr{m}.vecEr(:,:,ind).^2 ,3));
                                
                                hTs{i,m}(2) = shplot2(t,vecAv,vecEr,hA{m}(i));
                            % else
                            %     dbstack; error('double check that');
                            %     hTs{i,m} = plot(hA(i),roi(i).ts{m}.t,roi(i).ts{m}.vec);
                            % end
                            yLim{m}      = [];
                            yScale{m}    = 'linear';
                            yLabel{m}    = 'MR signal change';
                            cLim{m}      = [];
                            majorGrid{m} = 'on';
                            minorGrid{m} = 'off';
                        end
                        
                    otherwise
                        dbstack; error('double check metric');
                end

                

                % hTs{i}.UserData.labels = roi(i).ts.label;

            end
        end
    else
        dbstack; error('double check roi data format');
        dataMask = MRIread(rCond.volMt.runAv.psd.param.fMask); dataMask = dataMask.vol~=0;
        f       = rCond.volMt.runAv.psd.f;
        spec    = size(rCond.volMt.runAv.psd.PSD,1:8); spec(6) = length(roi); spec = zeros(spec);
        for i = 1:length(roi)
            vec2roi = roi{i}.mask(dataMask);
            spec(:,:,:,:,:,i,:,:) = mean(rCond.volMt.runAv.psd.PSD(:,:,:,:,:,vec2roi,:,:),6);
        end
    end

    % %% Plot PSD
    % for i = 1:length(roi)
    %     plot(hA{i},squeeze(f),squeeze(spec(:,:,:,:,:,i,:,:)),'k');
    % end


    if ~isempty(H)
        for m = 1:length(metric)
            %% Adjust axes
            for i = 1:length(H)
                hA{m}(i).XAxis.Color = H(i).XAxis.Color; hA{m}(i).XAxis.LineWidth = H(i).XAxis.LineWidth;
                hA{m}(i).YAxis.Color = H(i).YAxis.Color; hA{m}(i).YAxis.LineWidth = H(i).YAxis.LineWidth;
            end
            axis(hA{m},'tight');
            
            if ismember(metric{m},{'respPhs_peakBasePhsInDilate1'})
                for i = 1:length(roi)
                    set(hA{m}(i),...
                        'YLim',yLim{m}{i},...
                        'YScale',yScale{m},...
                        'XGrid','on','YGrid','on',...
                        'XMinorGrid',minorGrid{m},'YMinorGrid',minorGrid{m},...
                        'GridColor',[0.5 0.5 0.5],'MinorGridColor',[0.5 0.5 0.5]);
                end
            else
                vInd = [roi.anot_sig];
                if contains(metric{m},{'coh_' 'cohTrialGram_'})
                    vInd = true(size(vInd));
                end
                if ~exist('yLim','var') || isempty(yLim{m})
                    yLim{m} = get(hA{m}(vInd),'YLim'); yLim{m} = [min([yLim{m}{:}]) max([yLim{m}{:}])];
                end
                set(hA{m}(vInd),...
                'YLim',yLim{m},...
                'YScale',yScale{m},...
                'XGrid','on','YGrid','on',...
                'XMinorGrid',minorGrid{m},'YMinorGrid',minorGrid{m},...
                'GridColor',[0.5 0.5 0.5],'MinorGridColor',[0.5 0.5 0.5]);
            end
            if strcmp(metric{m},'psdTrialGram_dilate1_actQ') || strcmp(metric{m},'cohTrialGram_dilate1')
                if ~exist('cLim','var') || isempty(cLim{m})
                    cLim{m} = get(hA{m}(vInd),'CLim'); cLim{m} = [min([cLim{m}{:}]) max([cLim{m}{:}])];
                end
                set(hA{m}(vInd),'CLim',cLim{m});
                hCb = colorbar(hA{m}(end),'Location','manual');
                hCb.Position = [sum(hA{m}(end).Position([1 3])), hA{m}(end).Position(2), 0.1*hA{m}(end).Position(3), hA{m}(end).Position(4)];
                if strcmp(metric{m},'psdTrialGram_dilate1_actQ')
                    if length(hCb.Ticks) == 1
                        tick100   = round(cLim{m}/100  ); tick1000(tick100   ==0) = 1; tick100(  tick100>10)   = 9; tick100   = (tick100(  1):tick100(  end))*  100;
                        tick1000  = round(cLim{m}/1000 ); tick1000(tick1000  ==0) = 1; tick1000( tick1000>10)  = 9; tick1000  = (tick1000( 1):tick1000( end))* 1000;
                        tick10000 = round(cLim{m}/10000); tick10000(tick10000==0) = 1; tick10000(tick10000>10) = 9; tick10000 = (tick10000(1):tick10000(end))*10000;
                        tick = [tick100 tick1000 tick10000];
                        tick(tick<cLim{m}(1)) = []; tick(tick>cLim{m}(2)) = [];
                        hCb.Ticks = tick;
                    end
                    ylabel(hCb, 'PSD');
                else
                    ylabel(hCb, 'coherence');
                end
                
            end
            
            drawnow;
        end

        for m = 1:length(metric)
            for i = 1:length(roi)
                % Add text annotations for voxel counts
                if ismember(metric{m},{'resp_dilate1_actQ_actSgn' 'resp_original_actQ_actSgn'})
                    nTrial = round(max(roi(i).smr{m}.nTrial(2:end-1)));
                else
                    nTrial = roi(i).R*6;
                end
                switch metric{m}
                    case {'resp_dilate1_actQ_actSgn' 'resp_original_actQ_actSgn' 'respPhs_peakBasePhsInDilate1'}
                        % Bottom left corner - total ROI voxel count
                        if ~ismember(metric{m},{'respPhs_peakBasePhsInDilate1'})
                            text(hA{m}(i), min(hA{m}(i).XLim)+range(hA{m}(i).XLim)*0.01, min(hA{m}(i).YLim)+range(hA{m}(i).YLim)*0.01, ...
                                [num2str(roi(i).smr{m}.nVoxRoi) 'vox'], ...
                                'HorizontalAlignment', 'left', ...
                                'VerticalAlignment', 'bottom', ...
                                'FontSize', 8);
                        end
                        % Top right corner - positive and significant voxel count
                        if ismember(metric{m},{'resp_dilate1_actQ_actSgn' 'resp_original_actQ_actSgn'})
                            text(hA{m}(i),...
                            min(hA{m}(i).XLim)+range(hA{m}(i).XLim)*0.99, min(hA{m}(i).YLim)+range(hA{m}(i).YLim)*0.99, ...
                            [num2str(roi(i).smr{m}.nVox(ismember(roi(i).smr{m}.label,'pos'))) 'posVox'], ...
                            'HorizontalAlignment', 'right', ...
                            'VerticalAlignment', 'top', ...
                            'FontSize', 8);
                        end
                        % Bottom right corner - negative and significant voxel count
                        if ismember(metric{m},{'resp_dilate1_actQ_actSgn' 'resp_original_actQ_actSgn'})
                            text(hA{m}(i), min(hA{m}(i).XLim)+range(hA{m}(i).XLim)*0.99, min(hA{m}(i).YLim)+range(hA{m}(i).YLim)*0.01, ...
                            [num2str(roi(i).smr{m}.nVox(ismember(roi(i).smr{m}.label,'neg'))) 'negVox'], ...
                            'HorizontalAlignment', 'right', ...
                            'VerticalAlignment', 'bottom', ...
                            'FontSize', 8);
                        end
                        % top left corner - number of trials
                        text(hA{m}(i), min(hA{m}(i).XLim)+range(hA{m}(i).XLim)*0.01, min(hA{m}(i).YLim)+range(hA{m}(i).YLim)*0.99, ...
                        [num2str(nTrial) 'trials'], ...
                        'HorizontalAlignment', 'left', ...
                        'VerticalAlignment', 'top', ...
                        'FontSize', 8);
                    case {'psd_dilate1_actQ' 'psd_original_actQ'}
                        % Bottom left corner - number of significant voxels over total ROI voxel count
                        text(hA{m}(i), ...
                            0.01, ...
                            0.01, ...
                            [num2str(roi(i).smr{m}.nVox) '/' num2str(roi(i).smr{m}.nVoxRoi) 'vox'], ...
                            'Units', 'normalized', ...
                            'HorizontalAlignment', 'left', ...
                            'VerticalAlignment', 'bottom', ...
                            'FontSize', 8);
                        % top right corner - number of trials
                        text(hA{m}(i), ...
                            0.99, 0.99, ...
                            [num2str(nTrial) 'trials'], ...
                            'Units', 'normalized', ...
                            'HorizontalAlignment', 'right', ...
                            'VerticalAlignment', 'top', ...
                            'FontSize', 8);
                    case {'psdTrialGram_dilate1_actQ' 'psdTrialGram_original_actQ'}
                        % Top left corner - number of significant voxels over total ROI voxel count
                        text(hA{m}(i), ...
                            0.01, ...
                            0.99, ...
                            [num2str(roi(i).smr{m}.nVox) '/' num2str(roi(i).smr{m}.nVoxRoi) 'vox'], ...
                            'Units', 'normalized', ...
                            'HorizontalAlignment', 'left', ...
                            'VerticalAlignment', 'top', ...
                            'FontSize', 8);
                        % top right corner - number of trials
                        text(hA{m}(i), ...
                            0.99, 0.99, ...
                            [num2str(nTrial) 'trials'], ...
                            'Units', 'normalized', ...
                            'HorizontalAlignment', 'right', ...
                            'VerticalAlignment', 'top', ...
                            'FontSize', 8);
                    case {'coh_dilate1' 'coh_original' 'cohTrialGram_dilate1' 'cohTrialGram_original'}
                        yline(hA{m}(i),1/roi(i).mt.svd.K,':k');
                        % bottom right corner - significance
                        if roi(i).anot_sig
                            text(hA{m}(i), ...
                                0.99, 0.01, ...
                                '*', ...
                                'Units', 'normalized', ...
                                'HorizontalAlignment', 'right', ...
                                'VerticalAlignment', 'bottom', ...
                                'FontSize', 8);
                        end
                        % top right corner - number of trials
                        text(hA{m}(i), ...
                            0.99, 0.99, ...
                            [num2str(nTrial) 'trials'], ...
                            'Units', 'normalized', ...
                            'HorizontalAlignment', 'right', ...
                            'VerticalAlignment', 'top', ...
                            'FontSize', 8);
                        % Bottom left corner - total ROI voxel count
                        text(hA{m}(i), ...
                            0.01, ...
                            0.01, ...
                            [num2str(roi(i).smr{m}.nVoxRoi) 'vox'], ...
                            'Units', 'normalized', ...
                            'HorizontalAlignment', 'left', ...
                            'VerticalAlignment', 'bottom', ...
                            'FontSize', 8);
                    otherwise
                        dbstack; error('double check metric');
                end
            end
        end

    end
    

