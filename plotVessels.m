function [hFigCat,hFigRun] = plotVessels(runCond,catFlag,runFlag)

if ~exist('catFlag','var'); catFlag = []; end
if ~exist('runFlag','var'); runFlag = []; end
if isempty(catFlag); catFlag = 0; end
if isempty(runFlag); runFlag = 0; end

[runCond.volResp.thresh] = deal(runCond.volActCat);

param.KperRun = [1 4 4]; % [gram gramMD full]
param.win     = 25;

plotFlag = runFlag;
for R = 1:size(runCond.volTs,1)
    disp(['run' num2str(R) '/' num2str(size(runCond.volTs,1))])
    [hFigRun{R},artRoi{R},veiRoi{R},ambRoi{R},vesRoi{R},imRespEr{R},imRespN{R}] = ...
        doIt(runCond.volTs(R,1),runCond.dsgn,runCond.volAct(R,1),runCond.volResp(R,1),runCond.phys(R,1),runCond.volAnat,param,plotFlag);
    % [artRoi{R},veiRoi{R},ambRoi{R},vesRoi{R},imRespEr{R},imRespN{R}] = ...
    %     doIt(runCond.volTs(R,1),runCond.dsgn,runCond.volAct(R,1),runCond.volResp(R,1),runCond.volAnat,param,plotFlag);
end
if catFlag
    plotFlag = 1;
    hFigCat = doIt(runCond.volTs,runCond.dsgn,runCond.volActCat,runCond.volRespCat,runCond.phys,runCond.volAnat,param,plotFlag,catFlag,artRoi,veiRoi,ambRoi,vesRoi,imRespEr,imRespN);
end


% function [artRoi,veiRoi,ambRoi,vesRoi,imRespEr,imRespN] = doIt(volTs,dsgn,volAct,volResp,volAnat,param,plotFlag,catFlag,artRoi,veiRoi,ambRoi,vesRoi,imRespEr,imRespN)
function [hFig,artRoi,veiRoi,ambRoi,vesRoi,imRespEr,imRespN] = doIt(volTs,dsgn,volAct,volResp,physTs,volAnat,param,plotFlag,catFlag,artRoi,veiRoi,ambRoi,vesRoi,imRespEr,imRespN)
if ~exist('plotFlag','var'); plotFlag = []; end
if ~exist('imRespEr','var'); imRespEr = []; end
if ~exist('imRespN' ,'var');  imRespN = []; end
if ~exist('artRoi'  ,'var');   artRoi = []; end
if ~exist('veiRoi'  ,'var');   veiRoi = []; end
if ~exist('ambRoi'  ,'var');   ambRoi = []; end
if ~exist('vesRoi'  ,'var');   vesRoi = []; end
if ~exist('catFlag' ,'var');  catFlag = []; end
roi.art = artRoi;
roi.vei = veiRoi;
roi.amb = ambRoi;
roi.ves = vesRoi;
if isempty(plotFlag); plotFlag = 1; end
if isempty(catFlag  );   catFlag = 0; end
[volTs.dsgn] = deal(dsgn);

imBase = MRIload3(volAct.fs.fBaseAv,[],[],0); imBase = imBase.vol;
mask   = MRIload3(volAct.fs.fMask ,[],[],0); mask = logical(mask.vol);
imResp = MRIload3(volResp.fs.fRespTs,[],[],0); tr = imResp.tr/1000; nFrame = imResp.nframes; imResp  = imResp.vol;
if isfield(volResp,'thresh')
    %use volResp.thresh (e.g. volAct across all runs) for voxel selection for response timecourse visualization
    imP    = MRIload3(volResp.thresh.fs.fFullP,[],[],0); imP  = imP.vol;
    imQ    = ones(size(imP)); imQ(mask) = mafdr(imP(mask));
    imAct  = getSPMG2coef(volResp.thresh.fs.fCoef,volAnat.calcarineVessel.f);
    imAct(:,:,:,:,end+1) = imQ<0.05; % second layer of 5th dimension is alpha channel
    imPVoxSel   = imP  ;
    imActVoxSel = imAct;
else
    %use run-specific volAct for voxel selection for response timecourse visualization
    imPVoxSel   = [];
    imActVoxSel = [];
end
imP    = MRIload3(volAct.fs.fFullP,[],[],0); imP  = imP.vol;
imQ    = ones(size(imP)); imQ(mask) = mafdr(imP(mask));
imAct  = getSPMG2coef(volAct.fs.fCoef,volAnat.calcarineVessel.f);
imAct(:,:,:,:,end+1) = imQ<0.05; % second layer of 5th dimension is alpha channe
if isempty(imPVoxSel  );   imPVoxSel = imP  ; end
if isempty(imActVoxSel); imActVoxSel = imAct; end
tResp  = linspace(0,tr.*(nFrame-1),nFrame)';

if isempty(imRespEr)
    for e = 1:size(volResp.fs.fRespEr,1)
        mriRespEr = MRIload3(volResp.fs.fRespEr{e},[],[],0);
        if e==1; imRespEr = zeros([size(mriRespEr.vol) size(volResp.fs.fRespEr,1)]); end
        imRespEr(:,:,:,:,e) = mriRespEr.vol;
    end
else
    imRespEr = cat(5,imRespEr{:});
end

if isempty(imRespN)
    if ~iscell(volResp.fs.fRespN)
        volResp.fs.fRespN = {volResp.fs.fRespN};
    end
    for R = 1:size(volResp.fs.fRespN,1)
        mriRespN = load(volResp.fs.fRespN{R});
        if R==1
            imRespN = zeros(size(imRespEr,1:4));
        end
        imRespN = imRespN + permute(sum(mriRespN.n,3),[2 3 4 1]);
    end
else
    imRespN = sum(cat(5,imRespN{:}),5);
end

% Extract vessel roi
imVes = MRIload3(volAnat.calcarineVessel.f,[],[],0); imVes = imVes.vol;
cropSize = 8;
artRoi = getVesselRoi(imVes==902         ,{'base' 'act' 'actP' 'resp' 'respEr' 'respN' 'voxSelAct' 'voxSelP'},{imBase imAct imP imResp imRespEr imRespN imActVoxSel imPVoxSel},cropSize);
veiRoi = getVesselRoi(imVes==914         ,{'base' 'act' 'actP' 'resp' 'respEr' 'respN' 'voxSelAct' 'voxSelP'},{imBase imAct imP imResp imRespEr imRespN imActVoxSel imPVoxSel},cropSize);
ambRoi = getVesselRoi(imVes==30|imVes==62,{'base' 'act' 'actP' 'resp' 'respEr' 'respN' 'voxSelAct' 'voxSelP'},{imBase imAct imP imResp imRespEr imRespN imActVoxSel imPVoxSel},cropSize);
vesRoi = getVesselRoi(imVes~=0           ,{'base' 'act' 'actP' 'resp' 'respEr' 'respN' 'voxSelAct' 'voxSelP'},{imBase imAct imP imResp imRespEr imRespN imActVoxSel imPVoxSel},cropSize);
vesRoi = vesRoi(end);


if ~plotFlag
    hF    = [];
    hFrsp = [];
else
    % Plot activation maps and vessels
    [hF,hT,hUL] = plotUL(imBase,{artRoi veiRoi ambRoi vesRoi},{'r' 'b' 'y' ''},[0 800]);
    [~,b,~] = fileparts(fileparts(volAct.fs.fBaseAv));
    title(hUL{1},b,'interpreter','none')
    % title(hUL{1},subList{S})
    hOL = plotOL(hUL,imAct,{artRoi veiRoi ambRoi vesRoi},{'r' 'b' 'y' ''},[-150 150]);
    drawnow

    % Plot response timecourse
    hFrsp = figure('WindowStyle','docked');
    roi = [artRoi; veiRoi; ambRoi; vesRoi];
    hRsp = cell(size(hUL));
    for r = 1:length(roi)
        % if strcmp(roi(r).label,'vesselAll')
        %     keyboard
        % end
        anatMsk = roi(r).mask(roi(r).im.resp.y(1):roi(r).im.resp.y(2),roi(r).im.resp.x(1):roi(r).im.resp.x(2));
        posMsk  = roi(r).im.voxSelAct.im(:,:,:,:,1)>0;
        negMsk  = roi(r).im.voxSelAct.im(:,:,:,:,1)<0;
        actMsk  = false(size(anatMsk));
        actMsk(anatMsk)  = mafdr(roi(r).im.voxSelP.im(anatMsk),'BHFDR',true) < 0.05;
        % actMsk  = roi(r).im.act.im(:,:,:,:,2);

        hRsp{r+1} = axes('Position',get(hUL{r+1},'Position'));
        rsp   = permute(roi(r).im.resp.im  ,[4 5 1 2 3]);
        rspEr = permute(roi(r).im.respEr.im,[4 5 1 2 3]);
        rspN  = permute(roi(r).im.respN.im ,[4 5 1 2 3]);

        posRsp    = mean(rsp  (:,:,posMsk & anatMsk & actMsk),3);
        posRspEr  = mean(rspEr(:,:,posMsk & anatMsk & actMsk),3);
        posRspN   = mean(rspN (:,:,posMsk & anatMsk & actMsk),3);
        posRspSem = sqrt(sum(posRspEr.^2,2)./(posRspN-1)) ./ sqrt(posRspN);
        posRspSem(isnan(posRspSem)) = 0;
        % plot(hRsp{r+1},tResp,posRsp,'r'); hold on
        H = shplot(tResp,posRsp,posRspSem,'Color','r');
        H.patch.FaceAlpha = 0.5;
        hold on

        negRsp    = mean(rsp  (:,:,negMsk & anatMsk & actMsk),3);
        negRspEr  = mean(rspEr(:,:,negMsk & anatMsk & actMsk),3);
        negRspN   = mean(rspN (:,:,negMsk & anatMsk & actMsk),3);
        negRspSem = sqrt(sum(negRspEr.^2,2)./(negRspN-1)) ./ sqrt(negRspN);
        negRspSem(isnan(negRspSem)) = 0;
        H = shplot(tResp,negRsp,negRspSem,'Color','b');
        H.patch.FaceAlpha = 0.5;

        axis tight; grid on; grid minor;

        hRsp{r+1}.XAxis.Color = hUL{r+1}.XAxis.Color;
        hRsp{r+1}.YAxis.Color = hUL{r+1}.YAxis.Color;
        hRsp{r+1}.XAxis.LineWidth = hUL{r+1}.XAxis.LineWidth;
        hRsp{r+1}.YAxis.LineWidth = hUL{r+1}.YAxis.LineWidth;

        x = 0.5; y = 1;
        tx = [num2str(nnz(anatMsk & actMsk)) '/' num2str(nnz(anatMsk))];
        text(hRsp{r+1},x,y,tx,'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','Color','k');
        x = 0.5; y = 0;
        tx = [num2str(nnz(posMsk & anatMsk & actMsk)) '/' num2str(nnz(anatMsk & actMsk))];
        hTxPos = text(hRsp{r+1},x,y,tx,'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom','Color','r');
        x = 0.5; y = 0;
        tx = [num2str(nnz(negMsk & anatMsk & actMsk)) '/' num2str(nnz(anatMsk & actMsk))];
        hTxNeg = text(hRsp{r+1},x,y,tx,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','Color','b');
        linkprop([hUL{r+1} hRsp{r+1}],'Position');
    end
    yLim = get([hRsp{:}],'YLim'); yLim = [-1 1].*max(abs([yLim{:}]));
    set([hRsp{:}],'YLim',yLim);
    xLim = get([hRsp{:}],'XLim'); xLim = [min([xLim{:}]) max([xLim{:}])];
    set([hRsp{:}],'XLim',xLim);
    posi = get([hUL{2:end}],'Position'); posi = cat(1,posi{:});
    posi = [false; abs(posi(:,1)-min(posi(:,1)))<0.001];
    ylabel([hRsp{posi}],'signal change');
    set([hRsp{~posi}],'YTickLabel',[]);
    set([hRsp{:}],'Box','on');
end

hFig = [];
% return


% Plot trial-tiggered time-freq mt analysis
% MT spectral ana
if ~plotFlag
    hFtf = [];
else
    roiList = {artRoi veiRoi ambRoi vesRoi};
    K     = param.KperRun .* size(volTs,1);
    W     = [];
    win   = param.win; % in seconds [lenght, step]
    if catFlag
        volTsCat = volTs(1);
        volTsCat.mri.vol = [];
        volTsCat.mri.t   = [];
        volTsCat.dsgn.onsetList = [];
        volTsCat.dsgn.ondurList = [];
        volTsCat.dsgn.nullTrial = [];
        for R = 1:size(volTs,1)
            mri = MRIload3(volTs(R).mri,[],[],0);
            t = linspace(0,mri.tr/1000.*(mri.nFrame-1),mri.nFrame)';
            t = t + mri.tr/1000.*mri.nDummyRemoved;
            t0 = round(seconds(mri.acqTime)/(mri.tr/1000))*(mri.tr/1000);
            t = t + t0;
            volTsCat.mri.vol = cat(4,volTsCat.mri.vol,mri.vol);
            volTsCat.mri.t = cat(1,volTsCat.mri.t,t);
            volTsCat.mri.nFrame(R) = mri.nFrame;
            volTsCat.dsgn.onsetList = cat(2,volTsCat.dsgn.onsetList,dsgn.onsetList + t0);
            volTsCat.dsgn.ondurList = cat(2,volTsCat.dsgn.ondurList,dsgn.ondurList     );
            volTsCat.dsgn.nullTrial = cat(2,volTsCat.dsgn.nullTrial,dsgn.nullTrial     );
        end
        % volTsCat.mri.nFrame = size(volTsCat.mri.t,1);
        volTsCat.mri.nframes = size(volTsCat.mri.t,1);
        volTs = volTsCat; clear volTsCat
        tmp = strsplit(volTs.mri.fspec,'_run-'); tmp{2} = strsplit(tmp{2},'_'); tmp{2}{1} = 'cat'; tmp{2} = strjoin(tmp{2},'_');
        volTs.mri.fspec = strjoin(tmp,'_run-');

    end
    extra.catFlag = catFlag;
    extra.padTo = [nan 1025 2049]; % [gram gramMD full]
    verboseThis = 0;
    rrr = 0; nRoi = sum(cellfun('size',roiList,1));
    volPsd = cell(size(roiList));
    for r = 1:length(roiList)
        for rr = 1:length(roiList{r})
            rrr = rrr + 1;
            % if rrr==nRoi
            %     keyboard
            % end
            disp(['roi' num2str(rrr) '/' num2str(nRoi)]);
            % volPsd{r}(rr) = runFullMT3(volTs,W,K,win,[],[],roiList{r}(rr).mask,extra,[],[],verboseThis,[],[])';
            volPsd{r}(rr) = runFullMT4(volTs,W,K,win,[],[],roiList{r}(rr).mask,extra,[],[],verboseThis,[],[])';
            roiList{r}(rr).psd = volPsd{r}(rr).psd;
            roiList{r}(rr).svd = volPsd{r}(rr).svd;
            roiList{r}(rr).psdTrialGramMD = volPsd{r}(rr).psdTrialGramMD;
            roiList{r}(rr).svdTrialGramMD = volPsd{r}(rr).svdTrialGramMD;
        end
    end
    artRoi = roiList{1};
    veiRoi = roiList{2};
    ambRoi = roiList{3};
    vesRoi = roiList{4};

    metricLabelList = {'psdEPC' 'cohEPC'};
    timeLabel   = 'trialGramMD';
    for i = 1:length(metricLabelList)
        metricLabel = metricLabelList{i};
        hFtf(i) = figure('WindowStyle','docked');
        roi = [artRoi; veiRoi; ambRoi; vesRoi];
        hTf = cell(size(hUL));
        for r = 1:length(roi)
            switch timeLabel
                case 'trialGramMD'
                    switch metricLabel
                        case 'coh'
                            roiSpec = roi(r).svdTrialGramMD;
                            tf = permute(     roiSpec.vec.coh(:,:,:,:,:,:,:,1)   ,[5 7 1 2 3 4 6 8]);
                            t  = permute(mean(roiSpec.t         (:,1,:,:,:,:,:,:),1),[5 7 1 2 3 4 6 8]) - roiSpec.onsetList(1);
                            f  = permute(     roiSpec.f         (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                            label = 'coherence';
                            cScale = 'linear';
                        case 'psd'
                            roiSpec = roi(r).psdTrialGramMD;
                            tf = permute(mean(roiSpec.vec.psd(:,:,:,:,:,:,:,1),6),[5 7 1 2 3 4 6 8]);
                            t  = permute(mean(roiSpec.t        (:,1,:,:,:,:,:,:),1),[5 7 1 2 3 4 6 8]) - roiSpec.onsetList(1);
                            f  = permute(     roiSpec.f        (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                            label = 'psd';
                            cScale = 'log';
                        case 'cohEPC'
                            roiSpec = roi(r).svdTrialGramMD;
                            tf = permute(     roiSpec.vec.cohEPC(:,:,:,:,:,:,:,1)   ,[5 7 1 2 3 4 6 8]);
                            t  = permute(mean(roiSpec.t         (:,1,:,:,:,:,:,:),1),[5 7 1 2 3 4 6 8]) - roiSpec.onsetList(1);
                            f  = permute(     roiSpec.f         (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                            label = 'phase-locked coherence';
                            cScale = 'linear';
                        case 'psdEPC'
                            roiSpec = roi(r).psdTrialGramMD;
                            tf = permute(mean(roiSpec.vec.psdPC(:,:,:,:,:,:,:,1),6),[5 7 1 2 3 4 6 8]);
                            t  = permute(mean(roiSpec.t        (:,1,:,:,:,:,:,:),1),[5 7 1 2 3 4 6 8]) - roiSpec.onsetList(1);
                            f  = permute(     roiSpec.f        (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                            label = 'phase-locked psd';
                            cScale = 'log';
                        otherwise
                            dbstack; error('X');
                    end
                otherwise
                    dbstack; error('X');
            end

            hTf{r+1} = axes('Position',get(hUL{r+1},'Position'));
            hTf{r+1}.XAxis.Color = hUL{r+1}.XAxis.Color; hTf{r+1}.XAxis.LineWidth = hUL{r+1}.XAxis.LineWidth;
            hTf{r+1}.YAxis.Color = hUL{r+1}.YAxis.Color; hTf{r+1}.YAxis.LineWidth = hUL{r+1}.YAxis.LineWidth;
            hTf{r+1}.Box = 'on';


            hIm = imagesc(hTf{r+1},t,f,tf);
            hTf{r+1}.ColorScale = cScale;
            axis tight;

            linkprop([hUL{r+1} hTf{r+1}],'Position');
        end
        yLim = get([hTf{:}],'YLim'); yLim = [min([yLim{:}]) max([yLim{:}])];
        set([hTf{:}],'YLim',yLim);
        if exist('hRsp','var')
            xLim = get([hRsp{:}],'XLim'); xLim = [min([xLim{:}]) max([xLim{:}])];
            set([hTf{:}],'XLim',xLim);
        else
            xLim = get([hTf{:}],'XLim'); xLim = [min([xLim{:}]) max([xLim{:}])];
            set([hTf{:}],'XLim',xLim);
        end
        cLim = get([hTf{:}],'CLim'); cLim = [min([cLim{:}]) max([cLim{:}])];
        set([hTf{:}],'CLim',cLim);

        hCb = colorbar(hTf{2},'location','manual');
        hCb.Position(3) = 0.01;
        hCb.Position(1) = 0.99;
        hCb.AxisLocation = 'in';
        ylabel(hCb,label)

        addWin(hTf{end},roiSpec);
        addW(hTf{end},roiSpec);

        posi = get([hUL{2:end}],'Position'); posi = cat(1,posi{:});
        posi = [false; abs(posi(:,1)-min(posi(:,1)))<0.001];
        ylabel([hTf{posi}],'Hz');
        set([hTf{~posi}],'YTickLabel',[]);

    end
end


% Analysis phys
if ~isempty(physTs)


    [physRoi,physTsz] = getPhysRoi(physTs,[],dsgn);
    verboseThis = 9;
    K   = 20;
    W   = [];
    win = inf;
    testFlag = 0;
    tic
    physPsd = runFullMT4(physTs ,W,K,win,[],[],physRoi.mask            ,[]   ,[],[],verboseThis,[],[],testFlag)';
    toc

    [physRoiX,physTsX] = getPhysRoi(physTs,[],dsgn);
    physPsd
    physRoi

    vesRoi = getPhysRoi(imVes~=0           ,{'base' 'act' 'actP' 'resp' 'respEr' 'respN' 'voxSelAct' 'voxSelP'},{imBase imAct imP imResp imRespEr imRespN imActVoxSel imPVoxSel},cropSize);
    
    % chanInd = 1;
    % figure('WindowStyle','docked'); ax = {};
    % ht = tiledlayout(2,3);
    % 
    % ax{end+1} = nexttile([1 1]);
    % plot(squeeze(physPsd.psd.PSD(1,1,1,1,:,chanInd,:)),squeeze(physPsd.psd.f))
    % ylim([0 25])
    % ylabel('Hz')
    % ax{end}.YDir = 'reverse';
    % 
    % ax{end+1} = nexttile([1 2]);
    % imagesc(squeeze(mean(physPsd.psdGram.t,1)),squeeze(physPsd.psdGram.f),squeeze(physPsd.psdGram.PSD(1,1,1,1,:,chanInd,:)))
    % ax{end}.ColorScale = 'log';
    % ylabel('Hz'); xlabel('t (sec)'); ylim([0 25]); xlim(physPsd.t([1 end]));
    % 
    % ax{end+1} = nexttile([1 1]);
    % 
    % ax{end+1} = nexttile([1 2]);
    % imagesc(squeeze(mean(physPsd.psdTrialGramMD.t(:,1,:,:,:,:,:,:),1)),squeeze(physPsd.psdTrialGramMD.f),squeeze(physPsd.psdTrialGramMD.vec.psdPC(1,1,1,1,:,chanInd,:)))
    % ax{end}.ColorScale = 'log';
    % ylabel('Hz'); xlabel('t (sec)'); ylim([0 25]); xlim(physPsd.t([1 end]));


end


% Plot spectra
if ~plotFlag
    hFspec = [];
else
    metricLabelList = {'psd' 'coh'};
    for i = 1:length(metricLabelList)
        metricLabel = metricLabelList{i};
        hFspec(i) = figure('WindowStyle','docked');
        physRoi = [];
        roi = [artRoi; veiRoi; ambRoi; vesRoi; physPsd];
        hTf = cell(size(hUL));


        for r = 1:length(roi)
            if roi
                hTf{r+1} = axes('Position',get(hUL{r+1},'Position'));
                hTf{r+1}.XAxis.Color = hUL{r+1}.XAxis.Color; hTf{r+1}.XAxis.LineWidth = hUL{r+1}.XAxis.LineWidth;
                hTf{r+1}.YAxis.Color = hUL{r+1}.YAxis.Color; hTf{r+1}.YAxis.LineWidth = hUL{r+1}.YAxis.LineWidth;
                hTf{r+1}.Box = 'on'; hold on

                %%% Plot full time series spectrum
                switch metricLabel
                    case 'coh'
                        roiSpec = roi(r).svd;
                        tf = permute(     roiSpec.COH(:,:,:,:,:,:,:,1)   ,[5 7 1 2 3 4 6 8]);
                        f  = permute(     roiSpec.f  (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                        label = 'coherence';
                        cScale = 'linear';
                    case 'psd'
                        roiSpec = roi(r).psd;
                        tf = permute(mean(roiSpec.PSD(:,:,:,:,:,:,:,1),6),[5 7 1 2 3 4 6 8]);
                        f  = permute(     roiSpec.f  (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                        label = 'psd';
                        cScale = 'log';
                    otherwise
                        dbstack; error('X');
                end
                hFull = plot(hTf{r+1},f,tf,'k');
                roiSpec1 = roiSpec;


                %%% Plot spetrogram time points
                phaseAvList = {'PC'};
                for iii = 1:length(phaseAvList)
                    switch metricLabel
                        case 'coh'
                            roiSpec = roi(r).svdTrialGramMD;
                            tMid = squeeze(mean(roiSpec.t(:,1,:,:,:,:,:),1)) - roiSpec.onsetList(1);
                            [~,b1] = min(abs(tMid-5));
                            tMid(b1); % time window centerd on response peak
                            tStr = tMid - roiSpec.win(1)/2;
                            [~,b2] = min(abs(tStr));
                            tMid(b2); % time window starting at stim onset
                            b3 = length(tMid);
                            tMid(b3); % last time window, most devoid of direct response

                            switch phaseAvList{iii}
                                case 'PC'
                                    tf          = permute(     roiSpec.vec.cohEPC(:,:,:,:,:,:,[b1 b2 b3],1)   ,[5 7 1 2 3 4 6 8]);
                                    tf(:,end+1) = permute(mean(roiSpec.vec.cohEPC(:,:,:,:,:,:,:         ,1),7),[5 7 1 2 3 4 6 8]); % average over all time windows
                                case 'nPC'
                                    tf          = permute(     roiSpec.vec.coh   (:,:,:,:,:,:,[b1 b2 b3],1)   ,[5 7 1 2 3 4 6 8]);
                                    tf(:,end+1) = permute(mean(roiSpec.vec.coh   (:,:,:,:,:,:,:         ,1),7),[5 7 1 2 3 4 6 8]); % average over all time windows
                            end
                            f  = permute(     roiSpec.f         (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                            label = ['coherence_' phaseAvList{iii}];
                            cScale = 'linear';
                        case 'psd'
                            roiSpec = roi(r).psdTrialGramMD;
                            tMid = squeeze(mean(roiSpec.t(:,1,:,:,:,:,:),1)) - roiSpec.onsetList(1);
                            [~,b1] = min(abs(tMid-5));
                            tMid(b1); % time window centerd on response peak
                            tStr = tMid - roiSpec.win(1)/2;
                            [~,b2] = min(abs(tStr));
                            tMid(b2); % time window starting at stim onset
                            b3 = length(tMid);
                            tMid(b3); % last time window, most devoid of direct response

                            switch phaseAvList{iii}
                                case 'PC'
                                    tf          = permute(     mean(roiSpec.vec.psdPC(:,:,:,:,:,:,[b1 b2 b3],1),6)   ,[5 7 1 2 3 4 6 8]);
                                    tf(:,end+1) = permute(mean(mean(roiSpec.vec.psdPC(:,:,:,:,:,:,:         ,1),6),7),[5 7 1 2 3 4 6 8]); % average over all time windows
                                case 'nPC'
                                    tf          = permute(     mean(roiSpec.vec.psd  (:,:,:,:,:,:,[b1 b2 b3],1),6)   ,[5 7 1 2 3 4 6 8]);
                                    tf(:,end+1) = permute(mean(mean(roiSpec.vec.psd  (:,:,:,:,:,:,:         ,1),6),7),[5 7 1 2 3 4 6 8]); % average over all time windows
                            end
                            f           = permute(          roiSpec.f        (:,1,:,:,:,:,:         ,:)      ,[5 7 1 2 3 4 6 8]);
                            label = ['psd_' phaseAvList{iii}];
                            cScale = 'log';
                        otherwise
                            dbstack; error('X');
                    end
                    switch phaseAvList{iii}
                        case 'PC'
                            hGram_PC = plot(hTf{r+1},f,tf);
                            set(hGram_PC([1 2 4]),'Visible','off')
                            set(hGram_PC(3),'Color','b')
                        case 'nPC'
                            hGram_nPC = plot(hTf{r+1},f,tf);
                            set(hGram_nPC([1 2 4]),'Visible','off')
                            set(hGram_nPC(3),'Color','r')
                    end
                    % set(hGram(4),'Color','r')
                end
                roiSpec2 = roiSpec;

                hTf{r+1}.YScale = cScale;
                axis tight;
                switch metricLabel
                    case 'psd'
                        hL = findobj(hTf{r+1}.Children,'type','line');
                        vis = get(hL,'Visible');
                        hL = hL([vis{:}]==1);
                        yLim = ylim(hTf{r+1}); y = yLim(2);
                        for iiiii = 1:length(hL)
                            y = min(y,min(hL(iiiii).YData(hL(iiiii).XData>0.1)));
                        end
                        yLim(1) = y;
                        ylim(hTf{r+1},yLim)
                end
                grid on; grid minor;
                linkprop([hUL{r+1} hTf{r+1}],'Position');
            end
            yLim = get([hTf{:}],'YLim'); yLim = [min([yLim{:}]) max([yLim{:}])];
            set([hTf{:}],'YLim',yLim);
            posi = get([hUL{2:end}],'Position'); posi = cat(1,posi{:});
            posi = [false; abs(posi(:,1)-min(posi(:,1)))<0.001];
            ylabel([hTf{posi}],metricLabel)
            set([hTf{~posi}],'YTickLabel',[])

            hW = addW(hTf{end},roiSpec1);
            hW.Color = hFull.Color;
            switch metricLabel
                case 'coh'
                    hW.YData = mean(hTf{end}.YLim).*[1 1];
                case 'psd'
                    hW.YData = exp(mean(log(hTf{end}.YLim))).*[1 1];
            end
            hW(2) = addW(hTf{end},roiSpec2);
            hW(2).Color = hGram_PC(3).Color;
            switch metricLabel
                case 'coh'
                    hW(2).YData = mean(hTf{end}.YLim).*[1.15 1.15];
                case 'psd'
                    hW(2).YData = exp(mean(log(hTf{end}.YLim))).*[1.15 1.15];
            end
            hW(3) = addW(hTf{end},roiSpec2);
            hW(3).Color = hGram_nPC(4).Color;
            switch metricLabel
                case 'coh'
                    hW(3).YData = mean(hTf{end}.YLim).*[1.3 1.3];
                case 'psd'
                    hW(3).YData = exp(mean(log(hTf{end}.YLim))).*[1.3 1.3];
            end
            uistack(hW(1),'top')






            legend(hTf{end},[hFull; hGram_PC(3); hGram_nPC(3)],{'full' 'late-PC' 'late-nPC'},'box','off')




            % if exist('hRsp','var')
            %     xLim = get([hRsp{:}],'XLim'); xLim = [min([xLim{:}]) max([xLim{:}])];
            %     set([hTf{:}],'XLim',xLim);
            % else
            %     xLim = get([hTf{:}],'XLim'); xLim = [min([xLim{:}]) max([xLim{:}])];
            %     set([hTf{:}],'XLim',xLim);
            % end
            % cLim = get([hTf{:}],'CLim'); cLim = [min([cLim{:}]) max([cLim{:}])];
            % set([hTf{:}],'CLim',cLim);

            % hCb = colorbar(hTf{2},'location','manual');
            % hCb.Position(3) = 0.01;
            % hCb.Position(1) = 0.99;
            % hCb.AxisLocation = 'in';
            % ylabel(hCb,label)


            % physio
            % Get physio
            physTs.mri
            % physRoi.
            physRoi = getPhysRoi([],{'base' 'act' 'actP' 'resp' 'respEr' 'respN' 'voxSelAct' 'voxSelP'},{imBase imAct imP imResp imRespEr imRespN imActVoxSel imPVoxSel});
            vesRoi = vesRoi(end);





        end




    end
end

hFig = [hF hFrsp hFtf hFspec];
