function [hFigCat,hFigRun] = plotVessels3(runCond,catFlag,runFlag)

    if ~exist('catFlag','var'); catFlag = []; end
    if ~exist('runFlag','var'); runFlag = []; end
    if ~exist('subData','var'); subData = []; end
    if isempty(catFlag); catFlag = 0    ; end
    if isempty(runFlag); runFlag = 0    ; end
    if isempty(subData); subData = 'mag'; end

    
    
    
runFlag = 1;
catFlag = 1;


% [runCond.volResp.(subData).thresh] = deal(runCond.volActCat);

% param.KperRun = [1 7 5 5]; % [gram gramMD full]
% param.win     = 25;
param = [];

plotFlag = runFlag;
for R = 1:size(runCond.volTs,1)
    runCond.R = size(runCond.fPreprocList,1);
    runCond.r = R;
    disp(['run' num2str(R) '/' num2str(size(runCond.volTs,1))])
    
    [hFigRun{R},artRoi{R},veiRoi{R},ambRoi{R},vesRoi{R},imRespEr{R},imRespN{R}] = ...
        doIt(runCond.volTs(R,1),runCond.dsgn,runCond.volResp.(subData).actRun(R,1),runCond.volResp.(subData).respRun(R,1),runCond.volMt.run(R,1),runCond.phs,runCond.volAnat,param,plotFlag);
    % [artRoi{R},veiRoi{R},ambRoi{R},vesRoi{R},imRespEr{R},imRespN{R}] = ...
    %     doIt(runCond.volTs(R,1),runCond.dsgn,runCond.volAct(R,1),runCond.volResp(R,1),runCond.volAnat,param,plotFlag);
end
if catFlag
    plotFlag = 1;
    hFigCat = doIt(runCond.volTs,runCond.dsgn,runCond.volActCat,runCond.volRespCat,runCond.phys,runCond.volAnat,param,plotFlag,catFlag,artRoi,veiRoi,ambRoi,vesRoi,imRespEr,imRespN);
end


% function [artRoi,veiRoi,ambRoi,vesRoi,imRespEr,imRespN] = doIt(volTs,dsgn,volAct,volResp,volAnat,param,plotFlag,catFlag,artRoi,veiRoi,ambRoi,vesRoi,imRespEr,imRespN)
function [hFig,artRoi,veiRoi,ambRoi,vesRoi,imRespEr,imRespN] = doIt(volTs,dsgn,volAct,volResp,volMt,physTs,volAnat,param,plotFlag,catFlag,artRoi,veiRoi,ambRoi,vesRoi,imRespEr,imRespN)
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

fAct  = volAct.stats;  if isfield(fAct,'stats');   fAct = fAct.stats; end
fResp = volResp.stats; if isfield(fResp,'stats'); fResp = fResp.stats; end
imBase = MRIread(char(fAct.fPoly0Base)); imBase = imBase.vol;
mask   = MRIread(char(fAct.fMask));        mask = mask.vol;
imResp = MRIread(char(fResp.fResp)); trResp = imResp.tr/1000; nFrameResp = imResp.nframes; imResp = imResp.vol;

% imBase = MRIload3(volAct.fs.fBaseAv,[],[],0); imBase = imBase.vol;
% mask   = MRIload3(volAct.fs.fMask ,[],[],0); mask = logical(mask.vol);
% imResp = MRIload3(volResp.fs.fRespTs,[],[],0); tr = imResp.tr/1000; nFrame = imResp.nframes; imResp  = imResp.vol;
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
imP = MRIread(char(fAct.fCondF_pVal)); imP = imP.vol;
imQ = ones(size(imP)); imQ(mask) = mafdr(imP(mask));

[amp,hFig] = getSPMG2coef(fAct.fCondCoef,mask,plotFlag);

% imP    = MRIload3(volAct.fs.fFullP,[],[],0); imP  = imP.vol;
% imQ    = ones(size(imP)); imQ(mask) = mafdr(imP(mask));
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



% Analysis phys
if 1
    if ~isempty(physTs)
        
        % Catenate runs then compute psd using missing data tapers (not working because of )
        % [physRoi,physTs] = getPhysRoi(physTs,[],dsgn);
        % if size(physTs,1)>1
        %     %downsample
        %     t0 = mean(physRoi(1).ts.t0);
        %     for r = 1:size(physRoi,1)
        %         Fs = mean(physRoi(r).ts.Fs);
        %         n  = physRoi(r).ts.nFrame;
        %         ts = mean(physRoi(r).ts.t0) - t0;
        %         te = ts + (n-1)/Fs;
        % 
        %         Fs = 20;
        %         ts = round(ts*Fs)/Fs;
        %         te = round(te*Fs)/Fs;
        %         n  = (te-ts)*Fs;
        % 
        %         t = linspace(ts,te,(n-1))';
        %         physRoi(r).ts.vec = interp1(physRoi(r).ts.t,physRoi(r).ts.vec,t,[],'extrap');
        %         physRoi(r).ts.Fs = Fs;
        %         physRoi(r).ts.nFrame = n;
        %         physRoi(r).ts.nframes = n;
        %         physRoi(r).ts.t = t;
        %         physRoi(r).ts.tr = 1/physRoi(r).ts.Fs;
        %     end
        %     %catenate
        %     physRoiCat = physRoi(1);
        %     ts = [physRoi.ts];
        %     physRoiCat.ts.t0      =       cat(4,ts.t0     );
        %     physRoiCat.ts.Fs      = mean( cat(4,ts.Fs     ),4);
        %     physRoiCat.ts.info    =       cat(4,ts.info   );
        %     physRoiCat.ts.nFrame  = sum(  cat(4,ts.nFrame ),4);
        %     physRoiCat.ts.imMean  = mean( cat(4,ts.imMean ),4);
        %     physRoiCat.ts.t       =       cat(1,ts.t      );
        %     physRoiCat.ts.nDummy  =       cat(4,ts.nFrame );
        %     physRoiCat.ts.tr      = 1./physRoiCat.ts.Fs;
        %     physRoiCat.ts.nframes = sum(  cat(4,ts.nframes),4);
        %     physRoiCat.ts.vol     =       cat(4,ts.vol    );
        %     physRoiCat.ts.vol2vec;
        %     physRoiCat.ts.vec     =       cat(1,ts.vec    );
        %     physRoi = physRoiCat; clear ts physRoiCat
        % end

        % Compute psd run by run then combine
        for r = 1:size(physTs,1)
            [physRoi(r,1),~] = getPhysRoi(physTs(r,1),[],dsgn);
        end
        if length(physTs) > 1
            physTs


            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            nLim = 100000; % 224000

            t0 = cat(1,physTs.t0); t0 = mean(t0,2); t0 = t0-t0(1);
            nOrig  = cellfun('size',{physTs.vec},1)';
            FsOrig = cat(1,physTs.Fs); FsOrig = mean(FsOrig,2); FsOrig = mean(FsOrig,1);


            Fs = FsOrig;
            n  = nOrig;
            
            ind0 = round(t0*Fs);
            ind  = [];
            for i = 1:length(n)
                ind = cat(2,ind,ind0(i):(ind0(i)+n(i)-1));
            end
            t = ind/Fs;

            indOrig = ind;
            tOrig   = t';



            fac = nLim/sum(nOrig);
            nNew  = round(nOrig*fac);
            FsNew = FsOrig*fac;

            
            Fs = FsNew;
            n  = nNew;
            
            ind0 = round(t0*Fs);
            ind  = [];
            for i = 1:length(n)
                ind = cat(2,ind,ind0(i):(ind0(i)+n(i)-1));
            end
            t = ind/Fs;

            indNew = ind;
            tNew   = t';


            vec = cat(1,physTs.vec);
            Vq = interp1(tOrig,vec,tNew,'linear');
            Vq = interp1(tOrig,vec,tNew,'pchip');
            Vq = interp1(tOrig,vec,tNew,'cubic');
            Vq = interp1(tOrig,vec,tNew,'makima');
            Vq = interp1(tOrig,vec,tNew,'spline');

            figure('WindowStyle','docked');
            hold off
            plot(tOrig,vec(:,1))
            hold on
            plot(tNew,Vq(:,1))
            xlim(xLim)



            t   = tNew;
            Fs  = FsNew;

            K   = 5;
            tr  = 1/Fs;
            N   = length(t(:));
            pad = 1;
            tic
            [tp,eigs,tpNorm,Nx,padX] = getTapers(K,tr,N,t,pad);
            toc
            plot(t,tp,'.')






            






            nLim/3
            
            round(fac*sum(n))/sum(n)
            nLim/fac





            fac = 0.2;

            n  = cat(1,ts.nFrame);
            n  = mode(n).*ones(size(n));
            n  = round(n*fac);


            Fs = cat(1,ts.Fs); Fs = Fs(:,1);
            Fs = mean(Fs)*fac;

            ind0 = round(t0*Fs);
            ind  = [];
            for i = 1:length(n)
                ind = cat(1,ind,ind0(i):(ind0(i)+n(i)-1));
            end
            ind = ind';
            t   = ind/Fs;
            plot(t(:),t(:),'.')
            axis tight


            K = 5;
            tr = 1;
            N = length(t(:));
            pad = 1;
            [tp,eigs,tpNorm,Nx,padX] = getTapers(K,tr,N,t(:),pad);
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













            ts = [physRoi.ts]'; [physRoi.ts] = deal([]);
            % equalize number of frames
            nFrame = min([ts.nFrame]);
            for r = 1:length(ts)
                ts(r).nFrame  = nFrame;
                ts(r).t       = ts(r).t(1:nFrame,:);
                ts(r).nframes = nFrame;
                ts(r).vec     = ts(r).vec(1:nFrame,:);
            end
            % downsample
            Fs    = cat(4,ts.Fs); Fs = mean(Fs(:));
            FsNew = 10;
            t0    = mean(cat(4,ts.t0),2); t0 = t0 - t0(1);
            % t0 = cumsum([ts.nFrame]/Fs + 60) - (ts(1).nFrame/Fs + 60);
            t0new = round(t0*FsNew)/FsNew;
            
            
            for r = 1:length(ts)
                n    = ts(r).nFrame;
                nNew = round(n/Fs*FsNew);
                t    = ts(r).t + t0new(r);
                tNew = linspace(t(1),t(end) + 1/Fs - 1/FsNew,nNew)';
                ts(r).vec     = interp1(t,ts(r).vec,tNew);
                ts(r).t       = tNew;
                ts(r).Fs(:)   = FsNew;
                ts(r).tr      = 1/mean(ts(r).Fs);
                ts(r).nFrame  = nNew;
                ts(r).nframes = nNew;
            end
            %catenate
            ts(1).t0      =      cat(4,ts.t0     )   ; [ts(2:end).t0     ] = deal([]);
            % t0 = mean(ts(1).t0,2); t0 = t0-t0(1);
            ts(1).Fs      =      cat(4,ts.Fs     )   ; [ts(2:end).Fs     ] = deal([]);
            ts(1).nFrame  =  sum(cat(4,ts.nFrame ),4); [ts(2:end).nFrame ] = deal([]);
            ts(1).dsgn    =      cat(4,ts.dsgn   )   ; [ts(2:end).dsgn   ] = deal([]);
            ts(1).t       =      cat(1,ts.t      )   ; [ts(2:end).t      ] = deal([]);
            ts(1).tr      = mean(cat(4,ts.tr     ),4); [ts(2:end).tr     ] = deal([]);
            ts(1).nframes =  sum(cat(4,ts.nframes),4); [ts(2:end).nframes] = deal([]);
            ts(1).vec     =      cat(1,ts.vec    )   ; [ts(2:end).vec    ] = deal([]);
            ts     (2:end) = [];
            physRoi(2:end) = [];
            physRoi.ts = ts; clear ts;
        end

        verboseThis = 9;
        K   = param.KperRun(4)*length(physTs);
        % K   = 5*length(physTs);
        W   = [];
        win = inf;
        testFlag = 0;
        for r = 1:length(physRoi)
            physRoi(r,1).fs = runFullMT4(physRoi(r,1),W,K,win,[],[],[],[],1,[],verboseThis,[],[],testFlag)';
        end
        if size(physTs,1)>1
            PSD = [physRoi.fs]; PSD = [PSD.psd];
            physRoi(1).fs.psd.PSD = mean(cat(3,PSD.PSD),3);
            physRoi(1).fs.psd.dim(3) = size(physRoi(1).fs.psd.PSD,3);
            linePwr = [physRoi.fs]; linePwr = [linePwr.harm];
            physRoi(1).fs.harm.linePwr = mean(cat(3,linePwr.linePwr),3);
            physRoi(1).fs.harm.lineF = [];
            physRoi(1).fs.harm.lineP = [];
            physRoi(2:end) = [];
        end
    end
else
    physRoi = [];
end



if ~plotFlag
    hF    = [];
    hFrsp = [];
else
    % Plot activation maps and vessels
    [hF,hT,hUL] = plotUL2(imBase,{artRoi veiRoi ambRoi vesRoi physRoi},{'r' 'b' 'y' 'k' 'g'},[0 800]);
    [~,b,~] = fileparts(fileparts(volAct.fs.fBaseAv));
    title(hUL{1},b,'interpreter','none')
    % title(hUL{1},subList{S})
    hOL = plotOL(hUL,imAct,{artRoi veiRoi ambRoi vesRoi physRoi},{'r' 'b' 'y' 'k' 'g'},[-150 150]);
    drawnow

    % Plot response timecourse
    hFrsp = figure('WindowStyle','docked');
    roi = [artRoi; veiRoi; ambRoi; vesRoi; physRoi];
    hRsp = cell(size(hUL));
    for r = 1:length(roi)
        if strcmp(roi(r).label,'phys'); continue; end
        % if strcmp(roi(r).label,'phys')
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
        volTs.mri.dsgn = volTs.dsgn;
    end
    extra.catFlag = catFlag;
    extra.padTo = [nan 1025 2049]; % [gram gramMD full]
    verboseThis = 0;
    rrr = 0; nRoi = sum(cellfun('size',roiList,1));
    volPsd = cell(size(roiList));
    for r = 1:length(roiList)
        for rr = 1:length(roiList{r})
            rrr = rrr + 1;
            disp(['roi' num2str(rrr) '/' num2str(nRoi)]);
            roiList{r}(rr).fs = runFullMT4(volTs,W,K,win,[],[],roiList{r}(rr).mask,extra,[],[],verboseThis,[],[])';
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
        roi = [artRoi; veiRoi; ambRoi; vesRoi; physRoi];
        hTf = cell(size(hUL));
        for r = 1:length(roi)
            if strcmp(roi(r).label,'phys'); continue; end
            % hTf{r+1} = axes('Position',get(hUL{r+1},'Position'));
            % hTf{end}.Visible = 'off';

            switch timeLabel
                case 'trialGramMD'
                    switch metricLabel
                        case 'coh'
                            roiSpec = roi(r).fs.svdTrialGramMD;
                            tf = permute(     roiSpec.vec.coh(:,:,:,:,:,:,:,1)   ,[5 7 1 2 3 4 6 8]);
                            t  = permute(mean(roiSpec.t         (:,1,:,:,:,:,:,:),1),[5 7 1 2 3 4 6 8]) - roiSpec.onsetList(1);
                            f  = permute(     roiSpec.f         (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                            label = 'coherence';
                            cScale = 'linear';
                        case 'psd'
                            roiSpec = roi(r).fs.psdTrialGramMD;
                            tf = permute(mean(roiSpec.vec.psd(:,:,:,:,:,:,:,1),6),[5 7 1 2 3 4 6 8]);
                            t  = permute(mean(roiSpec.t        (:,1,:,:,:,:,:,:),1),[5 7 1 2 3 4 6 8]) - roiSpec.onsetList(1);
                            f  = permute(     roiSpec.f        (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                            label = 'psd';
                            cScale = 'log';
                        case 'cohEPC'
                            roiSpec = roi(r).fs.svdTrialGramMD;
                            tf = permute(     roiSpec.vec.cohEPC(:,:,:,:,:,:,:,1)   ,[5 7 1 2 3 4 6 8]);
                            t  = permute(mean(roiSpec.t         (:,1,:,:,:,:,:,:),1),[5 7 1 2 3 4 6 8]) - roiSpec.onsetList(1);
                            f  = permute(     roiSpec.f         (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                            label = 'phase-locked coherence';
                            cScale = 'linear';
                        case 'psdEPC'
                            roiSpec = roi(r).fs.psdTrialGramMD;
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


% Plot spectra
if ~plotFlag
    hFspec = [];
else
    metricLabelList = {'psd' 'coh'};
    for i = 1:length(metricLabelList)
        hFspec(i) = figure('WindowStyle','docked');
        roi = [artRoi; veiRoi; ambRoi; vesRoi; physRoi];
        hTf = cell(size(hUL));


        for r = 1:length(roi)
            if strcmp(roi(r).label,'phys')
                metricLabel = 'psd';
            else
                metricLabel = metricLabelList{i};
            end


            % if strcmp(roi(r).label,'phys'); keyboard; end
            hTf{r+1} = axes('Position',get(hUL{r+1},'Position'));
            hTf{r+1}.XAxis.Color = hUL{r+1}.XAxis.Color; hTf{r+1}.XAxis.LineWidth = hUL{r+1}.XAxis.LineWidth;
            hTf{r+1}.YAxis.Color = hUL{r+1}.YAxis.Color; hTf{r+1}.YAxis.LineWidth = hUL{r+1}.YAxis.LineWidth;
            hTf{r+1}.Box = 'on'; hold on

            %%% Plot full time series spectrum
            switch metricLabel
                case 'coh'
                    roiSpec = roi(r).fs.svd;
                    tf = permute(     roiSpec.COH(:,:,:,:,:,:,:,1)   ,[5 6 7 1 2 3 4 8]);
                    f  = permute(     roiSpec.f  (:,1,:,:,:,:,:,:)   ,[5 6 7 1 2 3 4 8]);
                    label = 'coherence';
                    cScale = 'linear';
                case 'psd'
                    roiSpec = roi(r).fs.psd;
                    if strcmp(roi(r).label,'phys')
                        tf = permute(roiSpec.PSD(:,:,:,:,:,:,:,1),[5 6 7 1 2 3 4 8]);
                    else
                        tf = permute(mean(roiSpec.PSD(:,:,:,:,:,:,:,1),6),[5 6 7 1 2 3 4 8]);
                    end
                    f  = permute(     roiSpec.f  (:,1,:,:,:,:,:,:)   ,[5 6 7 1 2 3 4 8]);
                    label = 'psd';
                    cScale = 'log';
                otherwise
                    dbstack; error('X');
            end
            if size(tf,2)==1
                hFull = plot(hTf{r+1},f,tf,'k');
            else
                hFull = plot(hTf{r+1},f,tf);
            end
            roiSpec1 = roiSpec;


            %%% Plot spetrogram time points
            if ~strcmp(roi(r).label,'phys')
                phaseAvList = {'PC'};
                for iii = 1:length(phaseAvList)
                    switch metricLabel
                        case 'coh'
                            roiSpec = roi(r).fs.svdTrialGramMD;
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
                                    % tf          = permute(     roiSpec.vec.cohEPC(:,:,:,:,:,:,[b1 b2 b3],1)   ,[5 7 1 2 3 4 6 8]);
                                    % tf(:,end+1) = permute(mean(roiSpec.vec.cohEPC(:,:,:,:,:,:,:         ,1),7),[5 7 1 2 3 4 6 8]); % average over all time windows
                                    tf          = permute(     roiSpec.vec.cohEPC(:,:,:,:,:,:,b3,1)   ,[5 7 1 2 3 4 6 8]);
                                case 'nPC'
                                    % tf          = permute(     roiSpec.vec.coh   (:,:,:,:,:,:,[b1 b2 b3],1)   ,[5 7 1 2 3 4 6 8]);
                                    % tf(:,end+1) = permute(mean(roiSpec.vec.coh   (:,:,:,:,:,:,:         ,1),7),[5 7 1 2 3 4 6 8]); % average over all time windows
                                    tf          = permute(     roiSpec.vec.coh   (:,:,:,:,:,:,b3,1)   ,[5 7 1 2 3 4 6 8]);
                            end
                            f  = permute(     roiSpec.f         (:,1,:,:,:,:,:,:)   ,[5 7 1 2 3 4 6 8]);
                            label = ['coherence_' phaseAvList{iii}];
                            cScale = 'linear';
                        case 'psd'
                            roiSpec = roi(r).fs.psdTrialGramMD;
                            tMid = squeeze(mean(roiSpec.t(:,1,:,:,:,:,:),1)) - roiSpec.onsetList(1);
                            % [~,b1] = min(abs(tMid-5));
                            % tMid(b1); % time window centerd on response peak
                            % tStr = tMid - roiSpec.win(1)/2;
                            % [~,b2] = min(abs(tStr));
                            % tMid(b2); % time window starting at stim onset
                            b3 = length(tMid);
                            tMid(b3); % last time window, most devoid of direct response

                            switch phaseAvList{iii}
                                case 'PC'
                                    % tf          = permute(     mean(roiSpec.vec.psdPC(:,:,:,:,:,:,[b1 b2 b3],1),6)   ,[5 7 1 2 3 4 6 8]);
                                    % tf(:,end+1) = permute(mean(mean(roiSpec.vec.psdPC(:,:,:,:,:,:,:         ,1),6),7),[5 7 1 2 3 4 6 8]); % average over all time windows
                                    tf          = permute(     mean(roiSpec.vec.psdPC(:,:,:,:,:,:,b3,1),6)   ,[5 6 7 1 2 3 4 8]);
                                case 'nPC'
                                    % tf          = permute(     mean(roiSpec.vec.psd  (:,:,:,:,:,:,[b1 b2 b3],1),6)   ,[5 7 1 2 3 4 6 8]);
                                    % tf(:,end+1) = permute(mean(mean(roiSpec.vec.psd  (:,:,:,:,:,:,:         ,1),6),7),[5 7 1 2 3 4 6 8]); % average over all time windows
                                    tf          = permute(     mean(roiSpec.vec.psd  (:,:,:,:,:,:,b3,1),6)   ,[5 6 7 1 2 3 4 8]);
                            end
                            f           = permute(          roiSpec.f        (:,1,:,:,:,:,:         ,:)      ,[5 6 7 1 2 3 4 8]);
                            label = ['psd_' phaseAvList{iii}];
                            cScale = 'log';
                        otherwise
                            dbstack; error('X');
                    end
                    switch phaseAvList{iii}
                        case 'PC'
                            hGram = plot(hTf{r+1},f,tf);
                            % set(hGram([1 2 4]),'Visible','off')
                            % set(hGram(3),'Color','b')
                        case 'nPC'
                            hGram = plot(hTf{r+1},f,tf);
                            % set(hGram([1 2 4]),'Visible','off')
                            % set(hGram(3),'Color','r')
                    end
                    set(hGram,'Color','b')
                    % set(hGram(4),'Color','r')
                end
                roiSpec2 = roiSpec;
            end

            hTf{r+1}.YScale = cScale;
            axis tight;
            if strcmp(roi(r).label,'phys')
                xlim([0 5]);
            end
            % switch metricLabel
            %     case 'psd'
            %         hL = findobj(hTf{r+1}.Children,'type','line');
            %         vis = get(hL,'Visible');
            %         hL = hL([vis{:}]==1);
            %         yLim = ylim(hTf{r+1}); y = yLim(2);
            %         for iiiii = 1:length(hL)
            %             y = min(y,min(hL(iiiii).YData(hL(iiiii).XData>0.1)));
            %         end
            %         yLim(1) = y;
            %         ylim(hTf{r+1},yLim)
            % end
            grid on; grid minor;
            linkprop([hUL{r+1} hTf{r+1}],'Position');
        end

        ind = find([false ~ismember({roi.label},'phys')],1,'last');

        xLim = hTf{ind}.XLim;
        nyq = xLim(2);
        xLim(2) = nyq*7;
        hTf{[false ismember({roi.label},'phys')]}.XLim = xLim;
        xline(hTf{[false ismember({roi.label},'phys')]},nyq,'k')
        ax = hTf{[false ismember({roi.label},'phys')]};
        ax.XTick = hTf{ind}.XTick(1):mean(diff(hTf{ind}.XTick)):ax.XLim(2);
        ax.TickLength = ax.TickLength/3;
        

        yLim = get([hTf{2:ind}],'YLim'); yLim = [min([yLim{:}]) max([yLim{:}])];
        set([hTf{2:ind}],'YLim',yLim);
        posi = get([hUL{2:end}],'Position'); posi = cat(1,posi{:});
        posi = [false; abs(posi(:,1)-min(posi(:,1)))<0.001];
        ylabel([hTf{posi}],metricLabelList{i});
        set([hTf{~posi' & [false ~ismember({roi.label},'phys')]}],'YTickLabel',[]);
        % set(hTf{[false ismember({roi.label},'phys')]},'YAxisLocation','right');
        ylabel([hTf{[false ismember({roi.label},'phys')]}],['physio psd']);

        hLine = hTf{ind}.Children;

        hW = addW(hTf{ind},roiSpec1);
        hW.Color = hLine(1).Color;
        switch metricLabel
            case 'coh'
                hW.YData = mean(hTf{ind}.YLim).*[1 1];
            case 'psd'
                hW.YData = exp(mean(log(hTf{ind}.YLim))).*[1 1];
        end
        hW(2) = addW(hTf{ind},roiSpec2);
        hW(2).Color = hLine(2).Color;
        switch metricLabel
            case 'coh'
                hW(2).YData = mean(hTf{ind}.YLim).*[1.15 1.15];
            case 'psd'
                hW(2).YData = exp(mean(log(hTf{ind}.YLim))).*[1.15 1.15];
        end
        % hW(3) = addW(hTf{ind},roiSpec2);
        % hW(3).Color = hGram(4).Color;
        % switch metricLabel
        %     case 'coh'
        %         hW(3).YData = mean(hTf{end}.YLim).*[1.3 1.3];
        %     case 'psd'
        %         hW(3).YData = exp(mean(log(hTf{end}.YLim))).*[1.3 1.3];
        % end
        
        legend(hTf{ind},hLine,{'late' 'full'},'box','off')
    end
end


hFig = [hF hFrsp hFtf hFspec];
