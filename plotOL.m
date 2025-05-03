function [hF,hAO,hIO] = plotOL(rCond,metric,roi,Hbase)
    if ~exist('roi'   ,'var');    roi = struct ; end    
    if ~exist('Hbase' ,'var');  Hbase = []     ; end
    if ~exist('metric','var'); metric = {}     ; end
    if isempty(metric);        metric = {'psd'}; end
    metric = cellstr(metric);

    %% Assert
    if iscell(roi)
        roi = [roi{:}];
    end
    if iscell(Hbase)
        Hbase = [Hbase{:}];
    end
    roi = roi(:);
    Hbase   = Hbase(:);
    % if all(size(roi) == flip(size(H))); roi = roi'; end
    if length(roi) ~= length(Hbase); dbstack; error('roi and H must have the same dimensions'); end

    %% Setup figure
    hF = figure('WindowStyle','docked');
    if isempty(roi)
    else
        hAO = cell(size(Hbase));
        for i = 1:length(Hbase)
            hAO{i} = axes(hF,'Position',Hbase(i).Position,'Box','on');
        end
        hAO = [hAO{:}]; hAO = hAO(:);
    end

    if isfield(roi(1),'im') && (isfield(roi(1).im,'resp') || isfield(roi(1).im,'act'))
        roiDataFlag = true;
    else
        roiDataFlag = false;
    end
    

    %% Plot activation pattern
    if roiDataFlag
        % lineStyle = {'-','-'};
        % lineColor = {[0 0 0],[0.5 0.5 0.5]};
        cLim = cell(1,length(metric));
        hIO = cell(size(roi));
        for i = 1:length(roi)
            hold(hAO(i),'on');
            if length(metric)>1; dbstack; error('accept only one metric, code that'); end
            for m = 1:length(metric)
                switch metric{m}
                    case 'coef'
                        x  = roi(i).im.act.x;
                        y  = roi(i).im.act.y;
                        im = roi(i).im.act.im(:,:,:,1);

                        % figure('WindowStyle','docked');
                        % imagesc(x,y,roi(i).im.act.im(:,:,:,1));
                        % figure('WindowStyle','docked');
                        % imagesc(roi(i).im.base.im(:,:,:,1));
                        
                        % f    = roi(i).vec.mt.psd.f;
                        % spec = mean(roi(i).vec.mt.psd.vec,6);
                    case 'psdPS'
                        f    = roi(i).vec.mt.psdTrialGram.f;
                        spec = mean(roi(i).vec.mt.psdTrialGram.vec(:,:,:,:,:,:,end),6);
                    otherwise
                        error('metric not found');
                end
                cLim{m}(end+1) = max(abs(im(:)));
                hIO{i} = imagesc(hAO(i),x,y,im,[-1 1].*cLim{m}(end));
                set(hAO(i),'YDir','reverse','YTick',[],'XTick',[],'XLim',x+[-0.5 0.5],'YLim',y+[-0.5 0.5],'DataAspectRatio',[1 1 1]);
            end
        end
        drawnow;
        hIO = [hIO{:}]';
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


    %% Add underlay
    drawnow;
    hAU = cell(size(Hbase));
    for i = 1:length(Hbase)
        hAU{i} = axes(hF);
        imU = findobj(Hbase(i).Children,'type','image');
        imagesc(hAU{i},imU.XData,imU.YData,imU.CData,Hbase(i).CLim);
        set(hAU{i},'Position',hAO(i).Position,'Box','on','Colormap',Hbase(i).Colormap,'XTick',[],'YTick',[],'YDir','reverse');
        uistack(hAU{i},'bottom');
    end
    drawnow;

    hAU = [hAU{:}]';

    for i = 1:length(hAU)
        linkaxes([hAU(i) hAO(i)]);
    end
    drawnow;


    %% Set transparency
    for i = 1:length(hAO)
        hIO(i).UserData.roi = roi(i);
        threshOL(hIO(i),'actQ_crop',0);
    end
    

    




    %% Adujst colormap
    maxClrCtrst = 0.75;
    minClrCtrst = 0.3;
    cMap = flip(multigradient(...
    [1 1-maxClrCtrst 1-maxClrCtrst; 1 1-minClrCtrst 1-minClrCtrst; 0.5 0.5 0.5; 1-minClrCtrst 1-minClrCtrst 1; 1-maxClrCtrst 1-maxClrCtrst 1],'pts',...
    [                            0                        0.5-eps          0.5                         0.5+eps                             1]));
    set(hAO,'Colormap',cMap)
    set(hAO,'CLim',[-1 1].*max([cLim{:}]))
    
    
    

    %% Add roi markings
    for i = 1:length(hAO)
        switch roi(i).class
            case 'artery'
                c = 'r';
                c2 = [1 0.7 0.7];
            case 'vein'
                c = 'b';
                c2 = [0.7 0.7 1];
            case 'unknown'
                c = 'y';
                c2 = [0.7 1 1];
                otherwise
                error('Unknown roi class');
        end
        
        % add roi contour
        hOtln(i,1) = plot(hAO(i),roi(i).poly(1),'FaceColor','none','EdgeColor',c);

        % axes outline
        hAO(i).XAxis.Color = c;     hAO(i).YAxis.Color     = hAO(i).XAxis.Color;
        hAO(i).XAxis.LineWidth = 2; hAO(i).YAxis.LineWidth = hAO(i).XAxis.LineWidth;

        % add roi center of mass
        hCom(i,1) = xline(hAO(i),roi(i).com(1),'w');
        hCom(i,2) = yline(hAO(i),roi(i).com(2),'w');

        % add roi id
        hId(i,1) = text(hAO(i),hAO(i).XLim(1),hAO(i).YLim(1),num2str(roi(i).id),'Color',c2,'FontSize',16,'HorizontalAlignment','left','VerticalAlignment','top');% text(tilingSubExtra.hA,tilingSub.hA{end}{end}.XLim(1),tilingSub.hA{end}{end}.YLim(1),num2str(roi(rc).id),'Color','w','HorizontalAlignment','left','VerticalAlignment','top');


    end



    
    %% Adjust axes
    for i = 1:length(hAU)
        switch roi(i).class
            case 'artery'
                c = 'r';
            case 'vein'
                c = 'b';
            case 'unknown'
                c = 'y';
            otherwise
                error('Unknown roi class');
        end
        % hAU(i).Position
        hAO(i).XAxis.Color = c;     hAO(i).YAxis.Color = hAO(i).XAxis.Color;
        hAO(i).XAxis.LineWidth = 2; hAO(i).YAxis.LineWidth = hAO(i).XAxis.LineWidth;
    end


    
    % for i = 1:length(Hbase)
    %     hAO(i).XAxis.Color = Hbase(i).XAxis.Color; hAO(i).XAxis.LineWidth = Hbase(i).XAxis.LineWidth;
    %     hAO(i).YAxis.Color = Hbase(i).YAxis.Color; hAO(i).YAxis.LineWidth = Hbase(i).YAxis.LineWidth;
    %     hAO(i).XLim = hAO(i).XLim;
    %     hAO(i).YLim = hAO(i).YLim;

    %     hAO(i).Position = hAU(i).Position;
    %     linkaxes([hAO(i) hAU(i)]);
    % end
    
    


    % figure('WindowStyle','docked');
    % hx = plot(cMap);
    % hx(1).Color = [1 0 0];
    % hx(2).Color = [0 1 0];
    % hx(3).Color = [0 0 1];


    
    % return


    % get([hAU{:}],'CLim')
    % set([hAU{:}],'CLim',[250 900]);
    
    % set([hAU{:}],'Visible','on');
    % for i = 1:length(Hbase)
    %     hAU{i}.Visible = 'off';
    % end

    % delete([hAU{:}]);


    % set(hAO,'XLimMode','manual');
    % set(hAO,'YLimMode','manual');
    % set([hIO{:}],'AlphaData',0);
    