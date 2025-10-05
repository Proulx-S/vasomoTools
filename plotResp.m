function [roi,hF,hA,hTs] = plotResp(rCond,metric,roi,H)
    if ~exist('roi','var');           roi = struct; end    
    if ~exist('H','var');               H = []    ; end
    if ~exist('metric','var');     metric = {}; end
    if isempty(metric);            metric = {'resp'}; end

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
        if usejava('desktop')
            hF = figure('WindowStyle','docked');
        else
            hF = figure('MenuBar','none','ToolBar','none');
        end
        
        if isempty(roi)
        else
            hA = cell(size(H));
            for i = 1:length(H)
                hA{i} = axes(hF,'Position',H(i).Position,'Box','on');
            end
            hA = [hA{:}]; hA = hA(:);
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
            if ~isempty(H)
                hold(hA(i),'on');
            end
            for m = 1:length(metric)
                metric1 = strsplit(metric{m},'_');
                if length(metric1)>1; metric2 = metric1{2}; else metric2 = ''; end; metric1 = metric1{1};
                switch metric1
                    case {'respSurVox' 'respPeakVox'}
                        % response within the surrounding voxel ROI
                        ts = mean(permute(roi(i).im.(metric{m}).vec,[2 1]),2);
                        t  = (0:size(ts,1)-1).*roi(i).im.resp.dt;
                        nVox = size(roi(i).im.(metric{m}).vec,1);
                        
                        roi(i).ts{m}.vec     = ts;
                        roi(i).ts{m}.nVox    = nVox;
                        roi(i).ts{m}.nVoxRoi = nVox;
                        roi(i).ts{m}.t       = t;
                        roi(i).ts{m}.label   = {};
                        roi(i).ts{m}.metric  = metric{m};

                        if ~isempty(H)
                            hTs{i,m} = plot(hA(i),roi(i).ts{m}.t,roi(i).ts{m}.vec,'Color',lineColor{1},'LineStyle',lineStyle{1});
                        end
                    case 'svTime'
                        % temporal singular vector from SVD transform
                        comp = 1;
                        if ~isempty(metric2) && all(ismember(metric2,'0123456789'))
                            comp = str2double(metric2);
                        end
                        ts = roi(i).svdResp.svTime(comp,:)';
                        t  = (0:size(ts,1)-1).*roi(i).im.resp.dt;
                        
                        roi(i).ts{m}.vec     = ts;
                        roi(i).ts{m}.nVox    = nnz(roi(i).svdResp.maskSVD);
                        roi(i).ts{m}.nVoxRoi = nnz(roi(i).svdResp.maskSVD);
                        roi(i).ts{m}.t       = t;
                        roi(i).ts{m}.label   = {};
                        roi(i).ts{m}.metric  = metric{m};

                        if ~isempty(H)
                            hTs{i,m} = plot(hA(i),roi(i).ts{m}.t,roi(i).ts{m}.vec,'Color',lineColor{2},'LineStyle',lineStyle{2});
                        end
                    case 'respArea'
                        % area change response
                        ts = permute(roi(i).im.respArea.vec,[2 1]);
                        t  = (0:size(ts,1)-1).*roi(i).im.resp.dt;
                        nVox = nnz(roi(i).im.respArea.maskResp.wMask|roi(i).im.respArea.maskResp.zMask);
                        
                        roi(i).ts{m}.vec     = ts;
                        roi(i).ts{m}.nVox    = nVox;
                        roi(i).ts{m}.nVoxRoi = nVox;
                        roi(i).ts{m}.t       = t;
                        roi(i).ts{m}.label   = {};
                        roi(i).ts{m}.metric  = metric{m};

                        if ~isempty(H)
                            hTs{i,m} = plot(hA(i),roi(i).ts{m}.t,roi(i).ts{m}.vec,'Color',lineColor{1},'LineStyle',lineStyle{1});
                        end
                    case 'respVel'
                        % velocity change response (just the peak voxel for now)
                        ts = permute(roi(i).im.respVel.vec,[2 1 3 4]);
                        t  = (0:size(ts,1)-1).*roi(i).im.resp.dt;
                        nVox = nnz(roi(i).im.respVel.im2vec);

                        roi(i).ts{m}.vec     = ts;
                        roi(i).ts{m}.nVox    = nVox;
                        roi(i).ts{m}.nVoxRoi = nVox;
                        roi(i).ts{m}.t       = t;
                        roi(i).ts{m}.label   = {};
                        roi(i).ts{m}.metric  = metric{m};

                        if ~isempty(H)
                            hTs{i,m} = plot(hA(i),roi(i).ts{m}.t,roi(i).ts{m}.vec,'Color',lineColor{1},'LineStyle',lineStyle{1});
                        end

                    case 'resp_dilate1_actQ_actSgn'
                        % response within the vessel ROI dilated by 1 voxel,
                        % including only active voxels based on SPMG2 activation detection,
                        % and segregated by sign of activation
                        im = permute(roi(i).im.resp.im,[4 1 2 3]);
                        indIm         = roi(i).polyMask{ismember(roi(i).polyLabel,'dilate1')};
                        indSig        = false(size(roi(i).im.actP.im));
                        indSig(indIm) = mafdr(roi(i).im.actP.im(indIm),'BHFDR',true)<0.05;
                        indNeg        = roi(i).im.act.im(:,:,:,1)<0;
                        indPos        = roi(i).im.act.im(:,:,:,1)>0;
                        
                        tsNeg = mean(im(:,indIm&indSig&indNeg),2);
                        tsPos = mean(im(:,indIm&indSig&indPos),2);
                        t     = (0:size(tsNeg,1)-1).*roi(i).im.resp.dt;
                        
                        roi(i).ts{m}.vec     = [tsNeg tsPos];
                        roi(i).ts{m}.nVox    = [nnz(indIm&indSig&indNeg) nnz(indIm&indSig&indPos)];
                        roi(i).ts{m}.nVoxRoi = [nnz(indIm)];
                        roi(i).ts{m}.t       = t;
                        roi(i).ts{m}.label   = {'neg','pos'};
                        roi(i).ts{m}.metric  = metric{m};


                        if ~isempty(H)
                            hTs{i,m} = plot(hA(i),roi(i).ts{m}.t,roi(i).ts{m}.vec);
                        end
                        
                    otherwise
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
        %% Adjust axes
        for i = 1:length(H)
            hA(i).XAxis.Color = H(i).XAxis.Color; hA(i).XAxis.LineWidth = H(i).XAxis.LineWidth;
            hA(i).YAxis.Color = H(i).YAxis.Color; hA(i).YAxis.LineWidth = H(i).YAxis.LineWidth;
        end
        axis(hA,'tight'); yLim = get(hA,'YLim'); yLim = [min([yLim{:}]) max([yLim{:}])];
        set(hA,...
        'YLim',yLim,...
        'YScale','linear',...
        'XGrid','on','YGrid','on',...
        'XMinorGrid','on','YMinorGrid','on',...
        'GridColor',[0.5 0.5 0.5],'MinorGridColor',[0.5 0.5 0.5]);
        drawnow;

        m = length(metric);
        if ~strcmp(metric{m},'respSurVox') && ~strcmp(metric{m},'respPeakVox')
            for i = 1:length(roi)
                % Add text annotations for voxel counts
                % Bottom left corner - total ROI voxel count
                text(hA(i), min(hA(i).XLim)+range(hA(i).XLim)*0.01, min(hA(i).YLim)+range(hA(i).YLim)*0.01, ...
                    [num2str(roi(i).ts{m}.nVoxRoi) 'vox'], ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 8);
                % Top right corner - positive and significant voxel count
                text(hA(i), min(hA(i).XLim)+range(hA(i).XLim)*0.99, min(hA(i).YLim)+range(hA(i).YLim)*0.99, ...
                    [num2str(roi(i).ts{m}.nVox(ismember(roi(i).ts{m}.label,'pos'))) 'posVox'], ...
                    'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', ...
                    'FontSize', 8);
                % Bottom right corner - negative and significant voxel count
                text(hA(i), min(hA(i).XLim)+range(hA(i).XLim)*0.99, min(hA(i).YLim)+range(hA(i).YLim)*0.01, ...
                    [num2str(roi(i).ts{m}.nVox(ismember(roi(i).ts{m}.label,'neg'))) 'negVox'], ...
                    'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 8);
            end
        end


        %% Adjust line color based on desktop theme
        drawnow;
        for i = 1:length(roi)
            hLine = findobj(hA(i).Children,'Type','Line');
            if isempty(hLine); continue; end
            if all(hA(i).Color ~= [1 1 1]) && length(hLine)==1
                hLine.Color = 'w';
            end
        end

    end
    

