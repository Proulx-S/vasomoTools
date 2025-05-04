function roi = smrRoi(rCond,metric,roi,H)
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
                switch metric{m}
                    case 'psd_dilate1_actQ'
                        % PSD within the vessel ROI dilated by 1 voxel,
                        % including only active voxels based on SPMG2 activation detection
                        indIm         = roi(i).polyMask{ismember(roi(i).polyLabel,'dilate1')};
                        indSig        = false(size(roi(i).im.actP.im));
                        indSig(indIm) = mafdr(roi(i).im.actP.im(indIm),'BHFDR',true)<0.05;
                        vec = mean(roi(i).mt.psd.vec(:,:,:,:,:,indIm&indSig,:,:),6);
                        f   = roi(i).mt.psd.f;
                        
                        roi(i).smr{m}.vec     = permute(vec,[5 6 1 2 3 4 7 8]);
                        roi(i).smr{m}.nVox    = [];
                        roi(i).smr{m}.nVoxRoi = [nnz(indIm)];
                        roi(i).smr{m}.x       = permute(f,[5 6 1 2 3 4 7 8]);
                        roi(i).smr{m}.label   = [];
                        roi(i).smr{m}.metric  = metric{m};
                        
                        
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
                        
                        roi(i).smr{m}.vec     = [tsNeg tsPos];
                        roi(i).smr{m}.nVox    = [nnz(indIm&indSig&indNeg) nnz(indIm&indSig&indPos)];
                        roi(i).smr{m}.nVoxRoi = [nnz(indIm)];
                        roi(i).smr{m}.x       = t;
                        roi(i).smr{m}.label   = {'neg','pos'};
                        roi(i).smr{m}.metric  = metric{m};

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
        for m = 1:length(metric)
            %% Adjust axes
            for i = 1:length(H)
                hA{m}{i}.XAxis.Color = H(i).XAxis.Color; hA{m}{i}.XAxis.LineWidth = H(i).XAxis.LineWidth;
                hA{m}{i}.YAxis.Color = H(i).YAxis.Color; hA{m}{i}.YAxis.LineWidth = H(i).YAxis.LineWidth;
            end
            axis(hA{m},'tight'); yLim = get(hA{m},'YLim'); yLim = [min([yLim{:}]) max([yLim{:}])];
            set(hA{m},...
            'YLim',yLim,...
            'YScale','linear',...
            'XGrid','on','YGrid','on',...
            'XMinorGrid','on','YMinorGrid','on',...
            'GridColor',[0.5 0.5 0.5],'MinorGridColor',[0.5 0.5 0.5]);
            drawnow;
        end

        for m = 1:length(metric)
            for i = 1:length(roi)
                % Add text annotations for voxel counts
                % Bottom left corner - total ROI voxel count
                text(hA{m}{i}, min(hA{m}{i}.XLim)+range(hA{m}{i}.XLim)*0.01, min(hA{m}{i}.YLim)+range(hA{m}{i}.YLim)*0.01, ...
                    [num2str(roi(i).ts{m}.nVoxRoi) 'vox'], ...
                    'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom', ...
                    'FontSize', 8);
                % Top right corner - positive and significant voxel count
                text(hA{m}{i},...
                min(hA{m}{i}.XLim)+range(hA{m}{i}.XLim)*0.99, min(hA{m}{i}.YLim)+range(hA{m}{i}.YLim)*0.99, ...
                [num2str(roi(i).ts{m}.nVox(ismember(roi(i).ts{m}.label,'pos'))) 'posVox'], ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'top', ...
                'FontSize', 8);
            % Bottom right corner - negative and significant voxel count
                text(hA{m}{i}, min(hA{m}{i}.XLim)+range(hA{m}{i}.XLim)*0.99, min(hA{m}{i}.YLim)+range(hA{m}{i}.YLim)*0.01, ...
                    [num2str(roi(i).ts{m}.nVox(ismember(roi(i).ts{m}.label,'neg'))) 'negVox'], ...
                    'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 8);
            end
        end

    end
    

