function plotSpec(rCond,metric,roi,H)
    if ~exist('roi','var');       roi = struct; end    
    if ~exist('H','var');           H = []    ; end
    if ~exist('metric','var'); metric = {}; end
    if isempty(metric);        metric = {'psd'}; end
    metric = cellstr(metric);

    %% Assert
    if iscell(roi)
        roi = [roi{:}];
    end
    if iscell(H)
        H = [H{:}];
    end
    roi = roi(:);
    H   = H(:);
    % if all(size(roi) == flip(size(H))); roi = roi'; end
    if length(roi) ~= length(H); dbstack; error('roi and H must have the same dimensions'); end

    %% Setup figure
    hF = figure('WindowStyle','docked');
    if isempty(roi)
    else
        hA = cell(size(H));
        for i = 1:length(H)
            hA{i} = axes(hF,'Position',H(i).Position,'Box','on');
        end
        hA = [hA{:}]; hA = hA(:);
    end

    if isfield(roi(1),'vec') && isfield(roi(1).vec,'mt')
        roiDataFlag = true;
    else
        roiDataFlag = false;
    end
    

    %% Plot spectra
    if roiDataFlag
        lineStyle = {'-','-'};
        lineColor = {[0 0 0],[0.5 0.5 0.5]};
        for i = 1:length(roi)
            hold(hA(i),'on');
            for m = 1:length(metric)
                switch metric{m}
                    case 'psd'
                        f    = roi(i).vec.mt.psd.f;
                        spec = mean(roi(i).vec.mt.psd.vec,6);
                    case 'psdPS'
                        f    = roi(i).vec.mt.psdTrialGram.f;
                        spec = mean(roi(i).vec.mt.psdTrialGram.vec(:,:,:,:,:,:,end),6);
                        otherwise
                end
                plot(hA(i),squeeze(f),squeeze(spec),lineStyle{m},'Color',lineColor{m});
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



    %% Adjust axes
    for i = 1:length(H)
        hA(i).XAxis.Color = H(i).XAxis.Color; hA(i).XAxis.LineWidth = H(i).XAxis.LineWidth;
        hA(i).YAxis.Color = H(i).YAxis.Color; hA(i).YAxis.LineWidth = H(i).YAxis.LineWidth;
    end
    axis(hA,'tight'); yLim = get(hA,'YLim'); yLim = [min([yLim{:}]) max([yLim{:}])];
    set(hA,'YLim',yLim,'YScale','log','XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on');

