function tiling = plotUL3(roi,UL,cLim,numRow)
    if isstruct(roi); roi = {roi}; end
    if ~exist('UL','var')  ;       UL = [     ]  ; end
    if ~exist('cLim','var');     cLim = [     ]  ; end
    if ~exist('numRow','var'); numRow = [     ]  ; end
    if isempty(cLim)       ;     cLim = [100 800]; end
    if isempty(UL) % extract full-size underlay image from roi
        UL = [roi{:}]; UL = [UL.im]; UL = [UL.base]; UL = unique({UL.fName}); if length(UL) ~= 1; dbstack; error('UL should be a single file name'); end; UL = char(UL);
    end
    if isempty(numRow); numRow = 4; end;

    %% Setup tiling
    tiling.main.gridSize      = [16 38];
    tiling.sub.left.gridSize  = [16 16];
    tiling.sub.left.axesSize  = [16 16];
    tiling.sub.right.gridSize = [16 22];
    tiling.sub.right.axesSize = [ 1  1].*floor(tiling.sub.right.gridSize(1)./numRow);
    tiling.sub.right.row0 = 1;
    tiling.sub.right.col0 = 1 + tiling.sub.left.gridSize(2);

    tiling.main.hF = figure('WindowStyle','docked');
    tiling.main.hT = tiledlayout(tiling.main.gridSize(1),tiling.main.gridSize(2)); tiling.main.hT.Padding = 'tight'; tiling.main.hT.TileSpacing = 'tight';


    %% Plot full-size underlay image on the left size of the figure
    if ~isnumeric(UL)
        UL = MRIread(UL); UL = UL.vol;
    end
    tiling.sub.left.hA = nexttile(tiling.sub.left.axesSize);
    imagesc(tiling.sub.left.hA,UL,cLim);
    colormap(tiling.sub.left.hA,'gray'); tiling.sub.left.hA.DataAspectRatio = [1 1 1];
    tiling.sub.left.hA.XAxis.Color = 'w'; tiling.sub.left.hA.XTick = []; tiling.sub.left.hA.XAxis.LineWidth = 1;
    tiling.sub.left.hA.YAxis.Color = 'w'; tiling.sub.left.hA.YTick = []; tiling.sub.left.hA.YAxis.LineWidth = 1;
    hold on

    %% Plot rois
    if isstruct(roi)
        roi = {roi};
    elseif iscell(roi)
    else
        error('Unknown roi type');
    end

    % set starting row
    tiling.sub.right.row = tiling.sub.right.row0;
    for rc = 1:length(roi)
        % reset starting column
        tiling.sub.right.col = tiling.sub.right.col0;
        % plot roi row
        [tiling.sub.right,tiling.sub.left] = plotRoiRow(tiling.main,tiling.sub.right,tiling.sub.left,roi{rc},cLim);
        % update row
        tiling.sub.right.row = tiling.sub.right.row+tiling.sub.right.axesSize(1);
    end

    drawnow;





function [tilingSub,tilingSubExtra] = plotRoiRow(tilingMain,tilingSub,tilingSubExtra,roi,cLim)

    % intiate axes handles
    if ~isfield(tilingSub,'hA');
        tilingSub.hA = {{}};
    else
        tilingSub.hA{end+1,1} = {};
    end

    % loop through rois
    for rc = 1:size(roi,1)
        switch roi(rc).class
            case 'phys'
                dbstack; error('code that')
                tiling.hA{end+1} = nexttile([3 3*7]);
            otherwise
                % update row if reaching last column
                if tilingSub.col + tilingSub.axesSize(2) > tilingMain.gridSize(2)
                    tilingSub.row = tilingSub.row+tilingSub.axesSize(1);
                    tilingSub.col = tilingSub.col0;
                end
                % create new axe
                tilingSub.hA{end}{1,end+1} = nexttile(tilenum(tilingMain.hT,tilingSub.row,tilingSub.col),tilingSub.axesSize);
        end

        % plot base image
        imagesc(roi(rc).im.base.x,roi(rc).im.base.y,roi(rc).im.base.im,cLim);
        colormap(tilingSub.hA{end}{end},'gray'); hold on
        switch roi(rc).class
            case 'phys'
                dbstack; error('code that')
                tiling.hA{end}.PlotBoxAspectRatio = [7 1 1];
            otherwise
                tilingSub.hA{end}{end}.DataAspectRatio = [1 1 1];
        end

        % roi class specific colors
        switch roi(rc).class
            case 'artery'
                c = 'r';
            case 'vein'
                c = 'b';
            case 'unknown'
                c = [0.5 0.5 0];
            case 'phys'
                c = 'k';
            otherwise
                error('Unknown roi class');
        end
        tilingSub.hA{end}{end}.XAxis.Color = c; tilingSub.hA{end}{end}.XTick = []; tilingSub.hA{end}{end}.XAxis.LineWidth = 1;
        tilingSub.hA{end}{end}.YAxis.Color = c; tilingSub.hA{end}{end}.YTick = []; tilingSub.hA{end}{end}.YAxis.LineWidth = 1;

        switch roi(rc).class
            case {'artery','vein','unknown'}
                % add roi contours
                plot(tilingSub.hA{end}{end},roi(rc).poly(1),'FaceColor','none','EdgeColor',c);
                if ~isempty(tilingSubExtra)
                    plot(tilingSubExtra.hA,roi(rc).poly(1),'FaceColor','none','EdgeColor',c);
                end
                if ~iscell(roi(rc).com)
                    % add roi center of mass
                    xline(tilingSub.hA{end}{end},roi(rc).com(1),'w');
                    yline(tilingSub.hA{end}{end},roi(rc).com(2),'w');
                    % add roi id
                    text(tilingSub.hA{end}{end},tilingSub.hA{end}{end}.XLim(1),tilingSub.hA{end}{end}.YLim(1),num2str(roi(rc).id),'Color','w','HorizontalAlignment','left','VerticalAlignment','top');
                    text(tilingSubExtra.hA,tilingSub.hA{end}{end}.XLim(1),tilingSub.hA{end}{end}.YLim(1),num2str(roi(rc).id),'Color','w','HorizontalAlignment','left','VerticalAlignment','top');
                end
            case 'phys'
            otherwise
                error('Unknown roi class');
        end

        % update column
        tilingSub.col = tilingSub.col+tilingSub.axesSize(2);

        % add roi info directly to graphics object
        tilingSub.hA{end}{end}.UserData = roi(rc);
    end

    % concatenate axes handles
    tilingSub.hA{end} = cat(2,tilingSub.hA{end}{:});

    drawnow;