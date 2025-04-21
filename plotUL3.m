function [hF,hT,hA] = plotUL3(UL,roi,cLim)

if ~exist('cLim','var'); cLim = [     ]; end
if isempty(cLim)       ; cLim = [0 800]; end
if ~isnumeric(UL)
    UL = MRIload3(UL,[],[],0); UL = UL.vol(:,:,:,end);
end



hF = figure('WindowStyle','docked');
hT = tiledlayout(16,38); hT.Padding = 'tight'; hT.TileSpacing = 'tight';
% hT.TileIndexing = 'columnmajor';
hA = {};


%% Plot full image
hA{end+1} = nexttile([16 16]);
imagesc(hA{end},UL,cLim); colormap(hA{end},'gray'); hA{end}.DataAspectRatio = [1 1 1];
% ax{end}.XAxis.Visible = 'off'; ax{end}.YAxis.Visible = 'off';
hA{end}.XAxis.Color = 'w'; hA{end}.XTick = [];
hA{end}.YAxis.Color = 'w'; hA{end}.YTick = [];
hold on

%% Plot rois
row = 1-3;
for rc = 1:length(roi)
    % row  = row+3;
    % col0 = 17;
    % col = col0;
    

    if 1
        switch roi(rc).class
            case 'phys'
                hA{end+1} = nexttile([3 3*7]);
            otherwise
                hA{end+1} = nexttile([3 3]);
        end

        imagesc(roi(rc).im.base.x,roi(rc).im.base.y,roi(rc).im.base.im,cLim);
        colormap(hA{end},'gray'); hold on

        switch roi(rc).class
            case 'phys'
                hA{end}.PlotBoxAspectRatio = [7 1 1];
            otherwise
                hA{end}.DataAspectRatio = [1 1 1];
        end

        switch roi(rc).class
            case 'artery'
                hA{end}.XAxis.Color = 'r'; hA{end}.XTick = [];
                hA{end}.YAxis.Color = 'r'; hA{end}.YTick = [];
                plot(hA{end},roi(rc).poly,'FaceColor','none','EdgeColor','r');
                plot(hA{1  },roi(rc).poly,'FaceColor','none','EdgeColor','r');
            case 'vein'
                hA{end}.XAxis.Color = 'b'; hA{end}.XTick = [];
                hA{end}.YAxis.Color = 'b'; hA{end}.YTick = [];
                plot(hA{end},roi(rc).poly,'FaceColor','none','EdgeColor','b');
                plot(hA{1  },roi(rc).poly,'FaceColor','none','EdgeColor','b');
            case 'unknown'
                hA{end}.XAxis.Color = [0.5 0.5 0]; hA{end}.XTick = [];
                hA{end}.YAxis.Color = [0.5 0.5 0]; hA{end}.YTick = [];
                plot(hA{end},roi(rc).poly,'FaceColor','none','EdgeColor',[0.5 0.5 0]);
                plot(hA{1  },roi(rc).poly,'FaceColor','none','EdgeColor',[0.5 0.5 0]);
            case 'phys'
                hA{end}.XAxis.Color = 'k'       ; hA{end}.XTick = [];
                hA{end}.YAxis.Color = 'k'       ; hA{end}.YTick = [];
            otherwise
                error('Unknown roi class');
        end


        % if ~strcmp(roi{rc}(r).label,'vesselAll') && ~skipThis
            xline(roi(rc).com(1),'w');
            yline(roi(rc).com(2),'w');
            hA{end}.XAxis.LineWidth = 1;
            hA{end}.YAxis.LineWidth = 1;
        % else
        %     hA{end}.XAxis.LineWidth = 3;
        %     hA{end}.YAxis.LineWidth = 3;
        % end
        text(hA{end},hA{end}.XLim(1),hA{end}.YLim(1),num2str(roi(rc).id),'Color','w','HorizontalAlignment','left','VerticalAlignment','top');
        text(hA{1},hA{end}.XLim(1),hA{end}.YLim(1),num2str(roi(rc).id),'Color','w','HorizontalAlignment','left','VerticalAlignment','top');

    else


        for r = 1:length(roi(rc))
            % if col>38-3; row = row+3; col = col0; end
            skipThis = ~isfield(roi{rc},'im') || isempty(roi{rc}(r).im);
            % hA{end+1} = nexttile(tilenum(hT,row,col),[3 3]);
            if strcmp(roi{rc}(r).label,'phys')
                hA{end+1} = nexttile([3 3*7]);
            else
                hA{end+1} = nexttile([3 3]);
            end
            if ~skipThis
                imagesc(roi{rc}(r).im.base.x,roi{rc}(r).im.base.y,roi{rc}(r).im.base.im,cLim); colormap(hA{end},'gray');
            end
            if strcmp(roi{rc}(r).label,'phys')
                hA{end}.PlotBoxAspectRatio = [7 1 1];
            else
                hA{end}.DataAspectRatio = [1 1 1];
            end
            if isempty(roiOpt{rc})
                hA{end}.XAxis.Color = 'k'       ; hA{end}.XTick = [];
                hA{end}.YAxis.Color = 'k'       ; hA{end}.YTick = [];
            else
                if strcmp(roiOpt{rc},'y')
                    hA{end}.XAxis.Color = [0.5 0.5 0]; hA{end}.XTick = [];
                    hA{end}.YAxis.Color = [0.5 0.5 0]; hA{end}.YTick = [];
                else
                    hA{end}.XAxis.Color = roiOpt{rc}; hA{end}.XTick = [];
                    hA{end}.YAxis.Color = roiOpt{rc}; hA{end}.YTick = [];
                end
            end

            hold on
            if ~strcmp(roi{rc}(r).label,'vesselAll') && ~skipThis
                xline(roi{rc}(r).com(1),'w');
                yline(roi{rc}(r).com(2),'w');
                hA{end}.XAxis.LineWidth = 1;
                hA{end}.YAxis.LineWidth = 1;
            else
                hA{end}.XAxis.LineWidth = 3;
                hA{end}.YAxis.LineWidth = 3;
            end

            if skipThis
                hA{end}.Visible = 'off';
            else
                if isempty(roiOpt{rc})
                    plot(roi{rc}(r).poly,'FaceColor','none','EdgeColor','k'       )
                else
                    plot(      roi{rc}(r).poly,'FaceColor','none','EdgeColor',roiOpt{rc})
                    plot(hA{1},roi{rc}(r).poly,'FaceColor','none','EdgeColor',roiOpt{rc})
                end
            end
            % col = col+3;
        end
    end
end