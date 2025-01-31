function [hOL,hCb] = plotOL(hUL,OL,roi,roiOpt,cLim)
if ~exist('cLim','var'); cLim = [        ]; end
if isempty(cLim)       ; cLim = [-150 150]; end
if ~isnumeric(OL)
    OL = MRIload3(OL,[],[],0);
end
if size(OL,5)>1
    OLalpha = OL(:,:,:,:,2);
    OL(:,:,:,:,2) = [];
else
    OLalpha = [];
end
cMap = [winter; flip(autumn,1)];
% cMap = brewermap([],'RdYlBu');
% cMap = jet;


% Create overlay axes
hOL = addOverlay([hUL{:}]);
for i = 1:length(hOL)
    hOL(i).XAxis.Color = hUL{i}.XAxis.Color;
end
ind0 = 0;

% Plot main map overlays
ind0 = ind0+1;
h = imagesc(hOL(ind0),OL,cLim);
if ~isempty(OLalpha)
    h.AlphaData = OLalpha;
end
hOL(ind0).Colormap = cMap;
set(findobj(hUL{ind0}.Children,'type','polygon'),'Visible','off')
hCb = colorbar(hOL(ind0),'location','manual');
hCb.Position(3) = 0.01;
hCb.Position(1) = 0.99;
hCb.AxisLocation = 'in';

% Plot roi overlays
for rc = 1:length(roi)
    for r = 1:length(roi{rc})
        ind0 = ind0+1;
        if strcmp(roi{rc}(r).label,'vesselAll') || strcmp(roi{rc}(r).label,'phys')
            continue
        end
        h = imagesc(hOL(ind0),roi{rc}(r).im.act.x,roi{rc}(r).im.act.y,roi{rc}(r).im.act.im(:,:,:,:,1),cLim);
        if ~isfield(roi{rc}(r).im,'voxSelP')
            %voxel selection based on fdr corrected p values from this data
            roi{rc}(r).im.voxSelQ       =       roi{rc}(r).im.actP       ; roi{rc}(r).im.voxSelQ.im = ones(size(roi{rc}(r).im.actP.im));
            roi{rc}(r).im.voxSelQ.im(:) = mafdr(roi{rc}(r).im.actP.im(:));            
        else
            %voxel selection based on voxSel fields, which can com from
            %other data
            roi{rc}(r).im.voxSelQ       =       roi{rc}(r).im.voxSelP       ; roi{rc}(r).im.voxSelQ.im = ones(size(roi{rc}(r).im.voxSelP.im));
            roi{rc}(r).im.voxSelQ.im(:) = mafdr(roi{rc}(r).im.voxSelP.im(:),'BHFDR',true);
        end
        h.AlphaData = roi{rc}(r).im.voxSelQ.im<0.05;
        hOL(ind0).Colormap = cMap;
        set(findobj(hUL{ind0}.Children,'type','polygon'),'visible','off')
        % set(findobj(hUL{ind0}.Children,'type','scatter'),'visible','off')
    end
end
set([hOL.XAxis],'Color','none')
set([hOL.YAxis],'Color','none')
