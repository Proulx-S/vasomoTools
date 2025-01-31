function plotAct(volAct)

hF = figure('WindowStyle','docked');
hT = tiledlayout(16,38); hT.Padding = 'tight'; hT.TileSpacing = 'tight';
ax = {};
cLimUL = [0 800];
cLimOL = [-150 150];

% plot underlay
imBase = MRIread(volAct.fs.fBaseAv); imBase = imBase.vol(:,:,:,end);
ax{end+1} = nexttile([16 16]);
imagesc(ax{end},imBase,cLimUL); colormap(ax{end},'gray'); ax{end}.DataAspectRatio = [1 1 1];
% ax{end}.XAxis.Visible = 'off'; ax{end}.YAxis.Visible = 'off';
ax{end}.XAxis.Color = 'w'; ax{end}.XTick = [];
ax{end}.YAxis.Color = 'w'; ax{end}.YTick = [];
hold on

fCoef = runCond{S}.(runCondAcq).(runCondStim).volActCat.fs.fCoef;
imAct = getSPMG2coef(fCoef,fVes);
imQ = MRIload3(runCond{S}.(runCondAcq).(runCondStim).volActCat.fs.fFullQ,[],[],0); imQ = imQ.vol;
imP = MRIload3(runCond{S}.(runCondAcq).(runCondStim).volActCat.fs.fFullP,[],[],0); imP = imP.vol;

% plot individual vessels
fVes  = volAnat.calcarineVessel.f;
imVes = MRIread(volAnat.calcarineVessel.f); imVes = imVes.vol;
%%% extract vessel roi
cropSize = 10;
artRoi = getVesselRoi(imVes==902         ,{'base' 'act' 'actP'},{imBase imAct imP},cropSize);
veiRoi = getVesselRoi(imVes==914         ,{'base' 'act' 'actP'},{imBase imAct imP},cropSize);
ambRoi = getVesselRoi(imVes==30|imVes==62,{'base' 'act' 'actP'},{imBase imAct imP},cropSize);
%%% plot vessel underlay
roi = artRoi;
for i = 1:length(roi)
    ax{end+1} = nexttile(tilenum(hT,1,17+3*(i-1)),[3 3]);
    imagesc(roi(i).im.base.x,roi(i).im.base.y,roi(i).im.base.im,cLim); colormap(ax{end},'gray'); ax{end}.DataAspectRatio = [1 1 1];
    % ax{end}.XAxis.Visible = 'off'; ax{end}.YAxis.Visible = 'off';
    ax{end}.XAxis.Color = 'r'; ax{end}.XTick = [];
    ax{end}.YAxis.Color = 'r'; ax{end}.YTick = [];

    hold on
    scatter(roi(i).com(1),roi(i).com(2),'.r')
    plot(roi(i).poly,'FaceColor','none','EdgeColor','r')
    plot(ax{1},roi(i).poly,'FaceColor','none','EdgeColor','r')
end
roi = veiRoi;
for i = 1:length(roi)
    ax{end+1} = nexttile(tilenum(hT,4,17+3*(i-1)),[3 3]);
    imagesc(roi(i).im.base.x,roi(i).im.base.y,roi(i).im.base.im,cLim); colormap(ax{end},'gray'); ax{end}.DataAspectRatio = [1 1 1];
    % ax{end}.XAxis.Visible = 'off'; ax{end}.YAxis.Visible = 'off';
    ax{end}.XAxis.Color = 'b'; ax{end}.XTick = [];
    ax{end}.YAxis.Color = 'b'; ax{end}.YTick = [];

    hold on
    scatter(roi(i).com(1),roi(i).com(2),'.b')
    plot(roi(i).poly,'FaceColor','none','EdgeColor','b')
    plot(ax{1},roi(i).poly,'FaceColor','none','EdgeColor','b')
end
roi = ambRoi;
for i = 1:length(roi)
    ax{end+1} = nexttile(tilenum(hT,7,17+3*(i-1)),[3 3]);
    imagesc(roi(i).im.base.x,roi(i).im.base.y,roi(i).im.base.im,cLim); colormap(ax{end},'gray'); ax{end}.DataAspectRatio = [1 1 1];
    % ax{end}.XAxis.Visible = 'off'; ax{end}.YAxis.Visible = 'off';
    ax{end}.XAxis.Color = 'y'; ax{end}.XTick = [];
    ax{end}.YAxis.Color = 'y'; ax{end}.YTick = [];

    hold on
    scatter(roi(i).com(1),roi(i).com(2),'.y')
    plot(roi(i).poly,'FaceColor','none','EdgeColor','y')
    plot(ax{1},roi(i).poly,'FaceColor','none','EdgeColor','y')
end

% Plot activation map overlay
axOl = addOverlay([ax{:}]);
ind0 = 0;

fCoef = runCond{S}.(runCondAcq).(runCondStim).volActCat.fs.fCoef;
imAct = getSPMG2coef(fCoef,fVes);

ind0 = ind0+1;
h = imagesc(axOl(ind0),imAct,cLimOL)
h.AlphaData = imQ<0.05
axOl(ind0).Colormap = jet;
set(findobj(ax{ind0}.Children,'type','polygon'),'Visible','off')

roi = artRoi;
for i = 1:length(roi)
    ind0 = ind0+1;
    h = imagesc(axOl(ind0),roi(i).im.act.x,roi(i).im.act.y,roi(i).im.act.im,cLimOL);
    roi(i).im.actQ = roi(i).im.actP; roi(i).im.actQ.im = ones(size(roi(i).im.actQ.im));
    roi(i).im.actQ.im(:) = mafdr(roi(i).im.actP.im(:));
    h.AlphaData = roi(i).im.actQ.im<0.05;
    axOl(ind0).Colormap = jet;
    set(findobj(ax{ind0}.Children,'type','polygon'),'visible','off')
    set(findobj(ax{ind0}.Children,'type','scatter'),'visible','off')
end
roi = veiRoi;
for i = 1:length(roi)
    ind0 = ind0+1;
    h = imagesc(axOl(ind0),roi(i).im.act.x,roi(i).im.act.y,roi(i).im.act.im,cLimOL);
    roi(i).im.actQ = roi(i).im.actP; roi(i).im.actQ.im = ones(size(roi(i).im.actQ.im));
    roi(i).im.actQ.im(:) = mafdr(roi(i).im.actP.im(:));
    h.AlphaData = roi(i).im.actQ.im<0.05;
    axOl(ind0).Colormap = jet;
    set(findobj(ax{ind0}.Children,'type','polygon'),'visible','off')
    set(findobj(ax{ind0}.Children,'type','scatter'),'visible','off')
end
roi = ambRoi;
for i = 1:length(roi)
    ind0 = ind0+1;
    h = imagesc(axOl(ind0),roi(i).im.act.x,roi(i).im.act.y,roi(i).im.act.im,cLimOL);
    roi(i).im.actQ = roi(i).im.actP; roi(i).im.actQ.im = ones(size(roi(i).im.actQ.im));
    roi(i).im.actQ.im(:) = mafdr(roi(i).im.actP.im(:));
    h.AlphaData = roi(i).im.actQ.im<0.05;
    axOl(ind0).Colormap = jet;
    set(findobj(ax{ind0}.Children,'type','polygon'),'visible','off')
    set(findobj(ax{ind0}.Children,'type','scatter'),'visible','off')
end


