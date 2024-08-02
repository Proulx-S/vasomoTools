function tryL2svd(volPsd)
method = 1;
adjustPhase = 1;
saveFlag = 1;

N = volPsd.svd.dim(1);
W = volPsd.svd.dim(2);
R = volPsd.svd.dim(3);
K = volPsd.svd.dim(4);
F = volPsd.svd.dim(5);
V = volPsd.svd.dim(6);
W = volPsd.svd.dim(7);
M = volPsd.svd.dim(8);

% subStr = strsplit(replace(volPsd.fspec,'.nii.gz',''),filesep); subStr = strsplit(subStr{end},'_');
subStr = strsplit(replace(volPsd.fspec,'.nii.gz',''),filesep); subStr = strsplit(subStr{end-1},'_');

spSV = permute(volPsd.svd.spSV,[5 6 8 1 2 3 4 7]); % [F V M]
switch method
    case 1
        spSV = reshape(spSV(:,:,1),[F V*1]); % [F V*M]
    case 2
        spSV = reshape(spSV,[F V*M]); % [F V*M]
end
spSV = permute(spSV,[2 1]); % [V*M F]
Ml2 = min(size(spSV));
[u,s,v] = svd(spSV,'econ','vector');

spSVl2 = permute(u,[2 1]); %[Ml2 V*M]
switch method
    case 1
        spSVl2 = reshape(spSVl2,[Ml2 V 1]); % [Ml2 V M];
    case 2
        spSVl2 = reshape(spSVl2,[Ml2 V M]); % [Ml2 V M];
end

frSVl2 = permute(v,[2 1]); %[Ml2 F*M]
frSVl2 = reshape(frSVl2,[Ml2 F]); %[Ml2 F]

GridSize = [2 3];
extraMode = [];

%% Spatial singular vector -- MAG
fMag = figure('WindowStyle','docked');
ht = tiledlayout(GridSize(1),GridSize(2)); ax = {};
ht.TileSpacing = 'tight'; ht.Padding = 'tight';
ax{end+1} = nexttile;
plot(s,'.-');
xline(prod(ht.GridSize)-1+0.5,'r')
xlabel('Ml2'); ylabel('singular value');
hold on
plot(extraMode,s(extraMode),'.r');
for ml2 = 1:(prod(ht.GridSize)-1)
    if ml2>(prod(ht.GridSize)-1-length(extraMode)); ml2 = extraMode(ml2-(prod(ht.GridSize)-1)+length(extraMode)); end
    ax{end+1} = nexttile;
    im = zeros(size(volPsd.vol2vec));
    im(volPsd.vol2vec) = mean(abs(spSVl2(ml2,:,:)),3);
    imagesc(im.*s(ml2))
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.Colormap = turbo;
    ax{end}.XAxis.Visible = 'off'; ax{end}.YAxis.Visible = 'off';
    ylabel(colorbar,'singular vector mag * singular value')
    title(['Level-2 mode ' num2str(ml2) '; method' num2str(method)])
end
axIm = [ax{2:end}];
cLim = get(axIm,'CLim');
cLim = [min([cLim{:}]) max([cLim{:}])];
set(axIm,'CLim',cLim);
title(ht,strjoin(subStr(1:3)))


% add underlay
drawnow
set(axIm,'Color','none');
ax = {};
for i = 1:length(axIm)
    ax{end+1} = axes('InnerPosition',axIm(i).Position,'DataAspectRatio',axIm(i).DataAspectRatio);
    imagesc(volPsd.imMean);
    ax{end}.Colormap = gray;
    ax{end}.XAxis.Visible = 'off'; ax{end}.YAxis.Visible = 'off';
    uistack(ax{end},'bottom')
end


axImMag = axIm;
axImMag_under = [ax{:}];


%% Spatial singular vector -- PHASE
% f = figure('SizeChangedFcn',@(src,evn) disp('Window resized'))
fPhase = figure('WindowStyle','docked');
% figure('WindowStyle','docked','KeyPressFcn',@Key_Down);
% figure('WindowStyle','docked');
ht = tiledlayout(GridSize(1),GridSize(2)); ax = {};
ht.TileSpacing = 'tight'; ht.Padding = 'tight';
ax{end+1} = nexttile;
plot(s,'.-');
xline(prod(ht.GridSize)-1+0.5,'r')
xlabel('Ml2'); ylabel('singular value');
hold on
plot(extraMode,s(extraMode),'.r');
for ml2 = 1:(prod(ht.GridSize)-1)
    if ml2>(prod(ht.GridSize)-1-length(extraMode)); ml2 = extraMode(ml2-(prod(ht.GridSize)-1)+length(extraMode)); end
    ax{end+1} = nexttile;
    im = zeros(size(volPsd.vol2vec));
    im(volPsd.vol2vec) = angle(spSVl2(ml2,:,1));
    alphaVal = mean(abs(spSVl2(ml2,:,:)),3);
    alphaVal = alphaVal - min(alphaVal);
    alphaVal = alphaVal ./ max(alphaVal);
    alphaVal = alphaVal .* (s(ml2)/s(1));
    alphaIm = zeros(size(volPsd.vol2vec));
    alphaIm(volPsd.vol2vec) = alphaVal;
    if adjustPhase
        [~,b] = max(alphaIm(:));
        im = wrapToPi(im - im(b));
    end
    hIm = imagesc(im);
    if adjustPhase
        [x,y] = ind2sub(size(alphaIm),b);
        xline(y,':w'); yline(x,':w')
    end
    hIm.AlphaData = alphaIm;
    % hIm.AlphaData = ones(size(alphaIm)).*0.5;
    ax{end}.DataAspectRatio = [1 1 1];
    ax{end}.Colormap = hsv;
    ax{end}.XAxis.Visible = 'off'; ax{end}.YAxis.Visible = 'off';
    ax{end}.Color = [0.5 0.5 0.5];
    ylabel(colorbar,'singular vector phase')
    title(['Level-2 mode ' num2str(ml2) '; method' num2str(method)])
end
axIm = [ax{2:end}];
cLim = [-pi pi];
set(axIm,'CLim',cLim);
title(ht,strjoin(subStr(1:3)))


% add underlay
drawnow
set(axIm,'Color','none');
ax = {};
for i = 1:length(axIm)
    ax{end+1} = axes('InnerPosition',axIm(i).Position,'DataAspectRatio',axIm(i).DataAspectRatio);
    imagesc(volPsd.imMean);
    ax{end}.Colormap = gray;
    ax{end}.XAxis.Visible = 'off'; ax{end}.YAxis.Visible = 'off';
    uistack(ax{end},'bottom')
end


axImPhase = axIm;
axImPhase_under = [ax{:}];


linkaxes([axImMag axImMag_under axImPhase axImPhase_under])


addlistener(axImMag  ,'XLim','PostSet',@(src,evn) set( axImMag_under'   , {'Position'} , get(axImMag  ,'Position') ));
addlistener(axImMag  ,'YLim','PostSet',@(src,evn) set( axImMag_under'   , {'Position'} , get(axImMag  ,'Position') ));
addlistener(axImPhase,'XLim','PostSet',@(src,evn) set( axImPhase_under' , {'Position'} , get(axImPhase,'Position') ));
addlistener(axImPhase,'YLim','PostSet',@(src,evn) set( axImPhase_under' , {'Position'} , get(axImPhase,'Position') ));

fMag.KeyPressFcn   = @toggleOverlay;
fPhase.KeyPressFcn = @toggleOverlay;

fMag.SizeChangedFcn   = @(src,evn) set( axImMag_under'     , {'Position'} , get(axImMag    ,'Position') );
fPhase.SizeChangedFcn = @(src,evn) set( axImPhase_under'   , {'Position'} , get(axImPhase  ,'Position') );


% switch event.Key; case 'leftarrow',  disp( 'Izquierda' ); case 'rightarrow', disp( 'Derecha' ); end
% addlistener(fMag,'YLim','PostSet',@(src,evn) set( axImPhase_under' , {'Position'} , get(axImPhase,'Position') ));
% % 'KeyPressFcn',@(src,evn) disp('key')


%% Frequency singular vector
fSV = figure('WindowStyle','docked');
ht = tiledlayout(GridSize(1),GridSize(2)); ax = {};
ht.TileSpacing = 'tight'; ht.Padding = 'tight';
ax{end+1} = nexttile;
plot(s,'.-');
xline(prod(ht.GridSize)-1+0.5,'r')
xlabel('Ml2'); ylabel('singular value');
hold on
plot(extraMode,s(extraMode),'.r');
for ml2 = 1:(prod(ht.GridSize)-1)
    if ml2>(prod(ht.GridSize)-1-length(extraMode)); ml2 = extraMode(ml2-(prod(ht.GridSize)-1)+length(extraMode)); end
    ax{end+1} = nexttile;
    x = squeeze(volPsd.svd.f);
    y = abs(frSVl2(ml2,:));
    plot(x,y.*s(ml2));
    xlabel('Hz')
    ylabel('singular vector mag * singular value')
    title(['Level-2 mode ' num2str(ml2) '; method' num2str(method)])
    grid on

    if isfield(volPsd.dsgn,'onsets')
        xline(1/mean(diff(volPsd.dsgn.onsets)),'--r')
    else
        xline(1/mean(diff(volPsd.dsgn.onsetList)),'--r')
    end
end
axIm = [ax{2:end}];
yLim = get(axIm,'YLim');
yLim = [min([yLim{:}]) max([yLim{:}])];
set(axIm,'YLim',yLim);

title(ht,strjoin(subStr(1:3)))



if saveFlag
    [a,b,~] = fileparts(fileparts(volPsd.fspec));
    saveas(fMag,  fullfile(a,[b '_L2SVD_spSvMag.fig']))
    saveas(fPhase,fullfile(a,[b '_L2SVD_spSvPhase.fig']))
    saveas(fSV,   fullfile(a,[b '_L2SVD_tSvMag.fig']))
end