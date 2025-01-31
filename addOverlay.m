function axOl = addOverlay(axUl)

F = axUl.Parent;
while ~strcmp(F.Type,'figure'); F = F.Parent; end

axOl = cell(size(axUl));
for i = 1:numel(axUl)
    drawnow
    axOl{i} = axes(F,'InnerPosition',axUl(i).Position,'DataAspectRatio',axUl(i).DataAspectRatio,'Color','none');
    axOl{i}.XAxis.Visible = axUl(i).XAxis.Visible;
    axOl{i}.YAxis.Visible = axUl(i).YAxis.Visible;
    axOl{i}.XTick = axUl(i).XTick;
    axOl{i}.YTick = axUl(i).YTick;
    axOl{i}.XLim = axUl(i).XLim;
    axOl{i}.YLim = axUl(i).YLim;
    axOl{i}.YDir = axUl(i).YDir;
    linkaxes([axOl{i} axUl(i)]);
    addlistener(axUl(i),'XLim','PostSet',@(src,evn) set(axOl{i},{'Position'},get(axUl(i),'Position')));
    addlistener(axUl(i),'YLim','PostSet',@(src,evn) set(axOl{i},{'Position'},get(axUl(i),'Position')));
    F.SizeChangedFcn   = @(src,evn) set(axOl{i},{'Position'},get(axUl(i),'Position') );
    hold(axOl{i},'on')
end
drawnow
axOl = reshape([axOl{:}],size(axUl));







% % % % addlistener(axImMag  ,'XLim','PostSet',@(src,evn) set( axImMag_under'   , {'Position'} , get(axImMag  ,'Position') ));
% % % % addlistener(axImMag  ,'YLim','PostSet',@(src,evn) set( axImMag_under'   , {'Position'} , get(axImMag  ,'Position') ));
% % % % addlistener(axImPhase,'XLim','PostSet',@(src,evn) set( axImPhase_under' , {'Position'} , get(axImPhase,'Position') ));
% % % % addlistener(axImPhase,'YLim','PostSet',@(src,evn) set( axImPhase_under' , {'Position'} , get(axImPhase,'Position') ));
% % % % 
F.KeyPressFcn   = @toggleOverlay;
% % % % fPhase.KeyPressFcn = @toggleOverlay;
% % % % 
% % % % fMag.SizeChangedFcn   = @(src,evn) set( axImMag_under'     , {'Position'} , get(axImMag    ,'Position') );
% % % % fPhase.SizeChangedFcn = @(src,evn) set( axImPhase_under'   , {'Position'} , get(axImPhase  ,'Position') );


function toggleOverlay(f,evn)

T  = findobj(f.Children,'type','TiledLayout');
ax = findobj(f.Children,'type','axes');
axUL = ax(ismember(ax,findobj(T.Children,'type','axes')));
axOL = ax(~ismember(ax,findobj(T.Children,'type','axes')));
imOL   = findobj(axOL,'type','image');
polyUL = findobj(axUL,'type','polygon');

switch evn.Key
    case 'o'
        % toggle overlay on and off
        if strcmp(imOL(1).Visible,'on')
            set(imOL  ,'Visible','off')
            set(polyUL,'Visible','on')
        elseif strcmp(imOL(1).Visible,'off')
            set(imOL  ,'Visible','on')
            set(polyUL,'Visible','off')
        else
            dbstack; error('X')
        end
    case 't'
        % toggle transparency (alpha chanel) on and off
        if numel(imOL(1).AlphaData) == 1 && imOL(1).AlphaData == 1 && ~isempty(imOL(1).UserData)
            [imOL.AlphaData] = deal(imOL.UserData);
        elseif numel(imOL(1).AlphaData) ~= 1
            [imOL.UserData] = deal(imOL.AlphaData);
            set(imOL,'AlphaData',1);
        end
end