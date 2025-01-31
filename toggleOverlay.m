function toggleOverlay(f,evn)

T   = findobj(f.Children,'type','TiledLayout');
if ~isempty(T)
    % overlay images are the ones not in the TiledLayout
    ax  = findobj(f.Children,'type','axes');
    im = findobj(ax(~ismember(ax,findobj(T.Children,'type','axes'))),'type','image');
else
    % overlay images are the ones with (or without) titles
    ax = findobj(f.Children,'type','axes');
    ind = get(ax,'Title'); ind = get([ind{:}],'String'); ind = ~cellfun('isempty',ind);
    im = findobj([ax(ind).Children],'Type','Image');
end

switch evn.Key
    case 'o'
        % toggle overlay on and off
        if strcmp(im(1).Visible,'on')
            set(im,'Visible','off')
        elseif strcmp(im(1).Visible,'off')
            set(im,'Visible','on')
        else
            dbstack; error('X')
        end
    case 't'
        % toggle transparency (alpha chanel) on and off
        if numel(im(1).AlphaData) == 1 && im(1).AlphaData == 1 && ~isempty(im(1).UserData)
            [im.AlphaData] = deal(im.UserData);
        elseif numel(im(1).AlphaData) ~= 1
            [im.UserData] = deal(im.AlphaData);
            set(im,'AlphaData',1);
        end
end