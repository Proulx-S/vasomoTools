function toggleOverlay(f,evn)

ax = findobj(f.Children,'type','axes');
ind = get(ax,'Title'); ind = get([ind{:}],'String'); ind = ~cellfun('isempty',ind);
im = findobj([ax(ind).Children],'Type','Image');


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