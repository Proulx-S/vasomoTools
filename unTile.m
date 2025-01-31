function F2 = unTile(F)

T = findobj(F.Children,'type','TiledLayout');
ax = findobj(T.Children,'type','axes');

F2 = figure('WindowStyle','docked');
