function tf = pickTimeFreq(hFig)
if ~exist('hFig','var') || isempty(hFig); hFig = gcf;
else; figure(hFig); end

%% User input
disp('Pick time-frequency coordinate on figure')
[tf.t, tf.f] = ginput(1);

%% Get info
h = findobj(hFig.Children,'type','Axes');
hS = findobj(h.Children,'type','Surface'); hS = flip(hS,1);
hL = findobj(h.Children,'type','Line'); hL = reshape(hL,[2 length(hL)/2])'; hL = flip(hL,1);
info = get(hS,'UserData'); info = [info{:}];

%% Find time and run
T = cat(1,info.t);
[~,tf.tInd] = min(abs(T - tf.t)); tf.t = T(tf.tInd);

tf.rInd = false(size(info));
for I = 1:length(info)
    tf.rInd(I) = any(ismember(info(I).t,tf.t));
end
tf.rInd = find(tf.rInd);

T = info(tf.rInd).t;
[~,tf.tInd] = min(abs(T - tf.t)); tf.t = T(tf.tInd);

%% Find frequency
F = info(tf.rInd).f;
[~,tf.fInd] = min(abs(F - tf.f)); tf.f = F(tf.fInd);


if isempty(hL(tf.rInd,1).UserData)
    %% Move the reference cross
    for i = 1:2; hL(tf.rInd,i).YData = hL(tf.rInd,i).YData - mean(hL(tf.rInd,i).YData) + tf.f; hL(tf.rInd,i).UserData = 1; end
    for i = 1:2; hL(tf.rInd,i).XData = hL(tf.rInd,i).XData - mean(hL(tf.rInd,i).XData) + tf.t; end
else
    %% Make a new reference cross
    % dbstack; error('code that');
    errorbar(tf.t,tf.f,info(tf.rInd).fBW/2,info(tf.rInd).fBW/2,info(tf.rInd).tBW/2,info(tf.rInd).tBW/2,'MarkerEdgeColor','none','MarkerFaceColor','none','CapSize',0,'Color','m','LineWidth',1.5);
end



