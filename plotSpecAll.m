function plotSpecAll(volPsd,volTs,volResp,respQthresh)
saveFlag = 1;

if ~exist('respQthresh','var'); respQthresh = []; end
if isempty(respQthresh);        respQthresh = 1; end % respQthresh = 1 means don't threshold

Fpsd = figure('WindowStyle','docked');
HtPsd = tiledlayout(8,1); HtPsd.TileSpacing = 'tight'; HtPsd.Padding = 'tight';
axPsd = {};
axPsd{end+1} = plotTs2(volTs,HtPsd);
axPsd{end+1} = plotSpec(volPsd,'psd',HtPsd,[],volTs);
axPsd{end+1} = plotSpecGram(volPsd,'gram'       ,'psd'   ,HtPsd,volTs,respQthresh);
axPsd{end+1} = plotSpecGram(volPsd,'trialGram'  ,'psd'   ,HtPsd,volTs,respQthresh);
axPsd{end+1} = plotSpecGram(volPsd,'trialGram'  ,'psdEPC',HtPsd,volTs,respQthresh);
axPsd{end+1} = nexttile;
axPsd{end+1} = plotSpecGram(volPsd,'trialGramMD','psd'   ,HtPsd,volTs,respQthresh);
axPsd{end+1} = plotSpecGram(volPsd,'trialGramMD','psdEPC',HtPsd,volTs,respQthresh);
% cLim = get([axPsd{3:5}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
% set([axPsd{3:5}],'CLim',cLim);
% cLim = get([axPsd{7:8}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
% set([axPsd{7:8}],'CLim',cLim);

[~,b,~] = fileparts(fileparts(volPsd.fspec));
% [~,b,~] = fileparts(replace(volPsd.fspec,'.nii.gz',''));
b = strsplit(b,'_');
title(HtPsd,strjoin(b(1:3)))


Fcoh = figure('WindowStyle','docked');
HtCoh = tiledlayout(8,1); HtCoh.TileSpacing = 'tight'; HtCoh.Padding = 'tight';
axCoh = {};
axCoh{end+1} = plotTs2(volTs,HtCoh,[],[],[],[],respQthresh);
if isfield(volPsd,'svd') && isfield(volPsd.svd,'COH') && ~isempty(volPsd.svd.COH)
    axCoh{end+1} = plotSpec(volPsd,'coh',HtCoh);
else
    axCoh{end+1} = nexttile;
end
if isfield(volPsd,'svdGram')        && isfield(volPsd.svdGram,'COH') && ~isempty(volPsd.svdGram.COH)
    axCoh{end+1} = plotSpecGram(volPsd,'gram'       ,'coh'   ,HtCoh);
else
    axCoh{end+1} = nexttile;
end
if isfield(volPsd,'svdTrialGram')   && isfield(volPsd.svdTrialGram,'vec')   && isfield(volPsd.svdTrialGram.vec,'coh')    && ~isempty(volPsd.svdTrialGram.vec.coh)
    axCoh{end+1} = plotSpecGram(volPsd,'trialGram'  ,'coh'   ,HtCoh);
else
    axCoh{end+1} = nexttile;
end
if isfield(volPsd,'svdTrialGram')   && isfield(volPsd.svdTrialGram,'vec')   && isfield(volPsd.svdTrialGram.vec,'cohEPC') && ~isempty(volPsd.svdTrialGram.vec.cohEPC)
    axCoh{end+1} = plotSpecGram(volPsd,'trialGram'  ,'cohEPC',HtCoh);
else
    axCoh{end+1} = nexttile;
end
if isfield(volPsd,'svdTrialGram')   && isfield(volPsd.svdTrialGram,'vec')   && isfield(volPsd.svdTrialGram.vec,'cohEK')  && ~isempty(volPsd.svdTrialGram.vec.cohEK)
    axCoh{end+1} = plotSpecGram(volPsd,'trialGram'  ,'cohEK' ,HtCoh);
else
    axCoh{end+1} = nexttile;
end
if isfield(volPsd,'svdTrialGramMD') && isfield(volPsd.svdTrialGramMD,'vec') && isfield(volPsd.svdTrialGramMD.vec,'coh')     && ~isempty(volPsd.svdTrialGram.vec.coh)
    axCoh{end+1} = plotSpecGram(volPsd,'trialGramMD','coh'   ,HtCoh);
else
    axCoh{end+1} = nexttile;
end
if isfield(volPsd,'svdTrialGramMD') && isfield(volPsd.svdTrialGramMD,'vec') && isfield(volPsd.svdTrialGramMD.vec,'cohEPC')  && ~isempty(volPsd.svdTrialGram.vec.cohEPC)
    axCoh{end+1} = plotSpecGram(volPsd,'trialGramMD','cohEPC',HtCoh);
else
    axCoh{end+1} = nexttile;
end

ind = 3:6; indX = false(size(ind)); for i = 1:length(ind); indX(i) = ~isempty(findobj(axCoh{ind(i)}.Children,'Type','Image')); end; ind = ind(indX);
if length(ind)>1
    cLim = get([axCoh{ind}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
    set([axCoh{ind}],'CLim',cLim);
end

cLim = get([axCoh{7:8}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
set([axCoh{7:8}],'CLim',cLim);

% title(HtCoh,b{1})
title(HtCoh,strjoin(b(1:3)))


%% save
if saveFlag
    [a,b,~] = fileparts(fileparts(volPsd.fspec));
    saveas(Fpsd,fullfile(a,[b '_tPsd.fig']))
    saveas(Fcoh,fullfile(a,[b '_tCoh.fig']))
end