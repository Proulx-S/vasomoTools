function plotSpecAll2(volPsd,volTs,volResp,dsgn,fMask,fUlay,Q)
saveFlag = 0;


if ~exist('volTs','var');       volTs = []; end
if ~exist('Q','var');       Q = []; end
if ~exist('dsgn','var'); dsgn = []; end
if isempty(Q);              Q = 1; end % respQthresh = 1 means don't threshold
if isempty(dsgn)
    if isfield(volPsd,'dsgn')
                         dsgn = volPsd.dsgn;
    end
end

if ~isempty(volTs)
    if ~isMRI(volTs)
        volTs = volTs.mri;
    end
end



Fpsd = figure('WindowStyle','docked');
HtPsd = tiledlayout(8,1); HtPsd.TileSpacing = 'tight'; HtPsd.Padding = 'tight';
axPsd = {};
axPsd{end+1} = plotTs3(      HtPsd,volTs,volResp                ,dsgn,fMask        ,Q   );
axPsd{end+1} = plotSpec3(    HtPsd,volPsd,'psd'                 ,dsgn,fMask,volResp,Q,[]);
axPsd{end+1} = plotSpecGram3(HtPsd,volPsd,'gram'       ,'psd'   ,dsgn,fMask,volResp,Q   );
axPsd{end+1} = plotSpecGram3(HtPsd,volPsd,'trialGram'  ,'psd'   ,dsgn,fMask,volResp,Q   );
axPsd{end+1} = plotSpecGram3(HtPsd,volPsd,'trialGram'  ,'psdEPC',dsgn,fMask,volResp,Q   );
axPsd{end+1} = nexttile;
axPsd{end+1} = plotSpecGram3(HtPsd,volPsd,'trialGramMD','psd'   ,dsgn,fMask,volResp,Q   );
axPsd{end+1} = plotSpecGram3(HtPsd,volPsd,'trialGramMD','psdEPC',dsgn,fMask,volResp,Q   );
% axPsd{end+1} = plotSpecGram(volPsd,'trialGramMD','psd'   ,HtPsd,volTs,respQthresh);
% axPsd{end+1} = plotSpecGram(volPsd,'trialGramMD','psdEPC',HtPsd,volTs,respQthresh);
% % cLim = get([axPsd{3:5}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
% % set([axPsd{3:5}],'CLim',cLim);
% % cLim = get([axPsd{7:8}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
% % set([axPsd{7:8}],'CLim',cLim);

[~,b,~] = fileparts(fileparts(volPsd.fspec));
% [~,b,~] = fileparts(replace(volPsd.fspec,'.nii.gz',''));
b = strsplit(b,'_');
title(HtPsd,strjoin(b(1:3)))

% ind = 4:5; indX = false(size(ind)); for i = 1:length(ind); indX(i) = ~isempty(findobj(axPsd{ind(i)}.Children,'Type','Image')); end; ind = ind(indX);
% if length(ind)>1
%     cLim = get([axPsd{ind}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
%     set([axPsd{ind}],'CLim',cLim);
% end

drawnow

if volPsd.param.tapers(2)~=1
    Fcoh = figure('WindowStyle','docked');
    HtCoh = tiledlayout(8,1); HtCoh.TileSpacing = 'tight'; HtCoh.Padding = 'tight';
    axCoh = {};
    axCoh{end+1} = plotTs3(      HtCoh,volTs,volResp                ,dsgn,fMask        ,Q   );
    axCoh{end+1} = plotSpec3(    HtCoh,volPsd,'coh'                 ,dsgn,fMask,volResp,Q,[]);
    axCoh{end+1} = plotSpecGram3(HtCoh,volPsd,'gram'       ,'coh'   ,dsgn,fMask,volResp,Q   );
    axCoh{end+1} = plotSpecGram3(HtCoh,volPsd,'trialGram'  ,'coh'   ,dsgn,fMask,volResp,Q   );
    axCoh{end+1} = plotSpecGram3(HtCoh,volPsd,'trialGram'  ,'cohEPC',dsgn,fMask,volResp,Q   );
    axCoh{end+1} = plotSpecGram3(HtCoh,volPsd,'trialGram'  ,'cohEK' ,dsgn,fMask,volResp,Q   );
    axCoh{end+1} = plotSpecGram3(HtCoh,volPsd,'trialGramMD','coh'   ,dsgn,fMask,volResp,Q   );
    axCoh{end+1} = plotSpecGram3(HtCoh,volPsd,'trialGramMD','cohEPC',dsgn,fMask,volResp,Q   );
end


% 
% 
% 
% Fcoh = figure('WindowStyle','docked');
% HtCoh = tiledlayout(8,1); HtCoh.TileSpacing = 'tight'; HtCoh.Padding = 'tight';
% axCoh = {};
% axCoh{end+1} = plotTs3(HtCoh,volTs,volResp,dsgn,fMask,0.01);
% if isfield(volPsd,'svd') && isfield(volPsd.svd,'COH') && ~isempty(volPsd.svd.COH)
%     axCoh{end+1} = plotSpec3(HtCoh,volPsd,'coh',fMask,volResp,0.01,dsgn,[]);
%     % axCoh{end+1} = plotSpec(volPsd,'coh',HtCoh);
% else
%     axCoh{end+1} = nexttile;
% end
% if isfield(volPsd,'svdGram')        && isfield(volPsd.svdGram,'COH') && ~isempty(volPsd.svdGram.COH)
%     axCoh{end+1} = plotSpecGram(volPsd,'gram'       ,'coh'   ,HtCoh);
% else
%     axCoh{end+1} = nexttile;
% end
% if isfield(volPsd,'svdTrialGram')   && isfield(volPsd.svdTrialGram,'vec')   && isfield(volPsd.svdTrialGram.vec,'coh')    && ~isempty(volPsd.svdTrialGram.vec.coh)
%     axCoh{end+1} = plotSpecGram(volPsd,'trialGram'  ,'coh'   ,HtCoh);
% else
%     axCoh{end+1} = nexttile;
% end
% if isfield(volPsd,'svdTrialGram')   && isfield(volPsd.svdTrialGram,'vec')   && isfield(volPsd.svdTrialGram.vec,'cohEPC') && ~isempty(volPsd.svdTrialGram.vec.cohEPC)
%     axCoh{end+1} = plotSpecGram(volPsd,'trialGram'  ,'cohEPC',HtCoh);
% else
%     axCoh{end+1} = nexttile;
% end
% if isfield(volPsd,'svdTrialGram')   && isfield(volPsd.svdTrialGram,'vec')   && isfield(volPsd.svdTrialGram.vec,'cohEK')  && ~isempty(volPsd.svdTrialGram.vec.cohEK)
%     axCoh{end+1} = plotSpecGram(volPsd,'trialGram'  ,'cohEK' ,HtCoh);
% else
%     axCoh{end+1} = nexttile;
% end
% if isfield(volPsd,'svdTrialGramMD') && isfield(volPsd.svdTrialGramMD,'vec') && isfield(volPsd.svdTrialGramMD.vec,'coh')     && ~isempty(volPsd.svdTrialGram.vec.coh)
%     axCoh{end+1} = plotSpecGram(volPsd,'trialGramMD','coh'   ,HtCoh);
% else
%     axCoh{end+1} = nexttile;
% end
% if isfield(volPsd,'svdTrialGramMD') && isfield(volPsd.svdTrialGramMD,'vec') && isfield(volPsd.svdTrialGramMD.vec,'cohEPC')  && ~isempty(volPsd.svdTrialGram.vec.cohEPC)
%     axCoh{end+1} = plotSpecGram(volPsd,'trialGramMD','cohEPC',HtCoh);
% else
%     axCoh{end+1} = nexttile;
% end



ind = 4:5;
tmp = [axCoh{ind}]; if iscell(tmp); tmp = [tmp{:}]; end
if ~any(cellfun('isempty',get(tmp,'Children')))
    indX = false(size(ind)); for i = 1:length(ind); indX(i) = ~isempty(findobj(axCoh{ind(i)}.Children,'Type','Image')); end; ind = ind(indX);
    if length(ind)>1
        cLim = get([axCoh{ind}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
        set([axCoh{ind}],'CLim',cLim);
    end
end


ind = 7:8;
tmp = [axCoh{ind}]; if iscell(tmp); tmp = [tmp{:}]; end
if ~any(cellfun('isempty',get(tmp,'Children')))
    cLim = get([axCoh{ind}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
    set([axCoh{ind}],'CLim',cLim);
end

% title(HtCoh,b{1})
title(HtCoh,strjoin(b(1:3)))

drawnow


%% save
if saveFlag
    [a,b,~] = fileparts(fileparts(volPsd.fspec));
    saveas(Fpsd,fullfile(a,[b '_tPsd.fig']))
    saveas(Fcoh,fullfile(a,[b '_tCoh.fig']))
end