function [volPsd,volTs,volResp,volAct] = plotSpecAll3(volPsd,volTs,volResp,volAct,dsgn,fMask,fUlay,Q)
saveFlag = 0;


if ~exist('volTs','var');     volTs = []; end
if ~exist('volResp','var'); volResp = []; end
if ~exist('volAct','var');   volAct = []; end
if ~exist('Q','var');             Q = []; end
if ~exist('dsgn','var');       dsgn = []; end
if ~exist('fMask','var');     fMask = []; end
if isempty(Q); Q = 0.05; end % respQthresh = 1 means don't threshold
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


if ~isempty(volAct)
    if ischar(volAct.fs.fFullQ)
        volAct.fs.fFullQ = MRIload2(volAct.fs.fFullQ);
    end
    thresh.map  = volAct.fs.fFullQ;
else
    if ischar(volResp.fs.fFullQ)
        volResp.fs.fFullQ = MRIload2(volResp.fs.fFullQ);
    end
    thresh.map  = volResp.fs.fFullQ;
end
thresh.val  = Q;
thresh.sign = -1; % -1->smaller than thresh passes; +1->greater than thresh passes



Fpsd = figure('WindowStyle','docked');
HtPsd = tiledlayout(8,1); HtPsd.TileSpacing = 'tight'; HtPsd.Padding = 'tight';
axPsd = {};
[axPsd{end+1},~,  volTs ,volResp] = plotTs4(      HtPsd,volTs,volResp                ,dsgn,fMask,thresh);
[axPsd{end+1},~,~,volPsd        ] = plotSpec4(    HtPsd,volPsd              ,'psd'   ,dsgn,fMask,thresh);
axPsd{end+1}                      = plotSpecGram4(HtPsd,volPsd,'gram'       ,'psd'   ,dsgn,fMask,thresh);
axPsd{end+1}                      = plotSpecGram4(HtPsd,volPsd,'trialGram'  ,'psd'   ,dsgn,fMask,thresh);
axPsd{end+1}                      = plotSpecGram4(HtPsd,volPsd,'trialGram'  ,'psdEPC',dsgn,fMask,thresh);
axPsd{end+1}                      = nexttile;
axPsd{end+1}                      = plotSpecGram4(HtPsd,volPsd,'trialGramMD','psd'   ,dsgn,fMask,thresh);
axPsd{end+1}                      = plotSpecGram4(HtPsd,volPsd,'trialGramMD','psdEPC',dsgn,fMask,thresh);
% axPsd{end+1} = plotSpecGram(volPsd,'trialGramMD','psd'   ,HtPsd,volTs,respQthresh);
% axPsd{end+1} = plotSpecGram(volPsd,'trialGramMD','psdEPC',HtPsd,volTs,respQthresh);
% % cLim = get([axPsd{3:5}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
% % set([axPsd{3:5}],'CLim',cLim);
% % cLim = get([axPsd{7:8}],'CLim'); cLim = cat(1,cLim{:}); cLim = [min(cLim(:)) max(cLim(:))];
% % set([axPsd{7:8}],'CLim',cLim);

[~,b,~] = fileparts(fileparts(volPsd.fspec));
b = strsplit(b,'_'); tmp = strsplit(b{1},'-'); if strcmp(tmp{1},'sub'); b = strjoin(b(1:3)); else; b = strsplit(volPsd.fspec,filesep); b = b{end}; end
title(HtPsd,b,'interpreter','none')

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
    [axCoh{end+1},~,  volTs ,volResp] = plotTs4(      HtCoh,volTs,volResp                ,dsgn,fMask,thresh);
    [axCoh{end+1},~,~,volPsd        ] = plotSpec4(    HtCoh,volPsd              ,'coh'   ,dsgn             );
    axCoh{end+1}                      = plotSpecGram4(HtCoh,volPsd,'gram'       ,'coh'   ,dsgn             );
    axCoh{end+1}                      = plotSpecGram4(HtCoh,volPsd,'trialGram'  ,'coh'   ,dsgn             );
    axCoh{end+1}                      = plotSpecGram4(HtCoh,volPsd,'trialGram'  ,'cohEPC',dsgn             );
    axCoh{end+1}                      = plotSpecGram3(HtCoh,volPsd,'trialGram'  ,'cohEK' ,dsgn             );
    axCoh{end+1}                      = plotSpecGram4(HtCoh,volPsd,'trialGramMD','coh'   ,dsgn             );
    axCoh{end+1}                      = plotSpecGram4(HtCoh,volPsd,'trialGramMD','cohEPC',dsgn             );
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

title(HtCoh,b,'interpreter','none')

drawnow


%% save
if saveFlag
    [a,b,~] = fileparts(fileparts(volPsd.fspec));
    saveas(Fpsd,fullfile(a,[b '_tPsd.fig']))
    saveas(Fcoh,fullfile(a,[b '_tCoh.fig']))
end