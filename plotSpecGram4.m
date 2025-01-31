function [ax,F] = plotSpecGram3(H,volPsd,timeLabel,metricLabel,dsgn,mask,thresh)
if ~exist('H','var'); H = [];                             end
if isempty(H);        H = figure('WindowStyle','docked'); end

if isfield(volPsd,'dsgn') && isfield(volPsd.dsgn,'f0'); f0 = volPsd.dsgn.f0; else; f0 = []; end

if ~exist('timeLabel','var');     timeLabel = [];                             end
if ~exist('metricLabel','var'); metricLabel = [];                             end
if ~exist('volResp','var');         volResp = []; end
if ~exist('thresh','var'); thresh = []; end
if ~exist('dsgn','var');               dsgn = []; end
if ~exist('mask','var');               mask = []; end
if isempty(timeLabel);            timeLabel = 'gram'; end % 'gram' 'trialGram' 'trialGramMD'
if isempty(metricLabel);        metricLabel = 'psd'; end % 'psd' 'psdEPC' 'coh' 'cohEPC' 'cohEK'
if isempty(dsgn)
    if isfield(volTs,'dsgn');       dsgn = volTs.dsgn;
    elseif isfield(volResp,'dsgn'); dsgn = volResp.dsgn;
    elseif isfield(volPsd,'dsgn');  dsgn = volPsd.dsgn;
    end
end
% threshAvFlag = ismember(metricLabel,{'psd' 'psdEPC'}) && ~isempty(volTs) && isfield(volTs,'resp') && ~isempty(respQthresh) && respQthresh~=inf && respQthresh~=1;

switch class(H)
    case 'matlab.graphics.layout.TiledChartLayout'
        F = H.Parent;
    case 'matlab.ui.Figure'
        F = H;
    case 'matlab.graphics.axis.Axes'
        switch class(H.Parent)
            case 'matlab.graphics.layout.TiledChartLayout'
                F = H.Parent.Parent;
            case 'matlab.ui.Figure'
                F = H.Parent;
        end

    otherwise
end

figure(F);
switch class(H)
    case 'matlab.graphics.axis.Axes'
        ax = {H};
    otherwise
        ax = {};
        ax{end+1} = nexttile;
end

%% Select approrpiate data
switch timeLabel
    case 'gram'
        switch metricLabel
            case 'psd'
                mt     = volPsd.psdGram;
                if isempty(mt); return; end
                spcGrm = mt.PSD;
                label  = 'spectrogram';
            case 'coh'
                mt     = volPsd.svdGram;
                if isempty(mt); return; end
                spcGrm = mt.COH;
                label  = 'coherogram';
            otherwise
                dbstack; error('code trhat')
        end
    case 'trialGram'
    %     case 'av'
    %     title(['evoked spectrogram ' paramStr])
    % case 'pc'
    %     title(['phase-locked evoked spectrogram ' paramStr])
    % case 'MD'
    %     title(['discontinuous-taper evoked spectrogram ' paramStr])
    % case 'MDpc'
    %     title(['discontinuous-taper phase-locked evoked spectrogram ' paramStr])
        switch metricLabel
            case 'psd'
                mt    = volPsd.psdTrialGram;
                if isempty(mt); return; end
                spcGrm   = mt.vec.psd;
                label = 'evoked spectrogram';
            case 'psdEPC'
                mt    = volPsd.psdTrialGram;
                if isempty(mt); return; end
                spcGrm   = mt.vec.psdPC;
                label = 'phase-locked evoked spectrogram';
            case 'coh'
                mt    = volPsd.svdTrialGram;
                if isempty(mt); return; end
                spcGrm   = mt.vec.coh;
                label = 'evoked coherogram';
            case 'cohEPC'
                mt    = volPsd.svdTrialGram;
                if isempty(mt); return; end
                spcGrm   = mt.vec.cohEPC;
                label = 'phase-locked evoked coherogram';
            case 'cohEK'
                mt    = volPsd.svdTrialGram;
                if isempty(mt); return; end
                spcGrm   = mt.vec.cohEK;
                label = 'evoked cross-trial coherogram';
            otherwise
                dbstack; error('code trhat')
        end
    case 'trialGramMD'
        switch metricLabel
            case 'psd'
                mt    = volPsd.psdTrialGramMD;
                if isempty(mt); return; end
                spcGrm   = mt.vec.psd;
                label = 'discontinuous-taper evoked spectrogram';
            case 'psdEPC'
                mt    = volPsd.psdTrialGramMD;
                if isempty(mt); return; end
                spcGrm   = mt.vec.psdPC;
                label = 'discontinuous-taper phase-locked evoked spectrogram';
            case 'coh'
                mt    = volPsd.svdTrialGramMD;
                if isempty(mt); return; end
                spcGrm   = mt.vec.coh;
                label = 'discontinuous-taper evoked coherogram';
            case 'cohEPC'
                mt    = volPsd.svdTrialGramMD;
                if isempty(mt); return; end
                spcGrm   = mt.vec.cohEPC;
                label = 'discontinuous-taper phase-locked evoked coherogram';                
            otherwise
                dbstack; error('code trhat')
        end
    otherwise
        dbstack; error('code trhat')
end


if ismember(metricLabel,{'coh' 'cohEPC' 'cohEK'})
    if ~isempty(mask); warning('the provided mask is replaced by the one used for coherence analysis'); end
    fMask = volPsd.((['svd' upper(timeLabel(1)) timeLabel(2:end)])).param.mask;
    if ~isempty(thresh); warning('the provided thresh is no applied since the coherence analysis is already done'); end
end
if contains(metricLabel,'psd')
    %% average across space
    %%% simple crop
    if isempty(mask)
        mask = getCropMask(volPsd);
    else
        if isMRI(mask)
            fMask = mask.fspec;
            mask = logical(mask.vol);
        elseif ischar(mask)
            fMask = mask;
            mask = MRIread(mask);
            mask = logical(mask.vol);
        else
            fMask = [];
        end

        if ~isempty(fMask)
            [~,b,~] = fileparts(replace(fMask,'.nii.gz',''));
            mskStr = ['voxSel:' b '==1'];
        else
            mskStr = 'voxSel:?';
        end
        

        mask = mask & getCropMask(volPsd);
    end

    %%% stat thresh
    if ~isempty(thresh)
        if ischar(thresh.map) || isempty(thresh.map.vol)
            thresh.map = MRIload3(thresh.map,[],[],0);
        end
        if thresh.sign<0
            signLabel = '<=';
            thresh.mask = thresh.map.vol<=thresh.val;
        else
            signLabel = '>=';
            thresh.mask = thresh.map.vol>=thresh.val;
        end
        mask = mask & thresh.mask;


        
        [~,b,~] = fileparts(replace(thresh.map.fspec,'.nii.gz',''));
        mskStr = [mskStr '&' b signLabel num2str(thresh.val,'%0.3f')];
    else
        mskStr = '?';
    end
    % if thresh==0; thresh = 0.05; end % respQthresh = 0 does not mean anything, reverts to default 0.05
    % if thresh~=inf && thresh~=1
    %     mask = mask & volResp.fs.fFullQ.vol<thresh;
    % end

    %%% apply
    if any(size(spcGrm,[1:4  8])~=1); dbstack; error('something unexpected here'); end
    spcGrm = spcGrm(:,:,:,:,:,mask(volPsd.vol2vec),:,:);
    spcGrm = permute(spcGrm,[5 7 6 1 2 3 4 8]);
    spcGrm = mean(spcGrm,3);

elseif contains(metricLabel,'coh')
    %% coherence already summarizes across space
    m = 1;
    if any(size(spcGrm,[1:4 6])~=1); dbstack; error('something unexpected here'); end
    spcGrm = spcGrm(:,:,:,:,:,:,:,m);
    spcGrm = permute(spcGrm,[5 7 6 1 2 3 4 8]);

    if ~isempty(fMask) && ischar(fMask)
        [~,b,~] = fileparts(replace(fMask,'.nii.gz',''));
        mskStr = ['voxSel:' b '==1'];
    else
        mskStr = 'voxSel:?';
    end
else
    dbstack; error('X');
end

% spcGrm1 = mt.vec.psd(:,:,:,:,:,mask(volPsd.vol2vec),:,:);
% spcGrm2 = mt.vec.psdPC(:,:,:,:,:,mask(volPsd.vol2vec),:,:);
% histogram(spcGrm1(:)-spcGrm2(:))
% nnz(spcGrm1(:)==spcGrm2(:))/numel(spcGrm2(:))
% histogram(spcGrm1(:))

% volPsd.vol = zeros([size(spc,5) size(volPsd.vol2vec,1:3)]);
% volPsd.vol(:,volPsd.vol2vec) = permute(spc,[5 6 1 2 3 4 7 8]);
% volPsd.vol = permute(volPsd.vol,[2 3 4 1]);
% volPsd = vol2vec(volPsd,mask,1);
% 
% 
% %% spatial average
% if threshAvFlag
%     spcGrm = permute(mean(spcGrm(:,:,:,:,:,volTs.resp.Fq.vol(volPsd.vol2vec)<=respQthresh,:,1),6),[5 7 2 8 1 3 4 6]);
% else
%     spcGrm = permute(mean(spcGrm(:,:,:,:,:,:,:,1),6),[5 7 2 8 1 3 4 6]);
% end
f   = permute(mt.f       ,[5 7 2 8 1 3 4 6]);
Fs  = mt.param.Fs;
t   = permute(mean(mt.t(:,1,:,:,:,:,:,1),1),[5 7 2 8 1 3 4 6]);
K   = mt.K;
T   = mt.T;
[TW,W] = K2W(T,K,0);

imagesc(t,f,spcGrm)

xlabel('time (s)')
ylabel('freq (Hz)')
switch metricLabel
    case {'psd' 'psdEPC'}
        ylabel(colorbar,'psd');
        ax{end}.ColorScale = 'log';
    case {'coh' 'cohEPC' 'cohEK'}
        ylabel(colorbar,'coherence');
    otherwise
        dbstack; error('code trhat')
end




% if contains(timeLabel,'trial')
%     mt.t(:,1,:,:,:,:,:)
% else
runDur = volPsd.t(end);
% end
xlim([0 runDur])

paramStr = {['K=' num2str(K)] ['2W=' num2str(W*2,'%0.4f') 'Hz'] ['T=' num2str(T,'%0.2f') 'sec'] ['TW=' num2str(TW)]};
paramStr{end+1} = mskStr;
% if contains(metricLabel,'psd')
    title([label ' (' strjoin(paramStr,'; ') ')'],'Interpreter','none')
    % if thresh~=inf && thresh~=1
    %     title([label ' ' paramStr ' cross-vox (Q<=' num2str(thresh,'%0.2f') ') mean'])
    % else
    %     title([label ' ' paramStr ' cross-vox mean'])
    % end
% else
%     title([label ' ' paramStr])
% end

cLim = [min(spcGrm(:)) max(spcGrm(:))];
if contains(metricLabel,'psd')
    cLim(1) = min(reshape(spcGrm(f>0.1,:),[numel(spcGrm(f>0.1,:)) 1])); % should define lower limit based on the lowest detectable frequency given the missing data taper used
elseif contains(metricLabel,'coh')
    if strcmp(metricLabel,'cohEK')
        cLim(1) = 1/(mt.K*mt.E);
    else
        cLim(1) = 1/mt.K;
    end
end
clim(cLim)





addWin([],mt)
addW([],mt)




% if isfield(volPsd,'psdTrialGram') && ~isempty(volPsd.psdTrialGram)
%     addFreq([],volPsd.psdTrialGram.onsetList,volPsd.psdTrialGram.durList,1)
% end
if ~isempty(dsgn)
    addFreq([],dsgn.onsetList,dsgn.ondurList)
end

% if ~isempty(f0)
%     yline(f0,'--b')
% end


addOnset([],dsgn.onsetList')



ax = [ax{:}];


%% %%%%%%%%%%%%%%%%%%
% Add to timeseries %
%%%%%%%%%%%%%%%%%% %%

axTs = findobj(allchild(F.Children),'type','axes'); ttl = get(axTs,'Title'); if ~iscell(ttl); ttl = {ttl}; end; ttl = get([ttl{:}],'String');
axTs = axTs(contains(ttl,'timeseries'));
if ~isempty(axTs)

    %%% delete previous window size visual elements (magentat lines)
    hLine = findobj(axTs.Children,'type','Line'); %hLine = {hLine(:)};
    mLine = get(hLine,'Color'); if iscell(mLine); mLine = cat(1,mLine{:}); end
    mLine = all(cat(1,mLine)==[1 0 1],2);
    delete(hLine(mLine));


    % % if length(hLine)==1
    % %     hLine = {hLine};
    % % end
    % %
    % %     mLine = get(hLine,'Color');
    % % else
    % %     mLine = {get([axTs.Children(:)],'Color')};
    % % end
    % mLine = get(hLine,'Color'); if ~iscell(mLine); mLine = {mLine}; end
    % ind = all(cat(1,mLine{:})==[1 0 1],2);
    % mLine(ind)
    %
    % mLine = hLine(all(cat(1,mLine{:})==[1 0 1],2));
    % delete(mLine);
    % % if length(axTs.Children)>1
    % %     mLine = get(hLine,'Color');
    % % else
    % %     mLine = {get([axTs.Children(:)],'Color')};
    % % end
    % % mLine = axTs.Children(all(cat(1,mLine{:})==[1 0 1],2));
    % % delete(mLine);

    %%% add window size
    addWin(axTs,mt)
end
