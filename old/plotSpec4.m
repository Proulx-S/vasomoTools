function [ax,F,hLine,volPsd] = plotSpec4(H,volPsd,metricLabel,dsgn,mask,thresh,tWin)
if ~exist('H','var'); H = [];                             end
if isempty(H);        H = figure('WindowStyle','docked'); end

% if isfield(volPsd,'dsgn') && isfield(volPsd.dsgn,'f0'); f0 = volPsd.dsgn.f0; else; f0 = []; end

if ~exist('tWin','var');                                     tWin = []; end
if ~exist('metricLabel','var');                       metricLabel = []; end
if ~exist('volTs','var');                                   volTs = []; end
if ~exist('thresh','var');                                 thresh = []; end
if ~exist('dsgn','var');                                     dsgn = []; end
if ~exist('mask','var');                                     mask = []; end
if isempty(metricLabel);                              metricLabel = 'psd'; end % 'psd' 'coh'
% threshAvFlag = ismember(metricLabel,{'psd' 'psdEPC'}) && ~isempty(volTs) && isfield(volTs,'resp') && ~isempty(respQthresh) && respQthresh~=inf && respQthresh~=1;

% Get volPsd
if ~isMRI(volPsd)
    if isfield(volPsd,'mri')
        volPsdOrig = rmfield(volPsd,'mri');
        volPsd = volPsd.mri;
    end
    if ~any(ismember(volPsd.dataType,'volPsd'))
        disp('no volPsd available')
        if any(ismember(volPsd.dataType,'volTs'))
            disp('volTs available, computing volPsd from it')
            info.K = 5;
            volPsd = volPsdFullMt3([],info,volPsd,mask);
        else
            dbstack; error('x');
        end
    end
end

if isempty(dsgn)
    if isfield(volTs,'dsgn');       dsgn = volTs.dsgn;
    elseif isfield(volResp,'dsgn'); dsgn = volResp.dsgn;
    elseif isfield(volPsd,'dsgn');  dsgn = volPsd.dsgn;
    end
end

yLim = [];
switch class(H)
    case 'matlab.graphics.layout.TiledChartLayout'
        F = H.Parent;
        figure(F);
        ax = {};
        ax{end+1} = nexttile;
    case 'matlab.ui.Figure'
        F = H;
        figure(F);
        ax = {};
        ax{end+1} = nexttile;
    case 'matlab.graphics.axis.Axes'
        F = H.Parent;
        switch class(F)
            case 'matlab.graphics.layout.TiledChartLayout'
                F = F.Parent;
            otherwise
                dbstack; error('X');
        end
        figure(F);
        ax = {H};
        axes(ax{end}); hold on
        yLim = ylim;
    otherwise
end

fieldToRemove = {'volAnat' 'volPsd' 'volResp' 'volTs' 'psd' 'psdGram' 'psdTrialGram' 'psdTrialGramMD' 'svd' 'svdGram' 'svdTrialGram' 'svdTrialGramMD'};
if iscell(ax{end}.UserData)
    ax{end}.UserData(end+1) = {''};
else
    ax{end}.UserData = {''};
end
ax{end}.UserData{end}.mri = rmfield(volPsd,fieldToRemove(ismember(fieldToRemove,fields(volPsd))));
hold on


%% Select approrpiate data
switch metricLabel
    case 'psd'
        mt    = volPsd.psd;
        spc   = mt.PSD;
        label = 'spatially averaged spectrum';
    case 'coh'
        mt    = volPsd.svd;
        spc   = mt.COH;
        label = 'coherence spectrum';
    otherwise
        dbstack; error('code trhat')
end



if ismember(metricLabel,{'coh' 'cohEPC' 'cohEK'})
    if ~isempty(mask); warning('the provided mask is replaced by the one used for coherence analysis'); end
    fMask  = volPsd.svd.param.mask;
end
switch metricLabel
    case 'psd'
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
            mask = mask & getCropMask(volPsd(1));

            if ~isempty(fMask)
                [~,b,~] = fileparts(replace(fMask,'.nii.gz',''));
                mskStr = ['voxSel:' b '==1'];
            else
                mskStr = 'voxSel:?';
            end
            
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
        end

        % if thresh==0; thresh = 0.05; end % respQthresh = 0 does not mean anything, reverts to default 0.05
        % if thresh~=inf && thresh~=1
        %     mask = mask & volResp.fs.fFullQ.vol<thresh;
        % end
        %%% apply
        if any(size(spc,[1:2 4 7 8])~=1); dbstack; error('something unexpected here'); end
        spc = spc(:,:,:,:,:,mask(volPsd.vol2vec),:,:);
        spc = mean(sqrt(spc),6).^2;
        f   = mt.f;
        % spc = permute(spc(:,:,:,:,:,mask(volPsd.vol2vec),:,:),[5 6 1 2 3 4 7 8]);
        % spc = mean(sqrt(spc),2).^2;
        % f = permute(mt.f,[5 1 2 3 4 6 7 8]);

        %% average across runs
        spcN = size(spc,3);
        if spcN>1
            spcEr = std(sqrt(spc),[],3).^2;
            spc   = mean(sqrt(spc),3).^2;

            spcEr = permute(spcEr,[5 1 2 3 4 6 7 8]);
        end

        spc = permute(spc,[5 1 2 3 4 6 7 8]);
        f   = permute(f  ,[5 1 2 3 4 6 7 8]);
    case 'coh'
        %% coherence is already summarizing space
        m = 1;
        spcN = size(spc,3);
        if any(size(spc,[1 2 4 6:7])~=1); dbstack; error('something unexpected here'); end
        spc = permute(mean(spc(:,:,:,:,:,:,:,m),3),[5 1 2 3 4 6 7 8]);
        f = permute(mt.f,[5 1 2 3 4 6 7 8]);

        if ~isempty(fMask) && ischar(fMask)
            [~,b,~] = fileparts(replace(fMask,'.nii.gz',''));
            mskStr = ['voxSel:' b '==1'];
        else
            mskStr = 'voxSel:?';
        end
    otherwise
        dbstack; error('code trhat')
end




%% Plot
if spcN>1 && exist('spcEr','var')
    hEr = shplot(f,spc,spcEr);
    delete(hEr.upper); hEr = rmfield(hEr,'upper'); delete(hEr.lower); hEr = rmfield(hEr,'lower');
    ax{end}.UserData{end}.data = hEr;  
    hLine = hEr.line;
else
    h = plot(f,spc,'k');
    ax{end}.UserData{end}.data = h;
    hLine = h;
end

grid on
axis tight
xlabel('f (Hz)')
switch metricLabel
    case 'psd'
        ylabel('psd')
        ax{end}.YScale = 'log';
    case 'coh'
        ylabel('coherence')
    otherwise
        dbstack; error('code trhat')
end



%% Add spectrum of fitted model timeseries if available
if exist('volResp','var') && isfield(volResp,'SPMG2') && isfield(volResp.SPMG2,'volPsd')
    switch metricLabel
        case 'psd'
            yyaxis right
            f_fit = squeeze(volResp.SPMG2.volPsd.psd.f);
            spc_fit = zeros([length(f_fit) size(volResp.SPMG2.volPsd.vol2vec)]);
            spc_fit(:,volResp.SPMG2.volPsd.vol2vec) = volResp.SPMG2.volPsd.psd.PSD;
            spc_fit = mean(spc_fit(:,mask),2);
            hold on
            plot(f_fit,spc_fit,'b')
            set(gca,'yscale','log')
            yyaxis left
        case 'coh'
            % f_fit = squeeze(volResp.SPMG2.volPsd.svd.f);
            % spc_fit = squeeze(volResp.SPMG2.volPsd.svd.COH(:,:,:,:,:,:,:,1));
            % hold on
            % plot(f_fit,spc_fit,'b')
        otherwise
            dbstack; error('code trhat')
    end
end



K = mt.K;
T = mt.T;
[TW,W] = K2W(T,K,0);


paramStr = {['K=' num2str(K)] ['2W=' num2str(W*2,'%0.4f') 'Hz'] ['T=' num2str(T,'%0.2f') 'sec'] ['TW=' num2str(TW)]};
paramStr{end+1} = mskStr;
% switch metricLabel
%     case 'psd'
        title(['spectrum (' strjoin(paramStr,'; ') ')'],'Interpreter','none')
        % title(['spectrum ' paramStr ' cross-vox (Q<=' num2str(thresh,'%0.2f') ') mean'])
        % if thresh~=inf && thresh~=1
        %     title(['spectrum ' paramStr ' cross-vox (Q<=' num2str(thresh,'%0.2f') ') mean'])
        % else
        %     title(['spectrum ' paramStr ' cross-vox mean'])
        % end
%     case 'coh'
%         title(['coherence spectrum ' paramStr])
% end

yLim = [yLim'; spc(f>0.01)];
yLim = [min(yLim) max(yLim)];
ylim(yLim)
if exist('volResp','var') && isfield(volResp,'SPMG2') && isfield(volResp.SPMG2,'volPsd') && strcmp(metricLabel,'psd')
    yyaxis right
    ylim(yLim)
    yyaxis left
end
xlim([0 f(end)])

h = addW([],volPsd.psd);
ax{end}.UserData{end}.W = h;


% if isfield(volPsd,'psdTrialGram') && ~isempty(volPsd.psdTrialGram)
%     addFreq([],volPsd.psdTrialGram.onsetList,volPsd.psdTrialGram.durList,2)
% end
if ~isempty(dsgn)
    h = addFreq([],dsgn.onsetList,dsgn.ondurList);
    ax{end}.UserData{end}.freqLine = h;
end

% if ~isempty(f0)
%     xline(f0,'--b')
% end

if isfield(volPsd,'psdTrialGramMD') && ~isempty(volPsd.psdTrialGramMD) && ~isempty(tWin)
    dbstack; error('double-check that');
    t   = volPsd.psdTrialGramMD.t(1,1,1,1,1,1,:,1);
    f   = volPsd.psdTrialGramMD.f(1,1,1,1,:,1,1,1);
    if tWin==inf
        psd = mean(mean(volPsd.psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,:,:),6),7);
    else
        [~,wInd] = min(abs(t-tWin));
        psd = mean(volPsd.psdTrialGramMD.vec.psdPC(:,:,:,:,:,:,wInd,:),6);
    end
    plot(squeeze(f),squeeze(psd),'--k')
end


ax = [ax{:}];




