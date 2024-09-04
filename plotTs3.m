function [ax,F] = plotTs3(H,volTs,volResp,dsgn,mask,respQthresh)

if ~exist('H','var');                      H = []; end
if isempty(H);                             H = figure('WindowStyle','docked'); end
if ~exist('onsets','var');            onsets = []; end
if ~exist('ondurs','var');            ondurs = []; end
if ~exist('roiInd','var');            roiInd = []; end
if ~exist('respQthresh','var');  respQthresh = []; end
if ~exist('volAnat','var');          volAnat = []; end
if ~exist('mask','var');                mask = []; end
if ~exist('dsgn','var');                dsgn = []; end
if isempty(respQthresh);         respQthresh = 1; end % respQthresh = 1 means don't threshold
if isempty(volAnat);                  volRoi = [];
                                     volMask = [];
elseif ~isempty(roiInd)
                                      volRoi = volAnat; clear volAnat
                                     volMask = [];
else
                                      volRoi = [];
                                     volMask = volAnat; clear volAnat
                             volMask.mri.vol = any(volMask.mri.vol,4);
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


% fieldToRemove = {'volAnat' 'volPsd' 'volResp' 'volTs' 'psd' 'psdGram' 'psdTrialGram' 'psdTrialGramMD' 'svd' 'svdGram' 'svdTrialGram' 'svdTrialGramMD'};
if iscell(ax{end}.UserData)
    ax{end}.UserData(end+1) = {''};
else
    ax{end}.UserData = {''};
end
% ax{end}.UserData{end}.mri = rmfield(volPsd,fieldToRemove(ismember(fieldToRemove,fields(volPsd))));
hold on

% if iscell(ax{end}.UserData)
%     ax{end}.UserData(end+1) = {''};
% else
%     ax{end}.UserData = {''};
% end
% hold on


%%% average or loop across space

%%%% simple crop
if isempty(mask)
    if ~isempty(volTs)
        mask = getCropMask(volTs);
    elseif ~isempty(volResp)
        mask = getCropMask(volResp(1).ts);
    else
        dbstack; error('X');
    end
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
        mask = logical(mask);
    end
    
    if ~isempty(volTs)
        if isMRI(volTs)
            mask = mask & getCropMask(volTs);
        else
            mask = mask & getCropMask(volTs.mri);
        end
    elseif ~isempty(volResp)
        mask = mask & getCropMask(volResp(1).ts);
    else
        dbstack; error('X');
    end
end

%%%% stat thresh
if respQthresh==0; respQthresh = 0.05; end % respQthresh = 0 does not mean anything, reverts to default 0.05
if respQthresh~=inf && respQthresh~=1
    mask = mask & volResp.Fq.vol<respQthresh;
end

legH = {};
if ~isempty(volTs)
    %%% average across runs
    volTs.vol = mean(volTs.vol,5);

    %%% average across voxels
    volTs = vol2vec(volTs,mask);
    t = volTs.t;
    ts = mean(volTs.vec,2);

    legH{end+1} = plot(t,ts,'k');
    grid on
    axis tight
    xlabel('t (sec)')
    ylabel('MR signal (a.u.)')

    legLabel = {'timeseries data'};
else
    legLabel = {};
end

%% Add response if available
if ~isempty(volResp)
    hold on
    if length(volResp)>1 && isstruct(volResp) && ~isMRI(volResp)
        for i = 1:length(volResp)
            mri = MRIread(volResp(i).ts.fspec);
            volResp(1).ts.vol = cat(6,volResp(1).ts.vol,mri.vol); clear mri
            % volResp(1).ts.vol(:,:,:,:,:,i) = mri.vol; clear mri
        end
        volResp(2:end) = [];
        if isfield(volResp.ts,'vec')
            volResp.ts = rmfield(volResp.ts,'vec');
        end

    else
        if isempty(volResp.ts.vol)
            mri = MRIread(volResp.ts.fspec);
            volResp.ts.vol = mri.vol; clear mri
            if isfield(volResp.ts,'vec')
                volResp.ts = rmfield(volResp.ts,'vec');
            end
        end
    end
    volResp.ts = vol2vec(volResp.ts,mask);
    tResp = volResp.ts.t; if ~isempty(dsgn) && ~isempty(volTs); tResp = tResp + dsgn.onsetList(1); end
    tsResp = mean(volResp.ts.vec,2);
    if size(volResp.ts.vec,4)>1
        tsRespEr = std(tsResp,[],4)./sqrt(size(volResp.ts.vec,4));
        tsResp   = mean(tsResp,4);
    end
    

    if ~isempty(volTs)
        mm = [min(tsResp) max(tsResp)];
        tsResp = tsResp-mean(mm);
        mm = [min(tsResp) max(tsResp)];
        tsResp = tsResp./mm(2)/2;
        mm = ylim;
        tsResp = tsResp.*diff(mm);
        tsResp = tsResp+mean(mm);
    end
    if exist('tsRespEr','var')
        hEr = shplot(tResp,tsResp,tsRespEr);
        delete(hEr.upper); hEr = rmfield(hEr,'upper'); delete(hEr.lower); hEr = rmfield(hEr,'lower');
        ax{end}.UserData{end}.data = hEr;
        hLine = hEr.line;
        hEr.patch.FaceAlpha = 0.1;
        hEr.patch.FaceColor = hLine.Color;
    else
        legH{end+1} = plot(tResp,tsResp,'r');
    end

    legLabel = [legLabel {'response time course'}];
end


%% Add modeled response fit if available
if ~isempty(volTs) && ~isempty(volResp) && isfield(volResp,'SPMG2')
    volTsFit = rmfield(volTs,'vec');
    hold on
    f   = replace(volResp.SPMG2.coef.fspec,'_coef.nii.gz','_fit.nii.gz');
    mri = MRIread(f);
    volTsFit.vol   = mri.vol;
    volTsFit.fspec = f; clear mri
    volTsFit = vol2vec(volTsFit);
    t   = volTsFit.t;
    fit = mean(volTsFit.vec,2);
    legH{end+1} = plot(t,fit,'b');

    legLabel = [legLabel {'SPMG2 (spm double gamma + derivative) fit'}];
end

legend([legH{:}],legLabel,'AutoUpdate','off')


if ~isempty(volTs)
    T = volTs.nframes.*volTs.tr/1000;
    xlim([0 T])
elseif ~isempty(volResp)
    T = volResp.ts.nframes.*volResp.ts.tr/1000;
    xlim([0 T-volResp.ts.tr/1000])
else
    dbstack; error('X');
end


paramStr = {['T=' num2str(T,'%0.2f') 'sec']};
if ~isempty(ondurs)
    paramStr{end+1} = ['dur=' num2str(mean(ondurs)) 'sec'];
end
if ~isempty(fMask)
    if exist(fMask,'file')
        [~,b] = fileparts(replace(fMask,'.nii.gz','')); b = strsplit(b,'_'); b = b{end};
        paramStr{end+1} = ['mask:' b];
    else
        paramStr{end+1} = ['mask:' strjoin(fMask,'+')];
    end
% elseif ~isempty(volRoi)
%     dbstack; error('code that');
end

if respQthresh~=inf && respQthresh~=1
    title(['timeseries (' strjoin(paramStr,'; ') ') averaged across voxels with Q<=' num2str(respQthresh,'%0.2f')])
else
    title(['timeseries (' strjoin(paramStr,'; ') ')'])
end


if isempty(dsgn)
    if isempty(onsets) && isfield(volTs,'dsgn') && isfield(volTs.dsgn,'onsets')
        onsets = volTs.dsgn.onsets;
    elseif isempty(onsets) && isfield(volTs,'dsgn') && isfield(volTs.dsgn,'onsetList')
        onsets = volTs.dsgn.onsetList';
    end
    if isempty(ondurs) && isfield(volTs,'dsgn') && isfield(volTs.dsgn,'ondurs')
        ondurs = volTs.dsgn.ondurs;
    elseif isempty(ondurs) && isfield(volTs,'dsgn') && isfield(volTs.dsgn,'ondurList')
        ondurs = volTs.dsgn.ondurList';
    end
else
    onsets = dsgn.onsetList;
    ondurs = dsgn.ondurList;
end

if isempty(volTs)
    onsets = 0;
    ondurs = ondurs(1);
end

if ~isempty(onsets) && isempty(ondurs)
    addOnset([],onsets)
elseif ~isempty(ondurs)
    addOndur([],onsets,ondurs)
end




ax = [ax{:}];