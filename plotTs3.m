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


switch class(H)
    case 'matlab.graphics.layout.TiledChartLayout'
        F = H.Parent;
    case 'matlab.ui.Figure'
        F = H;
    otherwise
end

figure(F);
ax = {};
ax{end+1} = nexttile;

%%% average or loop across space

%%%% simple crop
if isempty(mask)
    mask = getCropMask(volTs);
else
    if isMRI(mask)
        mask = logical(mask.vol);
    else
        mask = MRIread(mask);
        mask = logical(mask.vol);
    end
    mask = mask & getCropMask(volTs);
end

%%%% stat thresh
if respQthresh==0; respQthresh = 0.05; end % respQthresh = 0 does not mean anything, reverts to default 0.05
if respQthresh~=inf && respQthresh~=1
    mask = mask & volResp.Fq.vol<respQthresh;
end


volTs = vol2vec(volTs,mask);
t = volTs.t;
ts = mean(volTs.vec,2);

plot(t,ts,'k')
grid on
axis tight
xlabel('t (sec)')
ylabel('MR signal (a.u.)')


%% Add response if available
if ~isempty(volResp)
    hold on
    if isempty(volResp.ts.vol)
        mri = MRIread(volResp.ts.fspec);
        volResp.ts.vol = mri.vol; clear mri
        if isfield(volResp.ts,'vec')
            volResp.ts = rmfield(volResp.ts,'vec');
        end
    end
    volResp.ts = vol2vec(volResp.ts,mask);
    tResp = volResp.ts.t; if ~isempty(dsgn); tResp = tResp + dsgn.onsetList(1); end
    tsResp = mean(volResp.ts.vec,2);
    

    mm = [min(tsResp) max(tsResp)];
    tsResp = tsResp-mean(mm);
    mm = [min(tsResp) max(tsResp)];
    tsResp = tsResp./mm(2)/2;
    mm = ylim;
    tsResp = tsResp.*diff(mm);
    tsResp = tsResp+mean(mm);
    plot(tResp,tsResp,'r')
end







T = volTs.nframes.*volTs.tr/1000;
xlim([0 T])

paramStr = {['T=' num2str(T,'%0.2f') 'sec']};
if ~isempty(ondurs)
    paramStr{end+1} = ['dur=' num2str(mean(ondurs)) 'sec'];
end
if ~isempty(volMask)
    paramStr{end+1} = ['mask:' strjoin(volMask.label,'+')];
elseif ~isempty(volRoi)
    dbstack; error('code that');
end

if respQthresh~=inf && respQthresh~=1
    title(['timeseries (' strjoin(paramStr,'; ') ') averaged across voxels with Q<=' num2str(respQthresh,'%0.2f')])
else
    title(['timeseries (' strjoin(paramStr,'; ') ')'])
end



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



if ~isempty(onsets) && isempty(ondurs)
    addOnset([],onsets)
elseif ~isempty(ondurs)
    addOndur([],onsets,ondurs)
end


ax = [ax{:}];