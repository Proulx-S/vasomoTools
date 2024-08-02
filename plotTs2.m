function [ax,F] = plotTs2(volTs,H,onsets,ondurs,volAnat,roiInd,respQthresh)
if ~exist('H','var');                      H = []; end
if isempty(H);                             H = figure('WindowStyle','docked'); end
if ~exist('onsets','var');            onsets = []; end
if ~exist('ondurs','var');            ondurs = []; end
if ~exist('roiInd','var');            roiInd = []; end
if ~exist('respQthresh','var');  respQthresh = []; end
if ~exist('volAnat','var');          volAnat = []; end
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
if ~isempty(roiInd)
    curMask = volRoi.mri.vol(:,:,:,roiInd);
elseif ~isempty(volMask)
    curMask = volMask.mri.vol;
else
    curMask = [];
end
if isempty(curMask)
    volTs = vol2vec(volTs);
else
    volTs = vol2vec(volTs,curMask);
end
% volTs = vol2vec(volTs);
% if isfield(volTs,'resp')
%     % % % %%%% weighted average
%     % % % weights = volTs.resp.F.vol;
%     % % % tmp = weights(weights~=0); tmp = tmp - min(tmp); tmp = tmp./mean(tmp);
%     % % % weights(weights~=0) = tmp;
%     % % % ts = squeeze(mean(volTs.vec.*weights(volTs.vol2vec)',2));
%     %%%% thresholded average
%     q  = volTs.resp.Fq.vol(volTs.vol2vec);
%     ts = volTs.vec;
%     ts = mean(ts(:,q<0.001),2);
% else
if isfield(volTs,'resp') && (respQthresh~=inf || respQthresh~=1)
    ts = squeeze(mean(volTs.vec(:,volTs.resp.Fq.vol(volTs.vol2vec)<=respQthresh),2));
else
    ts = squeeze(mean(volTs.vec,2));
end
% end
plot(squeeze(volTs.t),ts,'k')
grid on
axis tight
xlabel('t (sec)')
ylabel('MR signal (a.u.)')


%% Add response if available

if ~isempty(volTs) && isfield(volTs,'resp')
    hold on
    ts = permute(volTs.resp.ts.vol,[4 1 2 3]);
    if respQthresh~=inf && respQthresh~=1
        ts = mean(ts(:,volTs.resp.Fq.vol<=respQthresh),2);
    else
        ts = mean(ts(:,volTs.vol2vec),2);
    end
    % %%%% weighted average
    % weights = volTs.resp.F.vol;
    % tmp = weights(weights~=0); tmp = tmp - min(tmp); tmp = tmp./mean(tmp);
    % weights(weights~=0) = tmp;
    % ts = permute(volTs.resp.ts.vol.*weights,[4 1 2 3]);
    % ts = mean(ts(:,:),2);
    % %%%% thresholded average
    % q  = volTs.resp.Fq.vol(volTs.vol2vec);
    % ts = permute(volTs.resp.ts.vol,[4 1 2 3]);
    % ts = ts(:,volTs.vol2vec);
    % ts = mean(ts(:,q<0.001),2);

    t = (0:1:(volTs.resp.ts.nframes-1)).*volTs.resp.ts.tr/1000;
    t = t + volTs.dsgn.onsets(1);


    mm = [min(ts) max(ts)];
    ts = ts-mean(mm);
    mm = [min(ts) max(ts)];
    ts = ts./mm(2)/2;
    mm = ylim;
    ts = ts.*diff(mm);
    ts = ts+mean(mm);
    plot(t,ts,'r')
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

if ~isempty(volTs) && isfield(volTs,'resp') && (respQthresh~=inf && respQthresh~=1)
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