function [physRoi,physTs] = getPhysRoi(physTs,chanList,dsgn)
if ~exist('chanList','var'); chanList = []             ; end
if isempty(chanList);        chanList = {'card' 'resp'}; end

% poly: [1×1 polyshape]
% mask: [400×400 logical]
% label: 'vesselAll'
% com: [200.9350 219.7329]
% im: [1×1 struct]
% psd: [1×1 struct]
% svd: [1×1 struct]
% psdTrialGramMD: [1×1 struct]
% svdTrialGramMD: [1×1 struct]

for r = 1:size(physTs,1)
    [physRoi(r,:),physTsX(r,:)] = doIt(physTs(r,:),chanList,dsgn); physTs(r,:).vec = []; physTs(r,:).mri = [];
end
physTs = physTsX;


function [physRoi,physTs] = doIt(physTs,chanList,dsgn)

% physRoi = rmfield(physTs,'vec');
physRoi.poly = [];
physRoi.mask      = contains(physTs.chanLabel,chanList) & ~isnan(physTs.vec(1,:));
physRoi.label = 'phys';
physRoi.com = [];
physRoi.im = [];

if ~isfield(physTs,'nFrame')
    physTs.nFrame = size(physTs.vec,1);
end
if ~isfield(physTs,'dsgn') && exist('dsgn','var')
    physTs.dsgn = dsgn;
end
if ~isfield(physTs,'vol2vec')
    physTs.vol2vec = physRoi.mask;
    physTs = vol2vec(vec2vol(physTs));
end
if isfield(physTs,'mri')
    physTs.info.mri = physTs.mri;
    physTs = rmfield(physTs,'mri');
end

physRoi.ts = physTs;
physRoi.ts.vec       = physTs.vec      (:,physRoi.mask,:,:);
physRoi.ts.t0        = physTs.t0       (:,physRoi.mask,:,:);
physRoi.ts.chanLabel = chanList                           ;
physRoi.ts.Fs        = physTs.Fs       (:,physRoi.mask,:,:);

return


%% 
verboseThis = 9;
K   = 10;
W   = [];
win = inf;
physTs.dsgn   = dsgn;
physTs.nFrame = size(physTs.vec,1);

physPsd = runFullMT3(physTs,W,K,win,[],[],[],[],[],[],verboseThis,[],[])';


physTs = [];
physTs(end+1).label = 'phys';
physTs(end).poly = [];
physTs(end).mask = [];
physTs(end).com = [];
physTs(end).im = [];
for metric = {'psd' 'coh' 'harmPwr' 'harmF' 'harmP'}
    % 'time/freq' 'vox' 'taper/mode' 'run'
    switch metric
        case 'psd'
            physTs(end).vec.(metric).vec = permute(physPsd.psd.PSD     ,[5 6 8 2 3 4 7 1]);
            physTs(end).vec.(metric).f   = permute(physPsd.psd.f       ,[5 6 8 2 3 4 7 1]);
        case 'coh'
            physTs(end).vec.(metric).vec = permute(physPsd.svd.COH     ,[5 6 8 2 3 4 7 1]);
            physTs(end).vec.(metric).f   = permute(physPsd.svd.f       ,[5 6 8 2 3 4 7 1]);
        case 'harmPwr'
            physTs(end).vec.(metric).vec = permute(physPsd.harm.linePwr,[5 6 8 2 3 4 7 1]);
            physTs(end).vec.(metric).f   = permute(physPsd.harm.f      ,[5 6 8 2 3 4 7 1]);
        case 'harmF'
            physTs(end).vec.(metric).vec = permute(physPsd.harm.lineF  ,[5 6 8 2 3 4 7 1]);
            physTs(end).vec.(metric).f   = permute(physPsd.harm.f      ,[5 6 8 2 3 4 7 1]);
        case 'harmP'
            physTs(end).vec.(metric).vec = permute(physPsd.harm.lineP  ,[5 6 8 2 3 4 7 1]);
            physTs(end).vec.(metric).f   = permute(physPsd.harm.f      ,[5 6 8 2 3 4 7 1]);
        otherwise
            dbstack; error('X');
    end
end


