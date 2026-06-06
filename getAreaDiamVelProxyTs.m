function [vessel,fAll] = getAreaDiamVelProxyTs(vessel,tValAvFlag)
    % Extract area, diameter and velocity proxy timeseries (per run) from the
    % raw vessel timeseries vessel.im.ts, mirroring how respArea and respVel
    % are computed from vessel.im.resp in getVesselResp.
    %
    % Unlike getVesselResp, the raw timeseries already carries absolute
    % signal levels (no deconvolved response), so the per-run baseline is
    % NOT added back before forming the area/velocity estimates.
    %
    % tsArea: area timeseries     [1 x time], one vec per run (cell)
    % tsVel : velocity timeseries [1 x time], one vec per run (cell)
    % tsD   : diameter timeseries [1 x time], one vec per run (cell),
    %         derived from the area assuming a circular cross-section:
    %         D = 2*sqrt(A/pi)
    if nargin<2; tValAvFlag = true; end

    for v = 1:length(vessel)
        [vessel(v).im.tsArea,vessel(v).im.tsVel,vessel(v).im.tsD,fAll] = doIt(vessel(v));
    end

    function [tsArea,tsVel,tsD,fAll] = doIt(vessel)
        fAll = {};

        wMask = vessel.polyMask{ismember(vessel.polyLabel,'peakVox')};
        zMask = vessel.polyMask{ismember(vessel.polyLabel,'dilate1p5')}; zMask(wMask) = false;
        tMask = vessel.polyMask{ismember(vessel.polyLabel,'tissue')};
        wN = nnz(wMask);
        zN = nnz(zMask);

        nRun = length(vessel.im.ts.im);

        % Area/velocity transform (per run)
        tsArea = vessel.im.ts;
        tsArea.fName = '';
        tsArea.maskResp.wMask = wMask;
        tsArea.maskResp.zMask = zMask;
        tsArea.maskResp.tMask = tMask;
        tsArea.im = [];
        tsArea.im2vec = [];
        tsArea.vec = cell(1,nRun);
        tsArea.info = 'vox x time, one cell per run';
        tsArea.info2 = ['(Nw*(Sw-St)+Nz*(Sz-St)) / (Sw-St)' newline...
                                   'Nw: number of intravascular voxels' newline...
                                   'Nz: number of surrounding voxels' newline...
                                   'Sw: mean signal in intravascular voxels' newline...
                                   'Sz: mean signal in surrounding voxels' newline...
                                   'St: mean signal in tissue voxels'   ];

        tsVel = vessel.im.ts;
        tsVel.fName = '';
        tsVel.maskResp.wMask = wMask;
        tsVel.maskResp.zMask = zMask;
        tsVel.maskResp.tMask = tMask;
        tsVel.im = [];
        tsVel.im2vec = wMask;
        tsVel.vec = cell(1,nRun);
        tsVel.info = 'vox x time, one cell per run';
        tsVel.info2 = 'Sw: mean over intravascular voxels (actually just the peak voxel for now)';

        tsD = vessel.im.ts;
        tsD.fName = '';
        tsD.maskResp.wMask = wMask;
        tsD.maskResp.zMask = zMask;
        tsD.maskResp.tMask = tMask;
        tsD.im = [];
        tsD.im2vec = [];
        tsD.vec = cell(1,nRun);
        tsD.info = 'vox x time, one cell per run';
        tsD.info2 = 'D = 2*sqrt(A/pi): vessel diameter from area, assuming a circular cross-section';

        for r = 1:nRun
            tsIm = permute(vessel.im.ts.im{r},[4 1 2 3]); % time x X x Y x Z
            wVal = mean(tsIm(:,wMask),2);
            zVal = mean(tsIm(:,zMask),2);
            tVal = mean(tsIm(:,tMask),2);
            if tValAvFlag
                tVal = mean(tVal,1);
            end
            % Variant 1:
            % AreaTs = ( wN.*(wVal-tVal) + zN.*(zVal-tVal) ) ./ (wVal-tVal);
            % Variant 2: allows bounding f from 0 to 1.
            f = (zVal-tVal)./(wVal-tVal);
            if any(f>1) || any(f<0)
                % fHandle{end+1} = figure;
                % histogram(f)
                warning('getAreaDiamVelProxyTs:fOutOfBounds', ...
                    'vessel %d, run %d: %d/%d surround-fraction samples outside [0,1] (%d<0, %d>1) and were clamped.', ...
                    v,r,nnz(f<0 | f>1),numel(f),nnz(f<0),nnz(f>1));
            end
            fAll{end+1} = f;
            f = min(max(f,0),1);
            AreaTs = wN + zN.*f;
            tsArea.vec{r} = permute(AreaTs,[2 1]);
            tsVel.vec{r}  = permute(mean(wVal,2),[2 1 3 4]);
            tsD.vec{r}    = 2*sqrt(tsArea.vec{r}/pi);
        end
    end
end
