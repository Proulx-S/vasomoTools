function vessel = getVesselResp(vessel)
    % Uresp: temporal singular vectors scaled by singular values [copmonents x time]
    % Vresp: spatial singular vectors scaled by singular values  [components x X x Y]
    % AreaResp: area response     [1 x time]
    % VelResp : velocity response [1 x time]
    % PeakVoxResp: peak voxel response     [1   x time]
    % SurrVoxResp: surround voxel response [vox x time]

    for v = 1:length(vessel)
        vessel(v).resp = doIt(vessel(v));
    end

    function resp = doIt(vessel)
        respIm = permute(vessel.im.resp.im,[4 1 2 3]);
        wMask = vessel.polyMask{ismember(vessel.polyLabel,'peakVox')};
        zMask = vessel.polyMask{ismember(vessel.polyLabel,'dilate1p5')};
        tMask = vessel.polyMask{ismember(vessel.polyLabel,'tissue')};
        d2Mask = vessel.polyMask{ismember(vessel.polyLabel,'dilate2')};
        
        % SVD transform
        [U,S,V] = svd(respIm(:,tMask|d2Mask),'econ','vector'); % time x vox (excluding those containing other vessels)
        Uresp = permute(U,[2 1]).*S;
        Vresp = zeros([size(S,1) size(respIm,[2 3 4])]);
        Vresp(:,tMask|d2Mask) = permute(V,[2 1]).*S;
        
        % Area/velocity transform
        base = permute(mean(cat(3,vessel.im.basePolyRun.im{:}),3),[4 1 2 3]);
        respIm = respIm + base;
        wVal = respIm(:,wMask);
        zVal = respIm(:,zMask);
        tVal = respIm(:,tMask);
        AreaResp = ( size(wVal,2).*mean(wVal-mean(tVal,2),2) + size(zVal,2).*mean(zVal-mean(tVal,2),2) ) ./ mean(wVal-mean(tVal,2),2);
        AreaResp = permute(AreaResp,[2 1 3 4]);
        VelResp  = permute(mean(wVal,2),[2 1 3 4]);

        % Peak voxel
        respIm = permute(vessel.im.resp.im,[4 1 2 3]);
        pMask = vessel.polyMask{ismember(vessel.polyLabel,'peakVox')};
        PeakVoxResp = permute(respIm(:,pMask),[2 1]);
        
        % Surround voxels
        respIm = permute(vessel.im.resp.im,[4 1 2 3]);
        sMask = vessel.polyMask{ismember(vessel.polyLabel,'dilate1')};
        sMask(vessel.polyMask{ismember(vessel.polyLabel,'original')}) = false;
        SurrVoxResp = permute(respIm(:,sMask),[2 1]);
        
        % Package output
        % vessel.resp.sv          = S;
        % vessel.resp.timeResp    = Uresp;
        % vessel.resp.spaceResp   = Vresp;
        % vessel.resp.areaResp    = AreaResp;
        % vessel.resp.velResp     = VelResp;
        % vessel.resp.peakVoxResp = PeakVoxResp;
        % vessel.resp.surrVoxResp = SurrVoxResp;
        % vessel.resp.info = 'conponent/vox x time';
        % % clear resp
        resp.sv          = S;
        resp.timeResp    = Uresp;
        resp.spaceResp   = Vresp;
        resp.areaResp    = AreaResp;
        resp.velResp     = VelResp;
        resp.peakVoxResp = PeakVoxResp;
        resp.surrVoxResp = SurrVoxResp;
        resp.info = 'conponent/vox x time';