function [fVessel,fNonVessel,fVesselFig,seg,labelList] = computeVesselness(fVol,fMask,force,verbose)
    if ~exist('force','var'); force = []; end
    if isempty(force);        force = 0 ; end
    if ~exist('verbose','var'); verbose = []; end
    if isempty(verbose);        verbose = 0 ; end
    
    %%% Read in data
    mri = MRIread(fVol);
    mask = MRIread(fMask); mask = logical(mask.vol);

    %%% Define gaussian mixture
    k = 3; labelList = {'vessel' 'brain' 'non-brain'};
    % k = 4; labelList = {'vessel' 'GM' 'WM' 'non-brain'};
    % k = 5; labelList = {'vessel1' 'vessel2' 'GM' 'WM' 'non-brain'};
    X = log(mri.vol(mask));

    %%% Plot intensity histogram
    hFig = figure('visible','off');
    % hFig = figure('WindowStyle','docked');
    hT = tiledlayout(2,2); hT.TileSpacing = "tight"; hT.Padding = 'tight';
    nexttile([2 1])
    h = histogram(X,'Normalization','pdf'); hold on
    
    %%% Fit gaussian mixture
    GMModel = fitgmdist(X,k);
    
    %%% Identify the component based on mean value
    % Sort mu in ascending order to get proper indices
    [~,b] = sort(GMModel.mu,'descend');
    b2 = zeros(size(b)); for i = 1:k; b2(i) = find(GMModel.mu(b)==GMModel.mu(i)); end
    labelList = labelList(b2);

    %%% Plot component pdfs
    binCent = h.BinEdges(1:end-1) - mean(diff(h.BinEdges))/2;
    plot(binCent',GMModel.pdf(binCent'))
    clear hPlot
    for i = 1:k
        n = makedist('normal',GMModel.mu(i),sqrt(GMModel.Sigma(i)));
        p = GMModel.ComponentProportion(i);
        hPlot(i) = plot(binCent,p.*pdf(n,binCent));
    end
    legend(hPlot,labelList)

    %%% Segment based on gaussian mixture
    idx = cluster(GMModel,X);
    seg = zeros(size(mri.vol));
    seg(mask) = idx;
    
    %%% Show base image and segmentation
    axIm  = nexttile;
    imagesc(mri.vol,[0 800]);
    ax = axIm; ax.PlotBoxAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off'; ax.Colormap = gray;
    axSeg = nexttile;
    imagesc(seg);
    ax = axSeg; ax.PlotBoxAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off'; ax.Colormap = jet;
    linkaxes([axIm axSeg]);
    
    %%% Write vessel mask
    fVessel = replace(fVol,'_volTs.nii.gz','_vesselMask.nii.gz');
    if force || ~exist(fVessel,'file')
        mri.vol = ismember(seg,find(ismember(labelList,'vessel')));
        MRIwrite(mri,fVessel);
    end

    %%% Write non-vessel mask
    fNonVessel = replace(fVol,'_volTs.nii.gz','_nonVesselMask.nii.gz');
    if force || ~exist(fNonVessel,'file')
        mri.vol = ismember(seg,find(~ismember(labelList,'vessel')));
        MRIwrite(mri,fNonVessel);
    end

    %%% Save figure
    fVesselFig = replace(fVessel,'.nii.gz','.fig');
    set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    savefig(hFig,fVesselFig,'compact')
    if verbose
        hFig.Visible = 'on';
        hFig.WindowStyle = 'docked';
        drawnow
    end

