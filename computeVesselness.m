function [fComp,fNonComp,labelList] = computeVesselness(fVol,fMask,force,verbose)
    if ~exist('force','var'); force = []; end
    if isempty(force);        force = 0 ; end
    if ~exist('verbose','var'); verbose = []; end
    if isempty(verbose);        verbose = 0 ; end
    
    %% Read and plot data
    %%% Read in data
    mri = MRIread(fVol);
    mask = MRIread(fMask); mask = logical(mask.vol);
    X = log(mri.vol(mask));


    % %%% Define gaussian mixture
    % k = 3; %labelList = {'vessel' 'brain' 'non-brain'};
    % % k = 4; labelList = {'vessel' 'GM' 'WM' 'non-brain'};
    % % k = 5; labelList = {'vessel1' 'vessel2' 'GM' 'WM' 'non-brain'};
    % X = log(mri.vol(mask));

    %%% Plot intensity histogram
    if verbose
        hFig = figure('WindowStyle','docked');
    else
        hFig = figure('visible','off');
    end
    
    hT = tiledlayout(2,2); hT.TileSpacing = "tight"; hT.Padding = 'tight';
    axDist = nexttile([2 1]);
    h = histogram(X,'Normalization','pdf'); hold on
    

    %% Fit and plot distribution
    %%% Fit gaussian mixture
    k = 3;
    GMModel = fitgmdist(X,k);
    
    % %%% Identify the component based on mean value
    % % Sort mu in ascending order to get proper indices
    % [~,b] = sort(GMModel.mu,'descend');
    % b2 = zeros(size(b)); for i = 1:k; b2(i) = find(GMModel.mu(b)==GMModel.mu(i)); end
    % labelList = labelList(b2);

    %%% Plot component pdfs
    binCent = h.BinEdges(1:end-1) - mean(diff(h.BinEdges))/2;
    % plot(binCent',GMModel.pdf(binCent'))
    clear hPlot
    for i = 1:k
        n = makedist('normal',GMModel.mu(i),sqrt(GMModel.Sigma(i)));
        p = GMModel.ComponentProportion(i);
        hPlot(i) = plot(binCent,p.*pdf(n,binCent));
    end
    axDist.YLim(2) = axDist.YLim(2)/4;
    

    %% Compute component probabilities and segmentation
    %%% component probability maps
    prob = posterior(GMModel, X);
    probMaps = zeros([k size(mri.vol)]);
    probMaps(:,mask) = permute(prob,[2 1]);
    probMaps = permute(probMaps,[2 3 4 1]);

    %%% segment
    idx = cluster(GMModel,X);
    segMap = zeros(size(mri.vol));
    segMap(mask) = idx;

    
    % axProb = {};
    % figure('WindowStyle','docked');
    % imagesc(log(mri.vol));
    % colormap(gray); colorbar;
    % title('image');
    % ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
    % axProb{end+1} = ax;
    % for i = 1:k
    %     figure('WindowStyle','docked');
    %     imagesc(probMaps(:,:,:,i),[0 1]);
    %     colormap(gray); colorbar;
    %     title(['Probability map for ' labelList{i} ' component']);
    %     ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
    %     axProb{end+1} = ax;
    % end
    % linkaxes([axProb{:}])
    % cLim = axProb{1}.CLim;
    % cLim(1) = 5; axProb{1}.CLim = cLim;


    %% Label components
    labelList = repmat({'?'},k,1);
    %%%% Narrowest component is brain
    [~,b] = min(GMModel.Sigma);
    labelList{b} = 'brain';
    
    %%%% Highest component is vessel
    [~,b] = max(GMModel.mu);
    labelList{b} = 'vessel';

    %%%% update legend
    legend(hPlot,labelList)

    

    %% Maps
    %%% Show base image and segmentation
    axIm  = nexttile;
    imagesc(mri.vol,[0 800]);
    ax = axIm; ax.PlotBoxAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off'; ax.Colormap = gray;
    axSeg = nexttile;
    imagesc(segMap);
    ax = axSeg; ax.PlotBoxAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
    cMap = get(hPlot,'Color'); cMap = cat(1,cMap{:}); cMap = cat(1,[0 0 0],cMap);
    ax.Colormap = cMap;
    linkaxes([axIm axSeg]);
    


    %% Write segmentation data
    fComp = cell(1,length(labelList));
    for i = 1:k
        if strcmp(labelList{i},'?'); continue; end
        fComp{i} = replace(fVol,'_volTs.nii.gz',['_' labelList{i} 'SegMask.nii.gz']);
        if force || ~exist(fComp{i},'file')
            mri.vol = ismember(segMap,find(ismember(labelList,labelList{i})));
            MRIwrite(mri,fComp{i});
        end
    end
    fNonComp = cell(1,length(labelList));
    for i = 1:k
        if strcmp(labelList{i},'?'); continue; end
        fNonComp{i} = replace(fVol,'_volTs.nii.gz',['_non' labelList{i} 'SegMask.nii.gz']);
        if force || ~exist(fNonComp{i},'file')
            mri.vol = ismember(segMap,find(~ismember(labelList,labelList{i})));
            MRIwrite(mri,fNonComp{i});
        end
    end


    % %%% Write vessel mask
    % fVessel = replace(fVol,'_volTs.nii.gz','_vesselMask.nii.gz');
    % if force || ~exist(fVessel,'file')
    %     mri.vol = ismember(segMap,find(ismember(labelList,'vessel')));
    %     MRIwrite(mri,fVessel);
    % end

    % %%% Write non-vessel mask
    % fNonVessel = replace(fVol,'_volTs.nii.gz','_nonVesselMask.nii.gz');
    % if force || ~exist(fNonVessel,'file')
    %     mri.vol = ismember(segMap,find(~ismember(labelList,'vessel')));
    %     MRIwrite(mri,fNonVessel);
    % end

    %%% Save figure
    fSegFig = replace(fVol,'_volTs.nii.gz','_seg.fig');
    if verbose
        % hFig.Visible = 'on';
        % hFig.WindowStyle = 'docked';
        % drawnow
    else
        set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    end
    savefig(hFig,fSegFig,'compact')


