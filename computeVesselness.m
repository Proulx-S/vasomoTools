function [fComp,fNonComp,fSegFig,fVol] = computeVesselness(fVol,fMask,force,verbose)
    if ~exist('force','var'); force = []; end
    if isempty(force);        force = 0 ; end
    if ~exist('verbose','var'); verbose = []; end
    if isempty(verbose);        verbose = 0 ; end
    fVol = cellstr(fVol);


    %% Main parameters and file names
    k = 3;
    tmpLabelList = {'vessel' 'brain' 'nonBrain'};

    fComp = cell(length(fVol),length(tmpLabelList));
    fNonComp = cell(length(fVol),length(tmpLabelList));
    for ii = 1:length(fVol)
        for i = 1:k
            fComp{ii,i}    = replace(fVol{ii},'_volTs.nii.gz',['_'    tmpLabelList{i} 'SegMask.nii.gz']);
            fNonComp{ii,i} = replace(fVol{ii},'_volTs.nii.gz',['_non' tmpLabelList{i} 'SegMask.nii.gz']);
        end
    end
    if length(fVol)>1
        fSegFig = strsplit(fVol{1},'_');
        fSegFig{contains(fSegFig,'run-')} = 'run-cat';
        fSegFig = strjoin(fSegFig,'_');
    else
        fSegFig = fVol{1};
    end
    fSegFig = replace(fSegFig,'_volTs.nii.gz','_seg.fig');
    if ~exist(fileparts(fSegFig),'dir')
        mkdir(fileparts(fSegFig));
    end






    if force || any(cellfun(@(x) ~exist(x,'file'),[fComp(:); fNonComp(:); cellstr(fSegFig)]))




    
    %% Read and plot data
    %%% Read in data
    mask = MRIread(fMask); mask = logical(mask.vol);
    X = cell(size(fVol));
    for R = length(fVol):-1:1
        mri = MRIread(fVol{R});
        X{R} = log(mri.vol(mask));
    end


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
    h = histogram(cat(1,X{:}),'Normalization','pdf'); hold on
    title(hT,fVol{1},'Interpreter','none');
    

    %% Fit and plot distribution
    % exception
    if strcmp(fVol{1},'/scratch/users/Proulx-S/doIt_generalPreproc/vsmDiamCenSur/prc/sub-vsmDiamCenSurP2/ses-1/acq-vfMRI_prsc-dflt/N4_av_cat_av_preproc_volTs.nii.gz')
        k = 4;
    end
    %%% Fit gaussian mixture
    rng(42, 'twister'); % Fixed seed for consistent GMM fitting
    GMModel = fitgmdist(cat(1,X{:}),k,'Options',statset('MaxIter',1000));
    
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
    cmap = get(hPlot,'Color');
    [~,b] = sort(GMModel.mu,'descend');
    % cmap = cmap(b);
    for i = 1:k
        set(hPlot(b(i)),'Color',cat(1,cmap{i}));
    end
    

    %% Compute component probabilities and segmentation
    %%% component probability maps
    prob = posterior(GMModel, cat(1,X{:}));
    idx  =   cluster(GMModel, cat(1,X{:}));

    prob = permute(reshape(permute(prob,[2 1]),[size(prob,2) nnz(mask) length(fVol)]),[2 1 3]);
    idx  = permute(reshape(permute(idx ,[2 1]),[1            nnz(mask) length(fVol)]),[2 1 3]);

    prob = mat2cell(prob,size(prob,1),size(prob,2),ones(1,length(fVol)));
    idx  = mat2cell(idx ,size(idx ,1),size(idx ,2),ones(1,length(fVol)));

    probMaps = cell(size(fVol));
    segMap   = cell(size(fVol));
    for R = 1:length(fVol)
        probMaps{R} = zeros([k size(mri.vol)]);
        probMaps{R}(:,mask) = permute(prob{R},[2 1]);
        probMaps{R} = permute(probMaps{R},[2 3 4 1]);
        segMap{R} = zeros(size(mri.vol));
        segMap{R}(mask) = idx{R};
    end

    
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

    % exception
    if strcmp(fVol,'/scratch/users/Proulx-S/doIt_generalPreproc/vsmDiamCenSur/prc/sub-vsmDiamCenSurP2/ses-1/acq-vfMRI_prsc-dflt/N4_av_cat_av_preproc_volTs.nii.gz')
        %%%% Second highest component is also vessel
        [~,b] = sort(GMModel.mu,'descend');
        labelList{b(2)} = 'vessel';
    end

    %%%% The rest
    % make sure labels are ok
    if ~all(ismember(labelList(~ismember(labelList,'?')),tmpLabelList))
        dbstack; error('labelList does not match tmpLabelList');
    end
    % fill in non-formally defined labels
    labelList(ismember(labelList,'?')) = tmpLabelList(~ismember(tmpLabelList,labelList));




    % %% Merge components
    % if strcmp(fVol,'/scratch/users/Proulx-S/doIt_generalPreproc/vsmDiamCenSur/prc/sub-vsmDiamCenSurP2/ses-1/acq-vfMRI_prsc-dflt/N4_av_cat_av_preproc_volTs.nii.gz')
    %     b = find(ismember(labelList,'vessel'));
    %     [~,bx] = max(GMModel.mu(b))
    %     bx = b(bx); % component to merge
    %     b = b(b~=bx); % component to merge to

    %     segMap(segMap==bx) = b;
    %     probMaps(:,:,:,b) = probMaps(:,:,:,b) + probMaps(:,:,:,bx);
    %     probMaps(:,:,:,bx) = nan;
    % end

    

    %%%% update legend
    legend(hPlot,labelList)

    

    %% Maps
    %%% Show base image and segmentation
    axIm  = nexttile;
    imagesc(mri.vol,[0 800]);
    ax = axIm; ax.PlotBoxAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off'; ax.Colormap = gray;
    axSeg = nexttile;
    imagesc(segMap{1});
    ax = axSeg; ax.PlotBoxAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
    cMap = get(hPlot,'Color'); cMap = cat(1,cMap{:}); cMap = cat(1,[0 0 0],cMap);
    ax.Colormap = cMap;
    % ax.Colormap = cat(1,cmap{:});
    linkaxes([axIm axSeg]);
    


    %% Write segmentation data
    for R = 1:size(fVol,1)
        for i = 1:size(fComp,2)
            if force || ~exist(fComp{R,i},'file') || ~exist(fNonComp{R,i},'file')
                mri.vol = ismember(segMap{R},find(ismember(labelList,tmpLabelList{i})));
                MRIwrite(mri,fComp{R,i});
                mri.vol = ismember(segMap{R},find(~ismember(labelList,tmpLabelList{i})));
                MRIwrite(mri,fNonComp{R,i});
            end
        end
    end

    %%% Save figure
    if verbose
        % hFig.Visible = 'on';
        % hFig.WindowStyle = 'docked';
        % drawnow
    else
        set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    end
    savefig(hFig,fSegFig,'compact')







    end
