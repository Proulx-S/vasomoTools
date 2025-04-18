function [fVessel,fNonVessel,fVesselFig,seg,labelList] = computeVesselness(fVol,fMask,force,verbose)
    if ~exist('force','var'); force = []; end
    if isempty(force);        force = 0 ; end
    if ~exist('verbose','var'); verbose = []; end
    if isempty(verbose);        verbose = 0 ; end
    
    %%% Read in data
    mri = MRIread(fVol);
    mask = MRIread(fMask); mask = logical(mask.vol);

    %%% Define gaussian mixture
    % k = 3; labelList = {'vessel' 'brain' 'non-brain'};
    % k = 4; labelList = {'vessel' 'brain' '?' '?'};
    % k = 5; labelList = {'?' '?' '?' '?' '?'};
    % k = 6; labelList = {'?' '?' '?' '?' '?' '?'};
    k = 7; labelList = {'?' '?' '?' '?' '?' '?' '?'};
    % X = log(mri.vol(mask));
    X = log(mri.vol(6:end-4,6:end-4,:));
    X = X(:);

    %%% Plot intensity histogram
    % hFig = figure('visible','off');
    hFig = figure('WindowStyle','docked');
    hT = tiledlayout(2,2); hT.TileSpacing = "tight"; hT.Padding = 'tight';
    nexttile([2 1])
    h = histogram(X,'Normalization','pdf'); hold on

    % %%% Fit a gaussian mixture with a rician noise background
    % [model, responsibility] = fitGaussianRicianMixture(X, k, 100, 1e-6);

    %%% Fit gaussian mixture
    GMModel = fitgmdist(X,k,'Options',statset('MaxIter',1000));
    GMModel = fitgmdist(X,k);


    %%%% Lowest component is noise
    [start.mu(1,1),b] = min(GMModel.mu);
    start.Sigma(1,1,1) = GMModel.Sigma(:,:,b);
    start.PComponents(1,1) = GMModel.ComponentProportion(b);

    %%%% Narrowest component is brain
    [start.Sigma(1,1,2),b] = min(GMModel.Sigma);
    start.mu(2,1) = GMModel.mu(b);
    start.PComponents(1,2) = GMModel.ComponentProportion(b);

    %%%% Highest component is vessel
    [start.mu(3,1),b] = max(GMModel.mu);
    start.Sigma(1,1,3) = GMModel.Sigma(:,:,b);
    start.PComponents(1,3) = GMModel.ComponentProportion(b);

    %%%% Closest component to brain is also brain
    [~,b] = sort(abs(start.mu(2,1) - GMModel.mu)); b = b(2);
    start.mu(4,1) = GMModel.mu(b);
    start.Sigma(1,1,4) = GMModel.Sigma(:,:,b);
    start.PComponents(1,4) = GMModel.ComponentProportion(b);

    GMModel = fitgmdist(X,4,'Start', start);



    %%% Identify the component based on mean value
    % Sort mu in ascending order to get proper indices
    [~,b] = sort(GMModel.mu,'descend');
    b2 = zeros(size(b)); for i = 1:k; b2(i) = find(GMModel.mu(b)==GMModel.mu(i)); end
    labelList = labelList(b2);





    % % Get the index of the component with the lowest mean
    % lowestMeanIdx = b(end);
    
    % % Extract the parameters of this component
    % lowestMeanMu = GMModel.mu(lowestMeanIdx);
    % lowestMeanSigma = GMModel.Sigma(lowestMeanIdx);
    % lowestMeanProp = GMModel.ComponentProportion(lowestMeanIdx);
    
    % % Create initial parameter estimates
    % start.mu = nan(k,1);
    % start.mu(end) = lowestMeanMu;

    % % For 1D data, Sigma should be initialized properly
    % start.Sigma = nan(1,1,k);
    % for i = 1:k-1
    %     start.Sigma(1,1,i) = 1; % Default positive value for unused components
    % end
    % start.Sigma(1,1,end) = max(lowestMeanSigma, 1e-6); % Ensure positive value

    % start.PComponents = [ones(1,k-1)*(1-lowestMeanProp)/(k-1), lowestMeanProp];

    % % Refit with constraints
    % options = statset('Display', 'final');
    % GMModel2 = fitgmdist(X, k, 'Options', options, 'RegularizationValue', 1e-6, ...
    %                     'Start', start);
    
    % % Update the GMModel with the new fit
    % GMModel = GMModel2;



    


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
    seg = zeros(size(mri.vol(6:end-4,6:end-4,:)));
    % seg(mask) = idx;
    seg(:) = idx;
    for i = 1:k
        seg(seg==i) = mean(X(seg==i));
    end
    
    %%% Show base image and segmentation
    axIm  = nexttile;
    imagesc(mri.vol,[0 800]);
    ax = axIm; ax.PlotBoxAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off'; ax.Colormap = gray;
    axSeg = nexttile;
    imagesc(seg);
    ax = axSeg; ax.PlotBoxAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off'; ax.Colormap = jet;
    linkaxes([axIm axSeg]);

    
    figure('WindowStyle','docked');
    imagesc(seg);
    ax = gca; ax.PlotBoxAspectRatio = [1 1 1]; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off'; ax.Colormap = gray;
    
    
    return

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

function [model, responsibility] = fitGaussianRicianMixture(X, k, maxIter, tol)
    % Initialize parameters
    [n, d] = size(X);
    k_total = k + 1; % k Gaussian components + 1 Rician
    
    % Initialize mixing coefficients
    pi_k = ones(k_total, 1) / k_total;
    
    % Initialize Gaussian parameters (means and covariances)
    % Use k-means for initialization
    [idx, centers] = kmeans(X, k);
    mu = centers;
    sigma = zeros(d, d, k);
    for i = 1:k
        cluster_points = X(idx == i, :);
        if ~isempty(cluster_points)
            sigma(:, :, i) = cov(cluster_points) + 1e-6 * eye(d); % Add small regularization
        else
            sigma(:, :, i) = eye(d); % Default if empty cluster
        end
    end
    
    % Initialize Rician parameters (nu, sigma)
    % For Rician, we need scale parameter (sigma) and non-centrality parameter (nu)
    rician_nu = mean(X(:));
    rician_sigma = std(X(:)) / 2;
    
    % Prepare for EM algorithm
    loglik = -inf;
    responsibility = zeros(n, k_total);
    
    % EM Algorithm
    for iter = 1:maxIter
        % E-step: Calculate responsibilities
        for i = 1:k
            % Gaussian density
            responsibility(:, i) = mvnpdf(X, mu(i, :), sigma(:, :, i)) * pi_k(i);
        end
        
        % Rician density (for the k+1 component)
        for j = 1:n
            x_mag = norm(X(j, :));
            responsibility(j, k+1) = ricianpdf(x_mag, rician_nu, rician_sigma) * pi_k(k+1);
        end
        
        % Normalize responsibilities
        responsibility = responsibility ./ sum(responsibility, 2);
        
        % M-step: Update parameters
        % Update mixing coefficients
        Nk = sum(responsibility);
        pi_k = Nk / n;
        
        % Update Gaussian parameters
        for i = 1:k
            resp_i = responsibility(:, i);
            % Update mean
            mu(i, :) = (resp_i' * X) / Nk(i);
            
            % Update covariance
            diff = X - mu(i, :);
            weighted_diff = diff .* sqrt(resp_i);
            sigma(:, :, i) = (weighted_diff' * weighted_diff) / Nk(i) + 1e-6 * eye(d);
        end
        
        % Update Rician parameters using MLE or moment matching
        resp_rician = responsibility(:, k+1);
        weighted_X = X .* resp_rician;
        % This is simplified - proper Rician parameter estimation requires numerical methods
        rician_nu = sum(sum(weighted_X)) / Nk(k+1);
        rician_sigma = sqrt(sum(sum((X - rician_nu).^2 .* resp_rician)) / (d * Nk(k+1)));
        
        % Check convergence
        new_loglik = sum(log(sum(responsibility, 2)));
        if abs(new_loglik - loglik) < tol
            break;
        end
        loglik = new_loglik;
    end
    
    % Package the model
    model.pi_k = pi_k;
    model.mu = mu;
    model.sigma = sigma;
    model.rician_nu = rician_nu;
    model.rician_sigma = rician_sigma;
    model.k = k;
    model.d = d;
    visualizeDistributions(X, model)


% Visualize the empirical distribution and fitted distributions
function visualizeDistributions(X, model)
    % Extract model parameters
    pi_k = model.pi_k;
    mu = model.mu;
    sigma = model.sigma;
    rician_nu = model.rician_nu;
    rician_sigma = model.rician_sigma;
    k = model.k;
    
    % Calculate magnitudes of data points
    magnitudes = sqrt(sum(X.^2, 2));
    
    % Create figure
    figure('Name', 'Distribution Fitting Results');
    
    % Plot histogram of data magnitudes
    histogram(magnitudes, 100, 'Normalization', 'pdf', 'DisplayName', 'Data');
    hold on;
    
    % Generate points for plotting the distributions
    x_range = linspace(0, max(magnitudes)*1.1, 1000);
    
    % Plot Gaussian components
    for i = 1:k
        % For each Gaussian, we need to project to 1D for visualization
        % We'll use the magnitude of the mean as the center
        mu_mag = norm(mu(i,:));
        sigma_1d = sqrt(trace(sigma(:,:,i))/size(sigma,1)); % Approximate 1D sigma
        
        y = pi_k(i) * normpdf(x_range, mu_mag, sigma_1d);
        plot(x_range, y, 'LineWidth', 2, 'DisplayName', sprintf('Gaussian %d', i));
    end
    
    % Plot Rician component
    rician_y = pi_k(k+1) * arrayfun(@(x) ricianpdf(x, rician_nu, rician_sigma), x_range);
    plot(x_range, rician_y, 'LineWidth', 2, 'DisplayName', 'Rician');
    
    % Plot combined model
    combined_y = zeros(size(x_range));
    for i = 1:k
        mu_mag = norm(mu(i,:));
        sigma_1d = sqrt(trace(sigma(:,:,i))/size(sigma,1));
        combined_y = combined_y + pi_k(i) * normpdf(x_range, mu_mag, sigma_1d);
    end
    combined_y = combined_y + pi_k(k+1) * arrayfun(@(x) ricianpdf(x, rician_nu, rician_sigma), x_range);
    plot(x_range, combined_y, 'k--', 'LineWidth', 2, 'DisplayName', 'Combined Model');
    
    % Add labels and legend
    xlabel('Magnitude');
    ylabel('Probability Density');
    title('Fitted Mixture Model vs. Empirical Distribution');
    legend('show', 'Location', 'best');
    grid on;
    hold off;


% Helper function for Rician PDF
function p = ricianpdf(x, nu, sigma)
    % Rician PDF implementation
    % x: magnitude value
    % nu: non-centrality parameter
    % sigma: scale parameter
    
    % Ensure x is non-negative
    x = max(x, 0);
    
    % Rice distribution formula
    p = (x / sigma^2) .* exp(-(x.^2 + nu^2) / (2 * sigma^2)) .* besseli(0, x * nu / sigma^2);
    
    % Handle numerical issues
    invalid = isnan(p) | isinf(p);
    p(invalid) = eps;


