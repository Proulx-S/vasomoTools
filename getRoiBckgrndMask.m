function [mask,f] = getRoiBckgrndMask(im,verbose)
% Find peaks in base image using iterative Gaussian fitting
if ~exist('verbose','var') verbose = []; end
if isempty(verbose);       verbose = 0 ; end

% Suppress warnings about clearing confidence bounds
warning('off', 'curvefit:fit:confBoundsCleared');

% Create coordinate grids for fitting
[X, Y] = meshgrid(1:size(im,2), 1:size(im,1));

% Initialize variables
imSmooth      = imgaussfilt(im, 1);
b             = mean(imSmooth(:)); % background
residualPrev  = im;
SSEfullPrev   = sum((residualPrev(:)-b).^2); % null model is just the mean of the image
DFfullPrev    = numel(im) - 1;
maxIterations = 5;


% Define 2D Gaussian function with rotation using fittype
gaussian2D = fittype('a*exp(-(((x-x0)*cos(theta) + (y-y0)*sin(theta))^2/(2*sx^2) + (-(x-x0)*sin(theta) + (y-y0)*cos(theta))^2/(2*sy^2))) + b', ...
    'independent', {'x', 'y'}, ...
    'dependent', 'z', ...
    'coefficients', {'a', 'x0', 'y0', 'sx', 'sy', 'theta', 'b'});

if verbose
    f{1} = figure('WindowStyle','docked');
    subplot(1,3,1);
    imagesc(im,[min(im(:)) max(im(:))]); axis image; colormap gray; colorbar;
    subplot(1,3,2);
    axis image; colormap gray; colorbar;
    subplot(1,3,3);
    axis image; colormap gray; colorbar;
    title(sprintf('peak %d',0));
else
    f = {};
end

% Iteratively fit peaks
for iter = 1:maxIterations
    imToFit       = residualPrev;
    imToFitSmooth = imgaussfilt(imToFit, 1);

    % Find the maximum intensity pixel in current residual
    [maxVal, maxIdx] = max(imToFitSmooth(:));
    [peakY, peakX] = ind2sub(size(im), maxIdx);

    if maxVal - b < std(residualPrev(:))
        break;
    end
    
    % Set initial parameter estimates
    a0     = maxVal - b;              % amplitude
    x0_0   = peakX;                   % x center
    y0_0   = peakY;                   % y center
    sx0    = 0.5;                     % x sigma
    sy0    = 0.5;                     % y sigma
    theta0 = 0;                       % rotation angle (radians)
    b0     = b;                       % background
    
    % Set parameter bounds
    lb = [0       , 1         , 1         , 0.25         , 0.25         , -pi/2, min(imSmooth(:))];     % lower bounds
    ub = [maxVal*2, size(im,2), size(im,1), size(im,1)/3, size(im,1)/3,  pi/2, max(imSmooth(:))]; % upper bounds
    
    % Fit the 2D Gaussian with rotation
    fit_result{iter} = fit([X(:), Y(:)], imToFit(:), gaussian2D, ...
        'StartPoint', [a0, x0_0, y0_0, sx0, sy0, theta0, b0], ...
        'Lower', lb, ...
        'Upper', ub);
    
    % Compute peak prediction
    imToFit_hat = fit_result{iter}(X, Y);
    peak_hat    = imToFit_hat - fit_result{iter}.b;

    % Compute F
    SSEfull    = sum((imToFit(:)-imToFit_hat(:)).^2);
    SSEreduced = SSEfullPrev;
    DFfull    = numel(im) - (7 * iter); % 7 parameters per peak
    DFreduced = DFfullPrev;
    F_stat = ((SSEreduced - SSEfull) / (DFreduced - DFfull)) / (SSEfull / DFfull);
    df1 = DFreduced - DFfull; % degrees of freedom for numerator
    df2 = DFfull; % degrees of freedom for denominator
    p_value = 1 - fcdf(F_stat, df1, df2);

    if verbose
        fit_result{iter}
        F_stat
        p_value
    end

    if p_value > 0.05
        fit_result(end) = [];
        break;
    end

    % Update
    b = fit_result{iter}.b;
    SSEfullPrev = SSEfull;
    DFfullPrev = DFfull;
    residualPrev = imToFit - peak_hat;

    % Visualize
    if verbose
        f{iter+1} = figure('WindowStyle','docked');
        subplot(1,3,1);
        imagesc(imToFit,[min(im(:)) max(im(:))]); axis image; colormap gray; colorbar;
        subplot(1,3,2);
        imagesc(imToFit_hat,[min(im(:)) max(im(:))]); axis image; colormap gray; colorbar;
        subplot(1,3,3);
        imagesc(imToFit-peak_hat,[min(im(:)) max(im(:))]); axis image; colormap gray; colorbar;
        title(sprintf('peak %d',iter));
    else
        f{iter} = [];
    end
end



%% Fit significant peaks simultaneously
numSignificantPeaks = length(fit_result);

if numSignificantPeaks > 0
    % Create a multi-peak Gaussian model with explicit coefficients
    % Build the fittype expression dynamically for up to maxIterations peaks
    if numSignificantPeaks > maxIterations
        warning('More than %d peaks detected. Limiting to first %d peaks for simultaneous fit.', maxIterations, maxIterations);
        numSignificantPeaks = maxIterations;
        fit_result = fit_result(1:maxIterations);
    end
    
    % Build the expression string
    expr = '';
    coeffs = {};
    
    for i = 1:numSignificantPeaks
        if i > 1
            expr = [expr ' + '];
        end
        expr = [expr sprintf('a%d*exp(-(((x-x%d)*cos(theta%d) + (y-y%d)*sin(theta%d))^2/(2*sx%d^2) + (-(x-x%d)*sin(theta%d) + (y-y%d)*cos(theta%d))^2/(2*sy%d^2)))', ...
            i, i, i, i, i, i, i, i, i, i, i)];
        
        % Add coefficients for this peak
        coeffs = [coeffs, {sprintf('a%d', i), sprintf('x%d', i), sprintf('y%d', i), ...
            sprintf('sx%d', i), sprintf('sy%d', i), sprintf('theta%d', i)}];
    end
    
    % Add background term
    expr = [expr ' + b'];
    coeffs = [coeffs, 'b'];
    
    % Create the fittype
    multiPeakGaussian = fittype(expr, ...
        'independent', {'x', 'y'}, ...
        'dependent', 'z', ...
        'coefficients', coeffs);
    
    % Extract initial parameters from individual fits
    % Note: In the fittype, parameters come in order: a1,x1,y1,sx1,sy1,theta1, a2,x2,y2,sx2,sy2,theta2, ..., b
    startParams = [];
    lb = [];
    ub = [];
    
    for i = 1:numSignificantPeaks
        fit_i = fit_result{i};
        startParams = [startParams, fit_i.a, fit_i.x0, fit_i.y0, fit_i.sx, fit_i.sy, fit_i.theta];
        lb = [lb, 0, 1, 1, 0.25, 0.25, -pi/2];
        ub = [ub, max(im(:))*2, size(im,2), size(im,1), size(im,1)/3, size(im,1)/3, pi/2];
    end
    
    % Add background parameter at the end
    startParams = [startParams, b];
    lb = [lb, min(im(:))];
    ub = [ub, max(im(:))];
    
    % Fit all peaks simultaneously
    finalFit = fit([X(:), Y(:)], im(:), multiPeakGaussian, ...
        'StartPoint', startParams, ...
        'Lower', lb, ...
        'Upper', ub);

    %% Store final parameters from simultaneous fit
    b = finalFit.b;
    final_params = cell(1, numSignificantPeaks);
    for i = 1:numSignificantPeaks
        final_params{i}.a     = finalFit.(sprintf('a%d', i));
        final_params{i}.x0    = finalFit.(sprintf('x%d', i));
        final_params{i}.y0    = finalFit.(sprintf('y%d', i));
        final_params{i}.sx    = finalFit.(sprintf('sx%d', i));
        final_params{i}.sy    = finalFit.(sprintf('sy%d', i));
        final_params{i}.theta = finalFit.(sprintf('theta%d', i));
        final_params{i}.b     = 0;
    end
    
    %% Get peak masks
    imBckgrnd = im - (finalFit(X, Y) - finalFit.b);
    thresh = std(imBckgrnd(:));

    mask = false(size(im));
    for i = 1:numSignificantPeaks
        % Calculate peak contribution directly using final parameters
        param = final_params{i};
        x_rot = (X - param.x0) * cos(param.theta) + (Y - param.y0) * sin(param.theta);
        y_rot = -(X - param.x0) * sin(param.theta) + (Y - param.y0) * cos(param.theta);
        peak_contribution = param.a * exp(-(x_rot.^2 / (2*param.sx^2) + y_rot.^2 / (2*param.sy^2)));
        
        mask = mask | (peak_contribution > param.a * 0.05);
    end

    % Dilate mask by 1 pixel in all directions (including diagonals)
    se = strel('disk', 1, 8);  % 8-connectivity includes diagonals
    mask = imdilate(mask, se);
    mask = ~mask;



    if verbose
        f{end+1} = figure('WindowStyle','docked');
        subplot(1,3,1);
        imagesc(im,[min(im(:)) max(im(:))]); axis image; colormap gray; colorbar;
        subplot(1,3,2);
        imagesc(finalFit(X,Y),[min(im(:)) max(im(:))]); axis image; colormap gray; colorbar;
        subplot(1,3,3);
        imagesc(mask); axis image; colormap gray; colorbar;
        title(sprintf('peak %s','all'));
    else
        f = {};
    end
    
else
    mask = false(size(im));
end

% Re-enable warnings
warning('on', 'curvefit:fit:confBoundsCleared');
