function mask = getRoiBckgrndMask(im,verbose)
% Find peaks in base image using iterative Gaussian fitting
if ~exist('verbose','var') verbose = []; end
if isempty(verbose);       verbose = 0 ; end

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

% Iteratively fit peaks
for iter = 1:maxIterations
    imToFit       = residualPrev;
    imToFitSmooth = imgaussfilt(imToFit, 1);

    % Find the maximum intensity pixel in current residual
    [maxVal, maxIdx] = max(imToFitSmooth(:));
    [peakY, peakX] = ind2sub(size(im), maxIdx);
    
    % Set initial parameter estimates
    a0     = maxVal - b;              % amplitude
    x0_0   = peakX;                   % x center
    y0_0   = peakY;                   % y center
    sx0    = 0.5;                     % x sigma
    sy0    = 0.5;                     % y sigma
    theta0 = 0;                       % rotation angle (radians)
    b0     = b;                       % background
    
    % Set parameter bounds
    lb = [0       , 1         , 1         , 0.5         , 0.5         , -pi/2, min(imSmooth(:))];     % lower bounds
    ub = [maxVal*2, size(im,2), size(im,1), size(im,1)/2, size(im,1)/2,  pi/2, max(imSmooth(:))]; % upper bounds
    
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
        break;
    end

    % Update
    b = fit_result{iter}.b;
    SSEfullPrev = SSEfull;
    DFfullPrev = DFfull;
    residualPrev = imToFit - peak_hat;

    % Visualize
    if verbose
        figure('WindowStyle','docked');
        subplot(1,3,1);
        imagesc(imToFit,[min(im(:)) max(im(:))]); axis image; colorbar;
        subplot(1,3,2);
        imagesc(imToFit_hat,[min(im(:)) max(im(:))]); axis image; colorbar;
        subplot(1,3,3);
        imagesc(imToFit-peak_hat,[min(im(:)) max(im(:))]); axis image; colorbar;
    end
end
fit_result(end) = [];


%% Fit significant peaks simultaneously
numSignificantPeaks = length(fit_result);

if numSignificantPeaks > 0
    % Create a multi-peak Gaussian model with explicit coefficients
    if numSignificantPeaks == 1
        multiPeakGaussian = fittype('a1*exp(-(((x-x1)*cos(theta1) + (y-y1)*sin(theta1))^2/(2*sx1^2) + (-(x-x1)*sin(theta1) + (y-y1)*cos(theta1))^2/(2*sy1^2))) + b', ...
            'independent', {'x', 'y'}, ...
            'dependent', 'z', ...
            'coefficients', {'a1', 'x1', 'y1', 'sx1', 'sy1', 'theta1', 'b'});
    elseif numSignificantPeaks == 2
        multiPeakGaussian = fittype('a1*exp(-(((x-x1)*cos(theta1) + (y-y1)*sin(theta1))^2/(2*sx1^2) + (-(x-x1)*sin(theta1) + (y-y1)*cos(theta1))^2/(2*sy1^2))) + a2*exp(-(((x-x2)*cos(theta2) + (y-y2)*sin(theta2))^2/(2*sx2^2) + (-(x-x2)*sin(theta2) + (y-y2)*cos(theta2))^2/(2*sy2^2))) + b', ...
            'independent', {'x', 'y'}, ...
            'dependent', 'z', ...
            'coefficients', {'a1', 'x1', 'y1', 'sx1', 'sy1', 'theta1', 'a2', 'x2', 'y2', 'sx2', 'sy2', 'theta2', 'b'});
    elseif numSignificantPeaks == 3
        multiPeakGaussian = fittype('a1*exp(-(((x-x1)*cos(theta1) + (y-y1)*sin(theta1))^2/(2*sx1^2) + (-(x-x1)*sin(theta1) + (y-y1)*cos(theta1))^2/(2*sy1^2))) + a2*exp(-(((x-x2)*cos(theta2) + (y-y2)*sin(theta2))^2/(2*sx2^2) + (-(x-x2)*sin(theta2) + (y-y2)*cos(theta2))^2/(2*sy2^2))) + a3*exp(-(((x-x3)*cos(theta3) + (y-y3)*sin(theta3))^2/(2*sx3^2) + (-(x-x3)*sin(theta3) + (y-y3)*cos(theta3))^2/(2*sy3^2))) + b', ...
            'independent', {'x', 'y'}, ...
            'dependent', 'z', ...
            'coefficients', {'a1', 'x1', 'y1', 'sx1', 'sy1', 'theta1', 'a2', 'x2', 'y2', 'sx2', 'sy2', 'theta2', 'a3', 'x3', 'y3', 'sx3', 'sy3', 'theta3', 'b'});
    else
        % For more than 3 peaks, use a simpler approach or limit to 3
        warning('More than 3 peaks detected. Limiting to first 3 peaks for simultaneous fit.');
        numSignificantPeaks = 3;
        fit_result = fit_result(1:3);
        multiPeakGaussian = fittype('a1*exp(-(((x-x1)*cos(theta1) + (y-y1)*sin(theta1))^2/(2*sx1^2) + (-(x-x1)*sin(theta1) + (y-y1)*cos(theta1))^2/(2*sy1^2))) + a2*exp(-(((x-x2)*cos(theta2) + (y-y2)*sin(theta2))^2/(2*sx2^2) + (-(x-x2)*sin(theta2) + (y-y2)*cos(theta2))^2/(2*sy2^2))) + a3*exp(-(((x-x3)*cos(theta3) + (y-y3)*sin(theta3))^2/(2*sx3^2) + (-(x-x3)*sin(theta3) + (y-y3)*cos(theta3))^2/(2*sy3^2))) + b', ...
            'independent', {'x', 'y'}, ...
            'dependent', 'z', ...
            'coefficients', {'a1', 'x1', 'y1', 'sx1', 'sy1', 'theta1', 'a2', 'x2', 'y2', 'sx2', 'sy2', 'theta2', 'a3', 'x3', 'y3', 'sx3', 'sy3', 'theta3', 'b'});
    end
    
    % Extract initial parameters from individual fits
    % Note: In the fittype, parameters come in order: a1,x1,y1,sx1,sy1,theta1, a2,x2,y2,sx2,sy2,theta2, ..., b
    startParams = [];
    lb = [];
    ub = [];
    
    for i = 1:numSignificantPeaks
        fit_i = fit_result{i};
        startParams = [startParams, fit_i.a, fit_i.x0, fit_i.y0, fit_i.sx, fit_i.sy, fit_i.theta];
        lb = [lb, 0, 1, 1, 0.5, 0.5, -pi/2];
        ub = [ub, max(im(:))*2, size(im,2), size(im,1), size(im,1)/2, size(im,1)/2, pi/2];
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

    %% Split in individual fit objects
    b = finalFit.b;
    for i = 1:numSignificantPeaks
        fit_result{i}.a     = finalFit.(sprintf('a%d', i));
        fit_result{i}.x0    = finalFit.(sprintf('x%d', i));
        fit_result{i}.y0    = finalFit.(sprintf('y%d', i));
        fit_result{i}.sx    = finalFit.(sprintf('sx%d', i));
        fit_result{i}.sy    = finalFit.(sprintf('sy%d', i));
        fit_result{i}.theta = finalFit.(sprintf('theta%d', i));
        fit_result{i}.b     = 0;
    end
    
    %% Get peak masks
    imBckgrnd = im - (finalFit(X, Y) - finalFit.b);
    thresh = std(imBckgrnd(:));

    mask = false(size(im));
    for i = 1:numSignificantPeaks
        mask = mask | (fit_result{i}(X, Y)>fit_result{i}.a*0.01);
    end

    % Dilate mask by 1 pixel in all directions (including diagonals)
    se = strel('disk', 1, 8);  % 8-connectivity includes diagonals
    mask = imdilate(mask, se);
end
