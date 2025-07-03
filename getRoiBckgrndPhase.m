function [peakInfo, watershedInfo] = getRoiBckgrndPhase(im)
% Find peaks in base image using Curve Fitting Toolbox with rotation
% Apply Gaussian smoothing to reduce noise
imSmooth = imgaussfilt(im, 1);

% Find the maximum intensity pixel
[maxVal, maxIdx] = max(imSmooth(:));
[peakY, peakX] = ind2sub(size(imSmooth), maxIdx);

% Create coordinate grids for fitting
[X, Y] = meshgrid(1:size(imSmooth,2), 1:size(imSmooth,1));

% Define 2D Gaussian function with rotation using fittype
gaussian2D = fittype('a*exp(-(((x-x0)*cos(theta) + (y-y0)*sin(theta))^2/(2*sx^2) + (-(x-x0)*sin(theta) + (y-y0)*cos(theta))^2/(2*sy^2))) + b', ...
    'independent', {'x', 'y'}, ...
    'dependent', 'z', ...
    'coefficients', {'a', 'x0', 'y0', 'sx', 'sy', 'theta', 'b'});

% Prepare data for fitting
x_data = X(:);
y_data = Y(:);
z_data = im(:);

% Set initial parameter estimates
a0 = maxVal - mean(im(:));  % amplitude
x0_0 = peakX;               % x center
y0_0 = peakY;               % y center
sx0 = 2;                    % x sigma
sy0 = 2;                    % y sigma
theta0 = 0;                 % rotation angle (radians)
b0 = mean(im(:));           % background

% Set parameter bounds
lb = [0, 1, 1, 0.5, 0.5, -pi/2, 0];     % lower bounds
ub = [maxVal*2, size(im,2), size(im,1), 10, 10, pi/2, maxVal]; % upper bounds

% Fit the 2D Gaussian with rotation
fit_result = fit([x_data, y_data], z_data, gaussian2D, ...
    'StartPoint', [a0, x0_0, y0_0, sx0, sy0, theta0, b0], ...
    'Lower', lb, ...
    'Upper', ub);

% Extract fitted parameters
a = fit_result.a;        % amplitude
x0 = fit_result.x0;      % x center
y0 = fit_result.y0;      % y center
sx = fit_result.sx;      % x sigma
sy = fit_result.sy;      % y sigma
theta = fit_result.theta; % rotation angle
b = fit_result.b;        % background

% Create peak mask based on fitted Gaussian
% Calculate rotated coordinates
x_rot = (X - x0) * cos(theta) + (Y - y0) * sin(theta);
y_rot = -(X - x0) * sin(theta) + (Y - y0) * cos(theta);

% Calculate the rotated Gaussian
gaussian_fit = a * exp(-(x_rot.^2 / (2*sx^2) + y_rot.^2 / (2*sy^2))) + b;

% Create mask (values above 1% of peak height)
threshold = a * 0.01 + b;
peakMask = gaussian_fit > threshold;

% Store parameters: [amplitude, x0, y0, sigma_x, sigma_y, theta, offset]
gaussianParams = [a, x0, y0, sx, sy, theta, b];

% Update peak location from fitted parameters
peakX = x0;
peakY = y0;
peaks = a + b;  % peak height

% Visualization
figure('WindowStyle','docked');
imagesc(im)
hold on;
plot(peakX, peakY, 'r*', 'MarkerSize', 10);
if length(gaussianParams) > 0
    % Plot fitted Gaussian contour
    threshold_contour = gaussianParams(1) * 0.01 + gaussianParams(7);
    contour(X, Y, gaussian_fit, [threshold_contour threshold_contour], 'r--', 'LineWidth', 2);
    
    % Plot rotation axes for visualization
    axis_length = max(sx, sy) * 2;
    x_axis_x = [x0 - axis_length*cos(theta), x0 + axis_length*cos(theta)];
    x_axis_y = [y0 - axis_length*sin(theta), y0 + axis_length*sin(theta)];
    y_axis_x = [x0 + axis_length*sin(theta), x0 - axis_length*sin(theta)];
    y_axis_y = [y0 - axis_length*cos(theta), y0 + axis_length*cos(theta)];
    
    plot(x_axis_x, x_axis_y, 'g-', 'LineWidth', 2);  % Major axis
    plot(y_axis_x, y_axis_y, 'b-', 'LineWidth', 2);  % Minor axis
end
hold off;
title('Peak Detection with Rotated 2D Gaussian Fitting');

end
