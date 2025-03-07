function pFit = fitDoubleGamma(t,y,paramMask,p0,lb_full,ub_full,h)
% fitDoubleGamma Fit a double gamma function to data
%   [pFit,HRFfit,HRFfit1,HRFfit2,p0,HRFguess,HRFguess1,HRFguess2] = fitDoubleGamma(t,y,paramMask,p0,lb_full,ub_full)
%
% Inputs:
%   t          - time vector
%   y          - data to fit
%   paramMask  - binary mask for which parameters to fit (1) or fix (0)
%   p0         - initial parameters [AMP1 PEAK1 FWHM1 AMP2 PEAK2 FWHM2]
%   lb_full    - lower bounds for all parameters
%   ub_full    - upper bounds for all parameters
%
% Outputs:
%   pFit       - fitted parameters
%   HRFfit     - fitted HRF
%   HRFfit1    - first gamma component of fitted HRF
%   HRFfit2    - second gamma component of fitted HRF
%   p0         - initial parameters used
%   HRFguess   - initial HRF guess
%   HRFguess1  - first gamma component of initial guess
%   HRFguess2  - second gamma component of initial guess

% Set default parameter mask if not provided
if nargin < 3 || isempty(paramMask)
    paramMask = [1 1 1 1 1 1];  % Default: fit all parameters
end

% Set default initial parameters if not provided
if nargin < 4 || isempty(p0)
    p0 = [-900 7 7 100 10 20];
end

% Set default bounds if not provided
if nargin < 5 || isempty(lb_full)
    lb_full = [-4000  0  0    0  5  0];
end
if nargin < 6 || isempty(ub_full)
    ub_full = [    0 15 20 1000 20 40];
end

% Extract bounds for free parameters only
lb = lb_full(paramMask == 1);
ub = ub_full(paramMask == 1);

% Extract free parameters for optimization
p0_free = p0(paramMask == 1);

% Set optimization options
options = optimoptions('lsqnonlin',...
    'Display','iter',...
    'FunctionTolerance',1e-12,...
    'StepTolerance',1e-12,...
    'OptimalityTolerance',1e-12,...
    'MaxFunctionEvaluations',100,...
    'MaxIterations',400);

% Get initial guess HRF and components
p0_1 = p0; p0_1(4) = 0;  % Zero out second gamma
p0_2 = p0; p0_2(1) = 0;  % Zero out first gamma
HRFguess = doubleGamma(t,p0);
HRFguess1 = doubleGamma(t,p0_1);
HRFguess2 = doubleGamma(t,p0_2);

% figure;
% plot(t,y,'k'); hold on;
% plotHandle = plot(t,HRFguess,'m--');
% grid on;

% Create masked objective function that only optimizes free parameters
if exist('h','var') && ~isempty(h)
    objFunMasked = @(p_free) objFunWrapper(p_free, p0, paramMask, t, y, h);
else
    objFunMasked = @(p_free) objFunWrapper(p_free, p0, paramMask, t, y);
end


% Perform the fit
[pFit_free,~,~,~,~] = lsqnonlin(objFunMasked,p0_free,lb,ub,options);

% Reconstruct full parameter set
pFit = p0;  % Start with initial values
pFit(paramMask == 1) = pFit_free;  % Update fitted parameters

% Get fitted HRF and components
pFit1 = pFit; pFit1(4) = 0;  % Zero out second gamma
pFit2 = pFit; pFit2(1) = 0;  % Zero out first gamma
HRFfit = doubleGamma(t,pFit);
HRFfit1 = doubleGamma(t,pFit1);
HRFfit2 = doubleGamma(t,pFit2);

end

function residual = objFunWrapper(p_free, p0, paramMask, t, y, plotHandle)
    % Reconstruct full parameter set
    p_full = p0;  % Start with initial values
    p_full(paramMask == 1) = p_free;  % Update free parameters
    % Calculate residual using full parameter set
    fit = doubleGamma(t,p_full);
    residual = fit - y;
    if nargin > 5 && ~isempty(plotHandle)
        plotHandle.YData = fit;
        drawnow;
    end 
end 