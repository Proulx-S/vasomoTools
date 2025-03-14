function [y,hrf_parameters] = doubleGamma(t,hrf_parameters)
% The hrf is modeled as the difference of two gamma density functions (Glover, G.H. (1999). "Deconvolution of impulse response in event-related BOLD fMRI." NeuroImage, 9:416-429). The parameters of the hrf are specified by a row vector whose elements are:
% 
% 1. AMP1:  amplitude of the first gamma density;
% 2. PEAK1: time to the peak of the first gamma density;
% 3. FWHM1: approximate FWHM of the first gamma density;
% 4. AMP2:  amplitude of the second gamma density;
% 5. PEAK2: time to the peak of the second gamma density;
% 6. FWHM2: approximate FWHM of the second gamma density;

% Heavily modified (simplified) from https://www.math.mcgill.ca/keith/fmristat/
% by Sébastien Proulx (proulxs@stanford.edu) March 3 2025



if ~exist('hrf_parameters','var') || isempty(hrf_parameters)
    hrf_parameters=[1 5.4 5.2 0.35 10.8 7.35];
    if ~exist('t','var') || isempty(t)
        y = hrf_parameters;
        return
    end
end
% Ensure t is a column vector
if size(t,2) > size(t,1)
    t = t';
end

amp1  = hrf_parameters(1);
peak1 = hrf_parameters(2);
fwhm1 = hrf_parameters(3); if fwhm1==0; fwhm1 = eps; end
amp2  = hrf_parameters(4);
peak2 = hrf_parameters(5);
fwhm2 = hrf_parameters(6); if fwhm2==0; fwhm2 = eps; end



y = 0;
if amp1~=0
    alpha1=peak1^2/fwhm1^2*8*log(2);
    beta1=fwhm1^2/peak1/8/log(2);
    gamma1=(t/peak1).^alpha1.*exp(-(t-peak1)./beta1);
    gamma1(isnan(gamma1)) = 0;
    gamma1 = gamma1./max(abs(gamma1));
    % y = y + amp1*gamma1/max(gamma1);
    y = y + amp1*gamma1;
end

if amp2~=0
    alpha2=peak2^2/fwhm2^2*8*log(2);
    beta2=fwhm2^2/peak2/8/log(2);
    gamma2=(t/peak2).^alpha2.*exp(-(t-peak2)./beta2);
    gamma2(isnan(gamma2)) = 0;
    gamma2 = gamma2./max(abs(gamma2));
    % y = y + amp2*gamma2/max(gamma2);
    y = y + amp2*gamma2;
end

% y = y./sum(y)

% if amp1~=0
%     y = y.*amp1;
% else
%     y = y.*amp2;
% end

