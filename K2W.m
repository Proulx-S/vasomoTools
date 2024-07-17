function [TW,W,K] = K2W(T,K,verbose)

if ~exist('verbose','var') || isempty(verbose)
    verbose = 1;
end

TW = (K+1)/2;
W = TW/T;
if verbose
    display(['K  (number of tapers): ' num2str(K)])
    display(['W  (halfwidth)       : ' num2str(W,'%0.5f ')])
    display(['TW (time-halfwidth)  : ' num2str(TW)])
end