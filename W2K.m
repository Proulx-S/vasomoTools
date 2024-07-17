function [TW,W,K] = W2K(T,W,verbose)

if ~exist('verbose','var') || isempty(verbose)
    verbose = 1;
end

TW = T*W;
K = round(2*TW-1); if K==0; K=1; end
TW = (K+1)/2;
W = TW/T;

if verbose
    display(['W  (halfwidth)       : ' num2str(W,'%0.5f ')])
    display(['TW (time-halfwidth)  : ' num2str(TW)])
    display(['K  (number of tapers): ' num2str(K)])
end