function [tp,eigs] = getTapers(K,tr,N,t)
if ~exist('t','var'); t = []; end

[TW,W,K] = K2W(N*tr,K,0);

if isempty(t)
    [tp,eigs] = dpsschk([TW K],N,1/tr); % check tapers
    eigs = permute(eigs,[2 1]);
    % t = permute(0:tr:((N-1)*tr),[2 1]);
else
    % t = t-t(1);
    [eigs,tp] = MDslepian(W,K,t,1/tr);
    % NB: eigs (the spectral concentration parameter) is much lower for
    % MDslepian than in dpsschk, suggesting that at least in some
    % situations, MDslepian does not perform well in terms of spectral
    % concentration. However, running MDslepian without any missing data
    % gives exactly the same taper as dpsschk, but their eigenvalues are
    % all 1/Fs times smaller...
end