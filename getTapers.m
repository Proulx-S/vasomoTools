function [tp,eigs,tpNorm,N,pad] = getTapers(K,tr,N,t,pad,MDslepianAlgo)
% Notes on taper behavior (only tested with k=5 tapers):
% (1) Temporal resolution does not matter, tapers will have the same shape
% from tStart to tEnd whether they are defined with a low or high sampling
% rate.
% (2) The above is not entirely true for missing data tapers. Low vs high
% sampling rate tapers will hava similar but not exactly the same shape.
% The difference gets worse with larger sampling rate differences.
% (3) The eigenvalue decomposition used to obtain missing data tapers fails
% at very high number of data points. This has something to do with the use
% of eigs.
% (4) The above point (3) is only true if there actually is missing
% data--using the missing data taper procedure without missing data gives
% the same result as dpsschk (which uses Matlab's implementation of regular
% slepian tapers).
% (5) For missing data tapers with actual missing data, results using eigs
% or eig are the same for small N, but diverge with larger N. So eig might
% not be a solution when eigs fails to converge.

if ~exist('t','var'); t = []; end
if ~exist('MDslepianAlgo','var'); MDslepianAlgo = ''; end
if isempty(MDslepianAlgo); forceMDslepian = 0; else forceMDslepian = 1; end
if isempty(MDslepianAlgo); MDslepianAlgo = 'eigs'; end % 'eigs' or 'eig'

Fs = 1/tr;
[TW,W,K] = K2W(N/Fs,K,0);

if (isempty(t) || max(abs(diff(t,2)))<1e-10) && ~forceMDslepian
    [tp,eigs] = dpsschk([TW K],N,Fs); % check tapers
    eigs = permute(eigs,[2 1]);
    % t = permute(0:tr:((N-1)*tr),[2 1]);
else
    % t = t-t(1);
    [eigs,tp] = MDslepian(W,K,t,Fs,MDslepianAlgo);
    if isnan(eigs)
        error(['failed to compute missing data taper' newline 'possibly to many data points'])
    end
    % NB: eigs (the spectral concentration parameter) is much lower for
    % MDslepian than in dpsschk, suggesting that at least in some
    % situations, MDslepian does not perform well in terms of spectral
    % concentration. However, running MDslepian without any missing data
    % gives exactly the same taper as dpsschk, but their eigenvalues are
    % all 1/Fs times smaller...
end

% orth = tp'*tp;
% orth(eye(size(tp,2),'logical')) = 0;
% figure('WindowStyle','docked');
% imagesc(orth); ax = gca; ax.DataAspectRatio = [1 1 1]; colorbar; clim(max(abs(clim)).*[-1 1]);

%% Normalization factor
% adapted from Jackob Duckworth (/autofs/space/takoyaki_001/users/proulxs/vsmU19/fmri_linespec_JD.m, line 231)
NFFT   = max(2^(nextpow2(N)+pad),N);
tpNorm = fft(tp,NFFT)/Fs; % freq x taper <- time x taper
tpNorm = tpNorm(1,:); % taper DC

