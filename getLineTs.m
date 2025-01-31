function ts = getLineTs(A,f,fInd,p,N,Fs)
% Heavily adapted from
% chronux_2_12/spectral_analysis/continuous/ftestc.m
% and
% chronux_2_12/spectral_analysis/continuous/fitlinesc.m
% to accept precomputed J.
% Sebastien Proulx Oct 24th, 2024

% !!!If p is specified, fInd becomes the p threshold alpha!!!



% A [time x trial x run x taper x freq x vox x window]
%   [N      E       R     K       F      V     W     ]
[~,~,~,~,~,V,~] = size(A);
if ~isempty(p)
    % Voxel specific frequencies, determined as p<alpha
    alpha = fInd;
    fInd = cell([1 1 1 1 1 V 1]);
    for v = 1:V
        fInd{v} = p(:,:,:,:,:,v)<alpha;
        ts(v,:) = ...
            permute(     A(:,:,:,:,fInd{v},v) ,[1 5 2 3 4 6 7])   *   exp( permute(f(:,:,:,:,fInd{v},:,:),[5 1 2 3 4 6 7]) *  1i*2*pi*(0:N-1) ./ Fs ) + ...
            permute(conj(A(:,:,:,:,fInd{v},v)),[1 5 2 3 4 6 7])   *   exp( permute(f(:,:,:,:,fInd{v},:,:),[5 1 2 3 4 6 7]) * -1i*2*pi*(0:N-1) ./ Fs ) ;
    end
else
    % Same frequencies for all voxels
    for v = 1:V
        ts(v,:) = ...
            permute(     A(:,:,:,:,fInd,v) ,[1 5 2 3 4 6 7])   *   exp( permute(f(:,:,:,:,fInd,:,:),[5 1 2 3 4 6 7]) *  1i*2*pi*(0:N-1) ./ Fs ) + ...
            permute(conj(A(:,:,:,:,fInd,v)),[1 5 2 3 4 6 7])   *   exp( permute(f(:,:,:,:,fInd,:,:),[5 1 2 3 4 6 7]) * -1i*2*pi*(0:N-1) ./ Fs ) ;
    end
end

