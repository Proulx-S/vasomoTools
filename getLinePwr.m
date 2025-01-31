function [A,F,p] = getLinePwr(J,tp,Fs)
% Heavily adapted from
% chronux_2_12/spectral_analysis/continuous/ftestc.m
% and
% chronux_2_12/spectral_analysis/continuous/fitlinesc.m
% to accept precomputed J.
% Sebastien Proulx Oct 24th, 2024

% J [time x trial x run x taper x freq x vox x window]
%   [N      E       R     K       F      V     W     ]
[~,~,~,K,~,~,~] = size(J);

Kodd=1:2:K;
Keven=2:2:K;

H0   = sum(tp(:,:,:,Kodd,:,:,:),1);
JH0  = sum(  J(:,:,:,Kodd,:,:,:) .* H0  ,4);
A    = JH0./sum(H0.^2,4);
Jhat = A.*H0;

num  = (K-1) .* (abs(A).^2) .* sum(H0.^2,4); % F-ratio numerator
den  = sum(abs(  J(:,:,:,Kodd,:,:,:) - Jhat  ).^2,4)  +  sum(abs( J(:,:,:,Keven,:,:,:) ).^2,4); % F-ratio denominator
F = num./den; % F-ratio

df1 = 2;
df2 = 2*K-2;
p = fcdf(F,df1,df2,'upper');


A = A.*Fs;
