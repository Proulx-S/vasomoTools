function [t2, d2, c2] = simTs(onsetsList,dursList,nframes,tr,dt,SNR,plotFlag,rndSeed)
if ~exist('SNR','var');         SNR = []; end
if ~exist('rndSeed','var'); rndSeed = []; end
if isempty(SNR);                SNR = inf; end
if isempty(rndSeed);        rndSeed = 0; end
if isempty(dursList);      dursList = ones(size(onsetsList)).*tr; end

t = 0:dt:nframes*tr;
d = zeros(size(t));
for i = 1:length(onsetsList)
    ind = t>=onsetsList(i,1) & t<onsetsList(i,1)+dursList(i,1);
    d(ind) = 1;
end

addpath(genpath('/space/takoyaki/1/users/proulxs/tools/spm12'))
c = conv(spm_hrf(dt),d); c(length(d)+1:end) = [];
rmpath(genpath('/space/takoyaki/1/users/proulxs/tools/spm12'))


if dt~=tr
    t2 = 0:tr:(nframes-1)*tr;
    c2 = interp1(t,c,t2,'nearest');
    d2 = interp1(t,d,t2,'nearest');
else
    t2 = t;
    c2 = c;
    d2 = d;
end
if SNR ~= inf
    rng(rndSeed)
    c2 = c2.*SNR + randn(size(c2)).*std(c2);
end

if plotFlag
    figure('WindowStyle','docked');
    plot(t,d,':k')
    hold on
    h1 = plot(t,c./max(c(:)),'k');

    plot(t2,d2,':r')
    hold on
    h2 = plot(t2,c2./max(c2(:)),'r');

    legend([h1 h2],{'dt' 'tr'})
end
