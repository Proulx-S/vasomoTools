function funPsdNorm = normPSD(funPsd,normMethod)
% normMethod: 'psdNoise_to1', 'psdAv_to1', 'psdNoise_toImAv', 'psdAv_toImAv' or 'none'

funPsd = vol2vec(funPsd);

%% Compute noise floor
noiseInd = logical(funPsd.psd.f>4 & funPsd.psd.f<funPsd.psd.f(round(0.95*end)));
funPsd.psd.noiseFloor = mean(funPsd.vec(noiseInd,:),1)';
funPsd.psd.noiseFloorThresh = prctile(funPsd.vec(noiseInd,:),95)';
funPsd.psd.noiseFloorInd = noiseInd';
funPsd.psd.aboveNoiseInd = funPsd.vec'>funPsd.psd.noiseFloorThresh;

%% Perform spatial normalization
switch normMethod
    case 'lowFreqTo1'
        lowFreqInd = funPsd.psd.f>0.01 & funPsd.psd.f<0.04;
        normFact = mean(funPsd.vec(lowFreqInd,:),1);
    case 'cardiacTo1'
        cardiacInd = funPsd.psd.f>0.7927 & funPsd.psd.f<1.091;
        normFact = mean(funPsd.vec(cardiacInd,:),1);
    case 'psdNoise_toImAv'
        error('double-check that for linear averaging')
        normFact = exp( mean(log(funPsd.vec(noiseInd,:)),1) );
        normFact = normFact .* exp(mean(log(normFact),2));
    case 'psdAv_toImAv'
        error('double-check that for linear averaging')
        normFact = exp( mean(log(funPsd.vec),1) );
        normFact = normFact .* exp(mean(log(normFact),2));
    case {'psdNoise_to1' 'psdNoiseTo1_thenLfSlope'}
        normFact = mean(funPsd.vec(noiseInd,:),1);
    case 'psdAv_to1'
        normFact = exp( mean(log(funPsd.vec),1) );
    case 'none'
        sz = size(funPsd.vec);
        normFact = ones([1 sz(2)]);
    otherwise
        error('')
end

%%% Apply normalization
% switch normMethod
%     case {'cardiacTo1' 'lowFreqTo1'}
%         funPsdNorm = funPsd;
%         funPsdNorm.vec = funPsdNorm.vec-1;
%         funPsdNorm.vec = funPsdNorm.vec./normFact;
%         funPsdNorm.vec = funPsdNorm.vec+1;
%     otherwise
        funPsdNorm = funPsd;
        funPsdNorm.vec = funPsd.vec./normFact;
% end

switch normMethod
    case 'psdNoiseTo1_thenLfSlope'
        f = funPsdNorm.psd.f';
        psd = funPsdNorm.vec;

        %%% Define the range within which to fit the low-frequency slope
        lfSlopeRange = [
%             0.01 0.1
            0    0.01
            0.06 0.08
%             0.17 0.2
            ];
        lfSlopeRange(lfSlopeRange==0) = f(find(f==0)+1); % make sure f=0 is not included
        ind = false(size(f));
        for i = 1:size(lfSlopeRange,1)
            ind(f>=lfSlopeRange(i,1) & f<=lfSlopeRange(i,2)) = true;
        end

        %%% Fit a linear slope in loglog space, for each voxel
        lfSlp = nan(size(psd));
        a = nan(size(f));
        b = nan(size(f));
        parfor ii = 1:size(psd,2)
            [lfSlp(:,ii),a(ii),b(ii)] = logloglin(f,psd(:,ii),ind,noiseInd);
        end

        %%% Define a new frequency- and voxel-specific normalization factor
        normFact2 = lfSlp;
        % cut the slope at its intersection with the noise floor
        normFact2(normFact2<1) = 1;

%         %%% Scale high frequency noise back to 1
%         normFact2 = normFact2./mean(normFact2(noiseInd,:),1);

% ii = 0;
% while 1
%     ii = ii+1
%         voxInd = false(1,size(psd,2));
%         voxInd(ii) = true;
        voxInd = true(1,size(psd,2));
        figure('WindowStyle','docked');
        plot(f,mean(psd(:,voxInd),2),'k'); hold on
        tmp = mean(normFact2(:,voxInd),2);
        plot(f,tmp,'--r')
        tmp(~ind) = nan;
        plot(f,tmp,'r')
        psd2 = psd./normFact2;
        hPlot = plot(f,mean(psd2(:,voxInd),2));
        plot(f,ones(size(f)),':k')
        tmp = mean(psd2,2);
        plot(f([2 end]),tmp([2 end]),':','Color',hPlot.Color)
        ax = gca;
        ax.YScale = 'log';
        ax.XScale = 'log';
%         keyboard
% end
        
        funPsdNorm.vec = funPsdNorm.vec./normFact2;
    otherwise
end

%%% Output
funPsdNorm.vecNorm = normFact;
funPsdNorm.psd.norm.method = normMethod;
funPsdNorm.psd.norm.fact = normFact;
switch normMethod
    case 'psdNoiseTo1_thenLfSlope'
        funPsdNorm.psd.norm.fact2 = normFact2;
        funPsdNorm.psd.norm.a = a;
        funPsdNorm.psd.norm.b = b;
    otherwise
end
funPsdNorm = setNiceFieldOrder(funPsdNorm,{'vol' 'vol2vec' 'vol2vecFlag' 'vec' 'vecNorm'});

