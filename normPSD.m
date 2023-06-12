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


%%% Output
funPsdNorm.vecNorm = normFact;
funPsdNorm.psd.norm.method = normMethod;
funPsdNorm.psd.norm.fact = normFact;
funPsdNorm = setNiceFieldOrder(funPsdNorm,{'vol' 'vol2vec' 'vol2vecFlag' 'vec' 'vecNorm'});

