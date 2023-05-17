function funPsdNorm = normPSD(funPsd,normMethod)
% normMethod: 'psdNoise_to1', 'psdAv_to1', 'psdNoise_toImAv', 'psdAv_toImAv' or 'none'

funPsd = vol2vec(funPsd);
%% Perform spatial normalization
switch normMethod
    case 'psdNoise_toImAv'
        noiseInd = funPsd.psd.f>4 & funPsd.psd.f<funPsd.psd.f(round(0.95*end));
        normFact = exp( mean(log(funPsd.vec(noiseInd,:)),1) );
        normFact = normFact .* exp(mean(log(normFact),2));
    case 'psdAv_toImAv'
        normFact = exp( mean(log(funPsd.vec),1) );
        normFact = normFact .* exp(mean(log(normFact),2));
    case 'psdNoise_to1'
        noiseInd = funPsd.psd.f>4 & funPsd.psd.f<funPsd.psd.f(round(0.95*end));
        normFact = exp( mean(log(funPsd.vec(noiseInd,:)),1) );
    case 'psdAv_to1'
        normFact = exp( mean(log(funPsd.vec),1) );
    case 'none'
        sz = size(funPsd.vec);
        normFact = ones([1 sz(2)]);
    otherwise
        error('')
end

%%% Apply normalization
funPsdNorm = funPsd;
funPsdNorm.vec = funPsd.vec./normFact;
funPsdNorm.psd.norm.method = normMethod;
funPsdNorm.psd.norm.fact = normFact;
