function funTs = normPSD3(funTs,noiseRange,modeToRemove)
if ~isfield(funTs,'imMean') || isempty(cat(6,funTs.imMean))
    for runInd = 1:length(funTs)
        funTs(runInd).imMean = mean(funTs(runInd).vol,4);
    end
end
funTs = vol2vec(funTs);
if ~exist("noiseRange",'var') || isempty(noiseRange); noiseRangeFlag = 0; else; noiseRangeflag = 1; end
if ~exist("modeToRemove",'var') || isempty(modeToRemove); modeToRemove = 1:10; end
param.tapers = [1 1]; param.Fs = 1/(funTs(1).tr/1000);
[~,f] = mtspectrumc(funTs(1).vec(:,1),param);
noiseRange(isinf(noiseRange)) = f(end);


%% Filter out noise using svd
disp('Normalization: computing')
for runInd = 1:length(funTs)
    disp(['Run' num2str(runInd) '/' num2str(length(funTs))])
    %%% Demean
    funTsX = funTs(runInd);
    funTsX.vec = funTsX.vec - mean(funTsX.vec,1);
    
    %%% Estimate physiological noise filter
    svdProulx = runMTsvdProulx(funTsX,noiseRange); clear funTsX
    % modeTs = reconMTsvdProulx_modeTs(svdProulx,funTs(runInd));
    % funPsd = runPSD4(funTs(runInd));
    % modePsd = reconMTsvdProulx_modePsd(svdProulx,funPsd);
    % for modeOrder = 1:10
    %     viewMTsvdProulx(svdProulx,modePsd,modeTs,modeOrder,funPsd,funTs(runInd));
    % end
    
    %%% Reconstruct data without physiological noise
    recTs = reconMTsvdProulx_spaceXtime(svdProulx,funTs(runInd),modeToRemove,[]);

    %%% Estimate noise floor from reconstructed data
    skipSVD = 1;
    verbose = 0;
    recPsd = runPSD6(recTs,[],[],[],[],skipSVD,[],verbose);
    % viewPSD3(recPsd,{noiseRange + [1 -1].*range(noiseRange).*0.001})
    noiseInd = logical(recPsd.psd.f>=noiseRange(1) & recPsd.psd.f<=noiseRange(2));
    noiseFloor = mean(conj(recPsd.vec(noiseInd,:,:,:)).*recPsd.vec(noiseInd,:,:,:),1);

    %%% Scale timeseries according to noise floor
    funTs(runInd).vec = funTs(runInd).vec ./ sqrt(noiseFloor);
    funTs(runInd).scale = nan(size(funTs(runInd).vol2vec));
    funTs(runInd).scale(funTs(runInd).vol2vec) = sqrt(noiseFloor);

    %%% Remove global mean
    shift = mean(funTs(runInd).vec(:));
    funTs(runInd).vec = funTs(runInd).vec - shift;
    funTs(runInd).shift = nan(size(funTs(runInd).vol2vec));
    funTs(runInd).shift(funTs(runInd).vol2vec) = shift;
end
disp('Normalization: done')
