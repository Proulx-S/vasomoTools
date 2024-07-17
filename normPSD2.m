function mriOut = normPSD2(funPsd,mriOut,noiseRange)
% normMethod: 'psdNoise_to1', 'psdAv_to1', 'psdNoise_toImAv', 'psdAv_toImAv' or 'none'
funPsd = vol2vec(funPsd);
if ~exist("noiseRange",'var') || isempty(noiseRange); noiseRangeFlag = 0; else; noiseRangeflag = 1; end
noiseRange(isinf(noiseRange)) = funPsd.psd.f(end);

%% Compute noise floor
noiseInd = logical(funPsd.psd.f>=noiseRange(1) & funPsd.psd.f<=noiseRange(2));
noiseFloor = mean(conj(funPsd.vec(noiseInd,:,:,:)).*funPsd.vec(noiseInd,:,:,:),1);
funPsd.psd.noiseFloor = sqrt(noiseFloor);

if ~exist('mriOut','var') || isempty(mriOut)
    %% Apply normalization to funPsd
    funPsd.normFac = funPsd.psd.noiseFloor;
    funPsd.vec = funPsd.vec./funPsd.normFac;
    mriOut = funPsd;
else
    %% Apply normalization to mriOut
    if ~isfield(mriOut,'imMean') || isempty(mriOut.imMean)
        mriOut.imMean = mean(mriOut.vol,4);
    end
    mriOut.normFac = nan(size(funPsd.vol2vec));
    mriOut.normFac(funPsd.vol2vec) = funPsd.psd.noiseFloor;
    mriOut.vol = mriOut.vol./mriOut.normFac;
end

return

if noiseRangeFlag && ~physRangeFlag
    noiseInd = logical(funPsd.psd.f>noiseRange(1) & funPsd.psd.f<noiseRange(2));
elseif ~noiseRangeFlag && physRangeFlag
    for i = 1:size(physRange,1)
    % i=2;
        curRange = physRange(i,:);
        while curRange(end,1) < funPsd.psd.f(end)
            curRange = cat(1 , curRange , curRange(end,:)+curRange(1,:) );
        end
        curRange(curRange(:,1) >= funPsd.psd.f(end),:) = [];



        delta = diff(physRange(i,:));
        cent = physRange(i,1)+delta/2;
    end
    
elseif noiseRangeFlag && physRangeFlag
    error('specify one of noiseRange OR physRange')
else
    noiseInd = logical(funPsd.psd.f>4 & funPsd.psd.f<funPsd.psd.f(round(0.95*end)));
end
% if ~exist('noiseRange','var') || isempty(noiseRange)
%     noiseInd = logical(funPsd.psd.f>4 & funPsd.psd.f<funPsd.psd.f(round(0.95*end)));
% else
%     noiseInd = logical(funPsd.psd.f>noiseRange(1) & funPsd.psd.f<noiseRange(2));
% end
funPsd.psd.noiseFloor = mean(funPsd.vec(noiseInd,:,:,:),1);
funPsd.psd.noiseFloorThresh = prctile(funPsd.vec(noiseInd,:,:,:),95);
funPsd.psd.noiseFloorInd = noiseInd';
funPsd.psd.aboveNoiseInd = funPsd.vec>funPsd.psd.noiseFloorThresh;

%% Perform spatial normalization
switch normMethod
    case 'lowFreqTo1'
        dbstack; error('double-check that');
        lowFreqInd = funPsd.psd.f>0.01 & funPsd.psd.f<0.04;
        normFact = mean(funPsd.vec(lowFreqInd,:),1);
    case 'cardiacTo1'
        dbstack; error('double-check that');
        cardiacInd = funPsd.psd.f>0.7927 & funPsd.psd.f<1.091;
        normFact = mean(funPsd.vec(cardiacInd,:),1);
    case 'psdNoise_toImAv'
        dbstack; error('double-check that for linear averaging')
        normFact = exp( mean(log(funPsd.vec(noiseInd,:)),1) );
        normFact = normFact .* exp(mean(log(normFact),2));
    case 'psdAv_toImAv'
        dbstack; error('double-check that for linear averaging')
        normFact = exp( mean(log(funPsd.vec),1) );
        normFact = normFact .* exp(mean(log(normFact),2));
    case {'psdNoise_to1' 'psdNoiseTo1_thenLfSlope'}
        normFact = mean(funPsd.vec(noiseInd,:,:,:),1);
%         normFact = sqrt(normFact);
    case 'psdAv_to1'
        dbstack; error('double-check that');
        normFact = exp( mean(log(funPsd.vec),1) );
    case 'none'
        dbstack; error('double-check that');
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
        dbstack; error('double-check that');
        f = funPsdNorm.psd.f';
        psd = funPsdNorm.vec;

        %%% Define the range within which to fit the low-frequency slope
%         lfSlopeRange = [
%             0    0.01
% %             0.17 0.2
%             0.13 0.2
%             ];
        lfSlopeRange = [
%             0.01 0.1
%             0    0.01
            0.013 0.019
%             0.06 0.08
            0.13 0.15
%             0.17 0.2
            ];


        lfSlopeRange(lfSlopeRange==0) = f(find(f==0)+1); % make sure f=0 is not included
        ind = false(size(f));
        for i = 1:size(lfSlopeRange,1)
            ind(f>=lfSlopeRange(i,1) & f<=lfSlopeRange(i,2)) = true;
        end

        %%% Fit a linear slope in loglog space, for each voxel
        lfSlp = nan(size(psd));
        a = nan(1,size(psd,2));
        b = nan(1,size(psd,2));
        parfor ii = 1:size(psd,2)
            [lfSlp(:,ii),a(ii),b(ii)] = logloglin(f,psd(:,ii),ind,noiseInd);
        end

        %%% Define a new frequency- and voxel-specific normalization factor
        normFact2 = lfSlp;
        % set the slope at the noise floor as soon as it intercepts it
        normFact2(a<0 & normFact2<1) = 1;
        normFact2(a>0 & normFact2>1) = 1;
        % set the slope at the noise floor if it never intercepts it
        normFact2( : , all(normFact2>1,1) | all(normFact2<1,1) ) = 1;
        % set f=0 to 0 as it is ill-defined
        normFact2(f==0,:) = 0;

        if exist('mask','var') && ~isempty(mask)
            voxInd = logical(mask(funPsdNorm.vol2vec));
        else
            voxInd = true(1,size(psd,2));
        end
        figure('WindowStyle','docked');
        plot(f,mean(psd(:,voxInd),2),'k'); hold on
        tmp = mean(normFact2(:,voxInd),2);
        plot(f,tmp,'--r')
        tmp(~ind) = nan;
        plot(f,tmp,'r','linewidt',3)
        psd2 = psd./normFact2;
        hPlot = plot(f,mean(psd2(:,voxInd),2));
        plot(f,ones(size(f)),':k')
        tmp = mean(psd2,2);
        plot(f([2 end]),tmp([2 end]),':','Color',hPlot.Color)
        ax = gca;
        ax.YScale = 'log';
        ax.XScale = 'log';
        
        %%% Apply the normalization
        funPsdNorm.vec = funPsdNorm.vec./normFact2;
    otherwise
end

%%% Output
volNorm = nan([size(normFact,[3 4]) size(funPsdNorm.vol2vec)]);
volNorm(:,:,funPsdNorm.vol2vec) = permute(normFact,[3 4 2 1]);
funPsdNorm.volNorm = permute(volNorm,[3 4 5 1 2]);
funPsdNorm.psd.norm.method = normMethod;
funPsdNorm.psd.norm.fact = normFact;
switch normMethod
    case 'psdNoiseTo1_thenLfSlope'
        dbstack; error('double-check that');
        funPsdNorm.psd.norm.fact2 = normFact2;
        funPsdNorm.psd.norm.a = a;
        funPsdNorm.psd.norm.b = b;
    otherwise
end
funPsdNorm = setNiceFieldOrder(funPsdNorm,{'vol' 'vol2vec' 'vol2vecFlag' 'vec' 'volNorm'});

