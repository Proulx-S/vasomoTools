function rec = reconMTsvdProulx_spaceXtime(svdStruct,funTs,modeOrderOut,modeOrderIn)
if ~exist('modeOrderOut','var') || isempty(modeOrderOut); outFlag = 0; else; outFlag = 1; end
if ~exist('modeOrderIn','var') || isempty(modeOrderIn); inFlag = 0; else; inFlag = 1; end
if all(~[outFlag inFlag]) || all([outFlag inFlag])
    error('specify one of modeOrderOut OR modeOrderIn')
elseif outFlag && ~inFlag
    modeOrder = true(size(svdStruct.sv));
    modeOrder(modeOrderOut) = false;
elseif ~outFlag && inFlag
    modeOrder = false(size(svdStruct.sv));
    modeOrder(modeOrderIn) = true;
end
if ~isfield(funTs,'imMean') || isempty(funTs.imMean)
    funTs.imMean = mean(funTs.vol(:,:,:,:),4);
end

%% Reconstruct data in taper domain
sp = squeeze(svdStruct.sp(:,:,modeOrder));
sv = diag(squeeze(svdStruct.sv(:,:,modeOrder)));
fm = squeeze(svdStruct.fm(:,:,modeOrder));
A = sp*sv*fm';
%% Project back to time domian
tp = permute(svdStruct.MTS.proj,[2 1 3]);
tsRec = real( tp*A' );

% funTs = vol2vec(funTs,svdStruct.vol2vec);
% funTs = dtrnd4psd(funTs);
% randInd = randperm(numel(data),1000);
% figure('WindowStyle','docked');
% scatter(tsRec(:),funTs.vec(:))
% grid on
% drawnow


% sv = svdStruct.sv(:,:,modeOrder); % spatial sv (space x 1 x mode)
% sp = svdStruct.sp(:,:,modeOrder); % spatial sv (space x 1 x mode)
% fm = svdStruct.fm(:,:,modeOrder); % taper sv (taper x 1 x mode)
% tp = svdStruct.MTS.proj; % tapers (taper x time)
% f = svdStruct.MTS.bandFreq; % frequencies (1 x freq)
% fInd = svdStruct.MTS.bandInd; % frequency indices (taper x 1)
% tsRec = permute(fm,[3 1 2])*tp; % WARNING! the ' operator would take the conjugate
% tsRec = tsRec .* squeeze(sv); % scale according to singular value
% tsRec = squeeze(sp)*tsRec; % project back to space time domain

%% Output
rec = funTs;
rec.vol2vec = svdStruct.vol2vec;
rec.vol = [];
rec.vec = tsRec;

