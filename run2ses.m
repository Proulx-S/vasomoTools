function mri = run2ses(mri)
if length(mri)<=1; warning('no run to catenate'); return; end
if ~isfield(mri,'volInfo'); [mri.volInfo] = deal(strjoin({'X' 'Y' 'Z' 'time' 'taper' 'run'},' x ')); end
if ~isfield(mri,'vecInfo'); [mri.vecInfo] = deal(strjoin({'time' 'vox' 'taper' 'run'},' x ')); end

%% Make sure all runs have the same mask
if isfield(mri,'vol2vec') && ~isempty(mri(1).vol2vec)
    tmp = permute(cat(4,mri.vol2vec),[4 1 2 3]);
    if diff(sum(tmp(:,:),2))
        tmp = permute(all(tmp,1),[2 3 4 1]);
        if isempty(mri(1).vec)
            [mri.vol2vec] = deal(tmp);
        elseif isempty(mri(1).vol)
            mri = vol2vec(vec2vol(mri),tmp,1);
        end
    end
end

%% Add time
mri = addTime(mri);

%% Concate runs
mri2 = mri(1);
if isfield(mri(1),'vol') && ~isempty(mri(1).vol)
    mri2.vol = cat(6,mri.vol); [mri.vol] = deal([]);
end
if isfield(mri(1),'vec') && ~isempty(mri(1).vec)
    mri2.vec = cat(4,mri.vec); [mri.vec] = deal([]);
end
if isfield(mri,'t')
    mri2.t = cat(4,mri.t);
end
if isfield(mri,'nruns') && ~isempty(mri(1).nruns)
    mri2.nruns = sum([mri.nruns]);
else
    mri2.nruns = length(mri);
end
if isfield(mri,'scale') && ~isempty(mri.scale)
    mri2.scale = cat(6,mri.scale);
else
    mri2.scale = [];
end
if isfield(mri,'shift') && ~isempty(mri.shift)
    mri2.shift = cat(6,mri.shift);
else
    mri2.shift = [];
end
if isfield(mri(1),'acqTime') && isdatetime(mri(1).acqTime)
    mri2.acqTime = cat(6,mri.acqTime);
end




if isfield(mri2,'vol2vec') && ~isempty(mri2.vol2vec) && ( (isfield(mri,'scale') && ~isempty(mri.scale)) || (isfield(mri,'shift') && ~isempty(mri.shift)) )
    dbstack; error('double-check that')
    mri2.vol2vec = all(~isnan(mri2.scale),6) && all(~isnan(mri2.shift),6);
end


if isfield(mri,'dtrnd')
    dbstack; error('double-check that')
end

if isfield(mri,'imMean')
    mri2.imMean = cat(6,mri.imMean);
elseif ~isempty(mri(1).vol)
    mri2.imMean = mean(mri2.vol,4);
% else
%     mri2.imMean = nan(size(mri2.vol2vec));
%     mri2.imMean(mri2.vol2vec) = mean(mri.vec,1);
end

mri = mri2; clear mri2
