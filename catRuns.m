function mriOut = catRuns(mriIn)

mriOut = mriIn(1);
if length(mriIn)==1; return; end

if ~isfield(mriOut,'vec') || isempty(mriOut.vec)
    mriOut.vol = cat(6,mriIn.vol);
    mriOut.volMean = cat(6,mriIn.volMean);
else
    mriOut.vec = cat(4,mriIn.vec);
    mriOut.vecMean = cat(4,mriIn.vecMean);
end
mriOut.nruns = length(mriIn);
mriOut.imMean = cat(6,mriIn.imMean);




