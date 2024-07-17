function [volPsdSim, volTsSim] = simPsd(onsetList,durList,volTs,info,dt,plotFlag)
volTs = vec2vol(volTs); volTs.vol = [];
if ~exist('plotFlag','var'); plotFlag = []; end
if isempty(plotFlag);        plotFlag = 0; end


nframes = volTs.nframes;
tr = volTs.tr/1000;
if ~exist('dt','var'); dt = volTs.tr/1000; end

SNR = 10;
rndSeed = 0;
[t, d, c] = simTs(onsetList,durList,nframes,tr,dt,SNR,plotFlag,rndSeed);
rndSeed = 1;
[tShift, dShift, cShift] = simTs(onsetList+2.5,durList,nframes,tr,dt,SNR,plotFlag,rndSeed);


volTsSim = vec2vol(volTs);
% volTsSim.t = t';
volTsSim.vol = permute(cat(1,c,cShift),[1 3 4 2]);

nVoxNull = 10;
volTsSim.vol = cat(1,volTsSim.vol,randn([nVoxNull 1 1 size(volTsSim.vol,4)]));
volTsSim.dsgn.onsetList = onsetList;
volTsSim.dsgn.durList = durList;

info.skipSvd = 0;
info.onsetList = onsetList;
info.ondurList = durList;
volPsdSim = volPsdFullMt2([],info,volTsSim);



