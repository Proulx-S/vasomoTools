function [volPsdSim, volTsSim] = simPsd2(onsetList,durList,volTs,info,dt,plotFlag)
volTs = vec2vol(volTs); volTs.vol = [];
if isfield(volTs,'resp'); volTs = rmfield(volTs,'resp'); end
if isfield(volTs,'dsgn'); volTs = rmfield(volTs,'dsgn'); end
if ~exist('plotFlag','var'); plotFlag = []; end
if isempty(plotFlag);        plotFlag = 0; end

ndummy  = volTs.dummy;
nframes = volTs.nframes + ndummy;
tr = volTs.tr/1000;
if ~exist('dt','var'); dt = volTs.tr/1000; end

SNR = 10;
rndSeed = 0;
[t, d, c] = simTs(onsetList,durList,nframes,tr,dt,SNR,plotFlag,rndSeed);
t(1:ndummy) = [];
d(1:ndummy) = [];
c(1:ndummy) = [];
rndSeed = 1;
[tShift, dShift, cShift] = simTs(onsetList+2.5,durList,nframes,tr,dt,SNR,plotFlag,rndSeed);
tShift(1:ndummy) = [];
dShift(1:ndummy) = [];
cShift(1:ndummy) = [];

% T = ceil((onsetList(end)+mean(diff(onsetList)))./tr);
% t = 0:tr:(T-1)*tr; % 5 seconds of data (time)
% f = 0.05;
% rng(rndSeed)
% phi = 0;
% c = sin(2*pi*f*t+phi).*SNR + randn(size(t));
% phi = 2*pi/16;
% cShift = sin(2*pi*f*t+phi).*SNR + randn(size(t));


volTsSim = vec2vol(volTs);
volTsSim.vol = permute(cat(1,c,cShift),[1 3 4 2]);

nVoxNull = 10;
volTsSim.vol = cat(1,volTsSim.vol,randn([nVoxNull 1 1 size(volTsSim.vol,4)]));
volTsSim.dsgn.onsets = onsetList;
volTsSim.dsgn.ondurs = durList;

volTsSim.t = t';
volTsSim.nframes = length(t);

info.skipSvd = 0;
info.onsetList = onsetList;
info.ondurList = durList;
volPsdSim = volPsdFullMt2([],info,volTsSim);



