function [volPsdSim, volTsSim] = simPsd3(dsgn,info,plotFlag)
% volTs = vec2vol(volTs); volTs.vol = [];
% if isfield(volTs,'resp'); volTs = rmfield(volTs,'resp'); end
% if isfield(volTs,'dsgn'); volTs = rmfield(volTs,'dsgn'); end
if ~exist('plotFlag','var'); plotFlag = []; end
if isempty(plotFlag);        plotFlag = 0; end
SNR = inf; %10;
nVoxNull = 0; %10;



% ndummy  = dsgn.dummy;
nframes = dsgn.nframes;
dt = dsgn.dt;

onsets = (0:1/dsgn.onsetFreq:nframes*dt)';
ondurs = repmat(1/dsgn.onFreq,size(onsets));
nframes = ceil((onsets(end) + mode(diff(onsets))) ./ dt);
onsets = onsets + 1/dsgn.onsetFreq/2;

rndSeed = 0;
[t, d, c] = simTs3(onsets,ondurs,nframes,dt,dt/100,SNR,plotFlag,rndSeed);
% t(1:ndummy) = [];
% d(1:ndummy) = [];
% c(1:ndummy) = [];
rndSeed = 1;
[tShift, dShift, cShift] = simTs3(onsets+1,ondurs,nframes,dt,dt,SNR,plotFlag,rndSeed);
% tShift(1:ndummy) = [];
% dShift(1:ndummy) = [];
% cShift(1:ndummy) = [];


% T = ceil((onsetList(end)+mean(diff(onsetList)))./tr);
% t = 0:tr:(T-1)*tr; % 5 seconds of data (time)
% f = 0.05;
% rng(rndSeed)
% phi = 0;
% c = sin(2*pi*f*t+phi).*SNR + randn(size(t));
% phi = 2*pi/16;
% cShift = sin(2*pi*f*t+phi).*SNR + randn(size(t));





% volTsSim = vec2vol(volTs);


volTsSim.vol = permute(cat(1,c,[]),[1 3 4 2]);
% volTsSim.vol = permute(cat(1,c,cShift),[1 3 4 2]);
volTsSim.vol = cat(1,volTsSim.vol,randn([nVoxNull 1 1 size(volTsSim.vol,4)]));



volTsSim.dsgn.onsets = onsets;
volTsSim.dsgn.ondurs = ondurs;
volTsSim.t = t';
volTsSim.nframes = length(t);
volTsSim.nvoxels = prod(size(volTsSim.vol,1:3));
volTsSim.tr = dt*1000;
info.win(1) = round(dsgn.winLength/info.tr);

% volTsSim.dsgn.onsets([1 end]) = [];
% volTsSim.dsgn.ondurs([1 end]) = [];
volTsSim.dsgn.onsets(end) = [];
volTsSim.dsgn.ondurs(end) = [];

volPsdSim = volPsdFullMt2([],info,volTsSim);



