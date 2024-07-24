function [volPsdSim, volTsSim] = simPsd3(dsgn,info,plotFlag)
if ~exist('plotFlag','var'); plotFlag = []; end
if isempty(plotFlag);        plotFlag = 0; end
if ~isfield(dsgn,'f0');       dsgn.f0 = []; end
SNR = 5; %10;
nVoxNull = 12; %10;


if ~isempty(dsgn.onsets) && all((size(dsgn.onsets)==1)==[1 0])
    dbstack; error('onsets must be a column vector')
end


% ndummy  = dsgn.dummy;
nframes = dsgn.nframes;
dt = dsgn.dt; % [sec]

if isempty(dsgn.f0)
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


else
    tS = 0; tE = (dsgn.nframes+dsgn.dummy-1)*dsgn.dt;
    t = tS:dt:tE;
    theta  = 0;
    amp    = 1;
    c      = real(   amp .* exp(1i.* (2*pi*dsgn.f0.*t - theta)  )   )  +  1;
    theta  = 1/8 .* (2*pi);
    amp    = 1;
    cShift = real(   amp .* exp(1i.* (2*pi*dsgn.f0.*t - theta)  )   )  +  1;
    % figure('WindowStyle','docked');
    % plot(t,real(c),'.-'); hold on
    % plot(t,real(cShift),'.-r');

    %% Add noise and null voxels
    volTsSim.vol = permute(cat(1,c,cShift),[1 3 4 2]);
    volTsSim.vol = cat(1,volTsSim.vol,zeros([nVoxNull 1 1 size(volTsSim.vol,4)])   +  1);
    if SNR~=inf
        rng(0)
        volTsSim.vol = volTsSim.vol.*SNR + randn(size(volTsSim.vol)).*std(c);
    end
    
    volTsSim.t = t';
    volTsSim.vol(:,:,:,1:dsgn.dummy) = [];
    volTsSim.t(1:dsgn.dummy) = [];

    volTsSim.dsgn = dsgn;
    volTsSim.nframes = dsgn.nframes;
    volTsSim.dummy = dsgn.dummy;
    volTsSim.nvoxels = prod(size(volTsSim.vol,1:3));
    volTsSim.tr = dt*1000;

    volPsdSim = [];
end



