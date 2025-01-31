function [pv21,pv21cent,pv21offCent,c] = pvCalc(a,r1,r2)
% Calculate tha maximum possible change in partial volume given voxel size and
% vessel radius before and after dilation.
% a  -> half voxel size
% r1 -> vessel radius before
% r2 -> vessel radius after

[~,b]   = sort([r1(1) r2(1)]);
if all(b==[2 1])
    r = r1;
    r1 = r2;
    r2 = r;
    clear r
end

%% Voxel centered on vessel
pv1 = calcCent(a,r1);
pv2 = calcCent(a,r2);
pv21cent = pv2-pv1;

%% Voxel overlapping the edge of the vessel (would be entirely intravascular if centered)
ind = r1>=a; %.*sqrt(2);
pv21offCent = nan(size(ind));
sqr = (2*a)^2;
theta  = asin(a./r1(ind));
smlArc = (theta / (2*pi))  .*  pi.*r1(ind).^2;
bigArc = (theta / (2*pi))  .*  pi.*r2(ind).^2;
arcSct = 2*bigArc - 2*smlArc;
arcSctTip = (theta ./ (2*pi))  .*  pi.*(r2(ind)-r1(ind)).^2;
arcSctOvr = arcSct - 2*arcSctTip;
pv21offCent(ind) = arcSctOvr / sqr;

c = nan(size(ind));
c(ind) = a+r1(ind) - (r1(ind) - sqrt(r1(ind).^2-a.^2));
c(~ind) = a+r1(~ind);
% pv1 = calcOffCent(a,r(1),c);
% pv2 = calcOffCent(a,r(2),c);
% pv2-pv1

%% Pick the max of the 2
pv21 = max(pv21cent,pv21offCent);
if all(b == [2 1])
    pv21 = -pv21;
end

% figure('WindowStyle','docked')
% plot(r1,pv21cent)
% hold on
% plot(r1,pv21offCent)


function pv = calcCent(a,r)
% https://math.stackexchange.com/questions/1450961/overlapping-area-between-a-circle-and-a-square
pv = nan(size(r));

ind1 = r<=a;
if any(ind1)
    % Unresolved vessel (diamter smaller than voxel size).
    % Maximum overlap will always be for the vessel
    % centered in the voxel.
    crc = pi.*r(ind1).^2;
    sq  = (2.*a).^2;
    pv(ind1) = crc ./ sq;
end
    
ind2 = r>=a*sqrt(2);
% Voxel that can be fully intravascular. Maximum will be for a vessel
% right next to a voxel and overlapping it half way.
if any(ind2)
    pv(ind2) = 1;
end

% Partially resolved vessel (no purely intravascular voxel)
ind3 = ~ind1 & ~ind2;
if any(ind3)
    theta = acos(a./r(ind3));
    crcOver = (pi - 4*theta).*r(ind3).^2 + 4.*a.*r(ind3).*sin(theta);
    sq = (2*a)^2;
    pv(ind3) = crcOver ./ sq;
end


function pv = calcOffCent(a,r,c,p)
if ~exist('p','var'); p = []; end
if isempty(p); p = 10000; end
x  = linspace(-a,a,p);
im = 0;
X = x-c;
    for ii = 1:length(x)
        Y = x(ii);
        [~,rho] = cart2pol(X,Y);
        im = im + double(rho<r);
    end
pv = sum(im)/p^2;

