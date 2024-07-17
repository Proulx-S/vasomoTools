
function [sxx, varargout] = MDmwps(tt,xx,varargin)
% Source: https://www.mathworks.com/matlabcentral/fileexchange/71909-mdmwps?s_tid=srchtitle_site_search_1_MDmwps
%
%multi-taper power spectrum estimation with f-test and reshaping
% for time series with missing data
%
%input arguments (variable)
%   tt -- real vector of time (required)
%   xx -- real vector of data (required)
%   bw -- bandwidth of estimate, 5/length(t) default
%   k -- number of Slepian tapers, must be <=2*bw*length(x), 2*bw*length(x)-1 default
%   nz -- zero padding factor, 0 default
%   alpha -- probability level for reshaping, 1 if none, 1 default
%output arguments (variable)
%   sxx -- power spectrum vector of length length(x)/2+1 (required)
%   nu1 -- degrees-of-freedom vector for sxx of length length(x)/2+1
%   ft -- f-test vector of length length(x)/2+1
%   il -- frequency indices of spectral lines
%   plin -- spectral line power vector
%   srr -- reshaped power spectrum vector of length length(x)/2+1 if alpha<1
%   nu2 -- degrees-of-freedom vector for srr of length length(x)/2+1
%
switch nargin
    case 2
        bw = 5/length(tt);
        k = 2*bw*length(tt) - 1;
        nz = 0;
        alpha = 1;
    case 3
        bw = varargin{1};
        k = 2*bw*length(tt) - 1;
        nz = 0;
        alpha = 1;
    case 4
        bw = varargin{1};
        k = varargin{2};
        nz = 0;
        alpha = 1;
    case 5
        bw = varargin{1};
        k = varargin{2};
        nz = varargin{3};
        alpha = 1;
    case 6
        bw = varargin{1};
        k = varargin{2};
        nz = varargin{3};
        alpha = varargin{4};
end
if alpha == 1 && nargout > 3
    switch nargout
        case 4
            varargout{4} = [];
        case 5
            varargout{4} = [];
            varargout{5} = [];
        case 6
            varargout{4} = [];
            varargout{5} = [];
            varargout{6} = [];
        case 7
            varargout{4} = [];
            varargout{5} = [];
            varargout{6} = [];
            varargout{7} = [];
    end
end
[n,p] = size(tt);
if n == 1
    t = tt';
else
    t = tt;
end
[n,p] = size(xx);
if n == 1
    x = xx';
else
    x = xx;
end
nfft = length(x);
if mod(nfft,2) ~= 0, nfft = nfft + 1; end
nfft = (nz + 1)*nfft;
nfft2 = nfft/2 + 1;
s2 = var(x,1);
e = zeros(nfft2,k,length(x));
sxx = zeros(1,nfft2);
[lambda,u] = MDslepian(bw,k,tt);

parfor i = 1:nfft2
    f = (i - 1)/nfft;
    e(i,:,:) = (u.*repmat(exp(complex(0,-2*pi*f*t)),1,k)).';
    sk = abs(squeeze(e(i,:,:))*x).^2;
	sxx(i) = fzero(@(y) Mtaper(y,lambda,sk,s2),(sk(1) + sk(2))/2);
end
sxx(2:nfft2-1) = 2*sxx(2:nfft2-1);
if nargout == 1, return, end
d = repmat(sqrt(lambda),1,nfft2).*repmat(sxx,k,1)./ ...
    (repmat(lambda,1,nfft2).*repmat(sxx,k,1) + ...
    s2*(ones(k,nfft2) - repmat(lambda,1,nfft2)));
nu1 = 2*lambda'*d.^2;
varargout{1} = nu1;
if nargout == 2, return, end
dpsw0 = squeeze(sum(e(1,:,:),3));
ak = zeros(nfft2,k);
for i = 1:nfft2
    ak(i,:) = squeeze(e(i,:,:))*x;
end
mu = ak(1:nfft2,:)*dpsw0.'/sum(dpsw0.^2);
num = (nu1' - 2).*abs(mu).^2*sum(dpsw0.^2);
denom = 2*sum(abs(ak(1:nfft2,:) - mu*dpsw0).^2,2);
ft = (num./denom)';
varargout{2} = ft;
if nargout == 3, return, end
if alpha == 1, return, end
il1 = find(ft >= finv(alpha,2,nu1-2));
varargout{3} = il1;
dpsw = squeeze(sum(e,3));
dpsw = fftshift([dpsw.' conj(dpsw(nfft2-1:-1:2,:)).'].',1);
zk = fftshift([ak.' conj(ak(nfft2-1:-1:2,:)).'].',1);
if length(il1) ~=0
    n = (nz + 1)*round(bw*length(x)) + 1;
    ii = 1;
    jj = 1;
    il = [];
    i = 2:length(il1);
    i1 = [find(il1(i) ~= il1(i-1) + 1) length(il1)];
    j = -n + 1:n - 1;
    for i = 1:length(i1)
        if i1(i) == ii
            m = nfft2 + il1(ii) + j - 1;
            mm = m > length(zk);
            m(mm) = m(mm) - length(zk);
            zk(m,1:k) = zk(m,1:k) - mu(il1(ii))*dpsw(nfft2+j,1:k);
            il(jj) = il1(ii);
            plin(jj) = mean(sum(abs(mu(il1(ii))*dpsw(nfft2+j,1:k)).^2));
            ii = ii + 1;
            jj = jj + 1;
        else
            i2 = find(ft == max(ft(il1(ii):il1(i1(i)))));
            m = nfft2 + i2 + j - 1;
            mm = m > length(zk);
            m(mm) = m(mm) - length(zk);
            zk(m,1:k) = zk(m,1:k) - mu(i2)*dpsw(nfft2+j,1:k);
            il(jj) = i2;
            plin(jj) = mean(sum(abs(mu(i2)*dpsw(nfft2+j,1:k)).^2));
            ii = i1(i) + 1;
            jj = jj + 1;
            end
        end
        varargout{3} = il;
end
if nargout == 5, return, end
zk = ifftshift(zk,1);
s2 = mean(abs(zk(1,:).^2 + 2*sum(abs(zk(2:nfft2-1,:)).^2) + ...
              abs(zk(nfft2,:).^2)))/nfft;
sk = abs(zk).^2;
srr = zeros(1,nfft2);
parfor i = 1:nfft2
    srr(i) = fzero(@(y) Mtaper(y,lambda,sk(i,:)',s2), ...
             (sk(i,1) + sk(i,2))/2);
end
srr(2:nfft2-1) = 2*srr(2:nfft2-1);
if length(il1) ~= 0
    if il(1) ~= 1
        plin = 2*plin;
    else
        plin(2:length(plin)) = 2*plin(2:length(plin));
    end
end
if length(il1) ~= 0
    varargout{4} = plin;
else
    varargout{4} = [];
end
if nargout == 5, return, end
varargout{5} = srr;
if nargout == 6, return, end
d = repmat(sqrt(lambda),1,nfft2).*repmat(srr,k,1)./ ...
    (repmat(lambda,1,nfft2).*repmat(srr,k,1) + ...
    s2*(ones(k,nfft2) - repmat(lambda,1,nfft2)));
nu2 = 2*lambda'*d.^2;
varargout{6} = nu2;
end
function Result = Mtaper(x,v,sk,s2)
    Result = sum(v.^2.*(x - sk)./(v*x + s2*(ones(size(v)) - v)).^2);
end
function [lambda u] = MDslepian(w,k,t)
%computes generalized slepian function for 1D missing data problem
%input variables
%     w = analysis half bandwidth
%     k = number of eigenvalue/vectors to compute
%     t = time vector
%output variables
%     lambda = eigenvalues
%     u = eigenvectors
% rng(2147483647);
rng('default');
n = length(t);
sigma = 'largestreal';
a = zeros(n,n);
a(1:n,1:n) = 2*w;
for i = 1:n
    j = i+1:n;
    a(i,j) = sin(2*pi*w*(t(i) - t(j)))./(pi*(t(i) - t(j)));
    a(j,i) = a(i,j);
end
[v,lambda] = eigs(a,k,sigma);
lambda = diag(lambda);
[lambda,i] = sort(lambda,'descend');
u = v(:,i);
for i = 1:2:k
    if mean(real(u(:,i))) < 0, u(:,i) = -u(:,i); end
end
for i = 2:2:k-1
    if real(u(2,i) - u(1,i)) < 0, u(:,i) = -u(:,i); end
end
end
