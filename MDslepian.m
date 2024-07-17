function [lambda,u] = MDslepian(w,k,t,Fs)
% Extracted from MDmwps
% (https://www.mathworks.com/matlabcentral/fileexchange/71909-mdmwps?s_tid=srchtitle_site_search_1_MDmwps)
% and modified to match the scale of the same tapers computed with Chronux
% and the eigenvalues of the same tapers computed with Matlab's dpss.m

%computes generalized slepian function for 1D missing data problem
%input variables
%     w = analysis half bandwidth
%     k = number of eigenvalue/vectors to compute (must be <=2*bw*length(x), 2*bw*length(x)-1 default)
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

%%%%%%%%%%%%%%%
% To match the scale of the same tapers obtain with Chronux
u = u*sqrt(Fs);
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
% To match eigenvalues of the same tapers obtain with Matlab's dpss.m
lambda = (lambda/Fs)';
%%%%%%%%%%%%%%%
end