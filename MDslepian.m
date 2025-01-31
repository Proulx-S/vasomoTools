function [lambda,v] = MDslepian(w,k,t,Fs,algo)
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
if ~exist('algo','var'); algo = ''; end
if isempty(algo); algo = 'eigs'; end % 'eigs' or 'eig'


rng('default');
n = length(t);
a = zeros(n,n);
disp('using parfor loop')
parfor i = 1:n
    a(i,:) = sin(2*pi*w*(t(i) - t))./(pi*(t(i) - t));
end
disp('parfor loop done')
a(eye(n,'logical')) = 2*w;
% a = zeros(n,n);
% a(1:n,1:n) = 2*w;
% for i = 1:n
%     j = i+1:n;
%     a(i,j) = sin(2*pi*w*(t(i) - t(j)))./(pi*(t(i) - t(j)));
%     a(j,i) = a(i,j);
% end

switch algo
    case 'eigs' % original
        disp('using eigs')
        [v,lambda,flag] = eigs(a,k,'largestreal','display',1);
        lambda = diag(lambda);
    case 'eig' % experimental
        disp('using eig')
        [v,lambda] = eig(a,'vector');
        flag = any(isnan(lambda));
end
[lambda,i] = sort(lambda,'descend');
lambda = lambda(1:k);
v = v(:,i(1:k));

if flag
    warning(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' newline ...
             '!! MDslepian tapers did not converge !!' newline ...
             '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!']);
end

% lambda = diag(lambda);
% [V,D,~] = eig(a,'vector');
% [D,b] = sort(D,'descend');
% V = V(:,b);
% V = V(:,1:k);
% D = D(1:k);
% figure('WindowStyle','docked')
% plot(lambda)
% hold on
% plot(D,'--r')
% for i = 1:k
%     figure('WindowStyle','docked')
%     plot(V(:,i)); hold on
%     plot(v(:,i))
% end
% 



for i = 1:2:k
    if mean(real(v(:,i))) < 0, v(:,i) = -v(:,i); end
end
for i = 2:2:k-1
    if real(v(2,i) - v(1,i)) < 0, v(:,i) = -v(:,i); end
end

%%%%%%%%%%%%%%%
% To match the scale of the same tapers obtain with Chronux
v = v*sqrt(Fs);
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
% To match eigenvalues of the same tapers obtain with Matlab's dpss.m
lambda = (lambda/Fs)';
%%%%%%%%%%%%%%%
end