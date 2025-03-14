function flobsFit(t,y,flobsBasis)

% Use default parameters if not provided
if ~exist('flobsBasis','var') || isempty(flobsBasis)
    flobsBasis.param = doubleGamma;
    flobsBasis.amp1 = 1;
    flobsBasis.peak1 = 1:1:10;
    flobsBasis.fwhm1 = 1:1:15;
    flobsBasis.amp2 = -2:0.2:0;
    flobsBasis.peak2 = flobsBasis.peak1 + 5;
    flobsBasis.fwhm2 = flobsBasis.fwhm1 + 5;
end

% Generate flobs basis functions if not provided
if isstruct(flobsBasis)
    flobsBasis.X = [];
    % flobsBasis.X = zeros(length(t),length(flobsBasis.amp1)*length(flobsBasis.peak1)*length(flobsBasis.fwhm1)*length(flobsBasis.amp2)*length(flobsBasis.peak2)*length(flobsBasis.fwhm2));
    for i = 1:length(flobsBasis.amp1)
        for j = 1:length(flobsBasis.peak1)
            for k = 1:length(flobsBasis.fwhm1)
                for l = 1:length(flobsBasis.amp2)
                    for m = 1:length(flobsBasis.peak2)
                        for n = 1:length(flobsBasis.fwhm2)
                            if flobsBasis.amp1(i) == -flobsBasis.amp2(l) && ...
                                flobsBasis.peak1(j) == flobsBasis.peak2(m) && ...
                                flobsBasis.fwhm1(k) == flobsBasis.fwhm2(n)
                                continue;
                            else
                                flobsBasis.X(:,end+1) = doubleGamma(t,[flobsBasis.amp1(i) flobsBasis.peak1(j) flobsBasis.fwhm1(k) flobsBasis.amp2(l) flobsBasis.peak2(m) flobsBasis.fwhm2(n)]);
                            end
                            % flobsBasis.X(:,i*j*k*l*m*n) = doubleGamma(t,[flobsBasis.amp1(i) flobsBasis.peak1(j) flobsBasis.fwhm1(k) flobsBasis.amp2(l) flobsBasis.peak2(m) flobsBasis.fwhm2(n)]);
                            % if any(isnan(flobsBasis.X(:,i*j*k*l*m*n)))
                            %     keyboard;
                            % end
                            % if all(flobsBasis.X(:,i*j*k*l*m*n)==0)
                            %     [flobsBasis.amp1(i) flobsBasis.peak1(j) flobsBasis.fwhm1(k) flobsBasis.amp2(l) flobsBasis.peak2(m) flobsBasis.fwhm2(n)]
                            %     keyboard;
                            % end
                        end
                    end
                end
            end
        end
    end
end

figure;
plot(t,flobsBasis.X(:,1:end));
flobsBasis.X(:,all(flobsBasis.X==0,1)) = [];
nnz(all(flobsBasis.X==0,1))

% SVD of flobsBasis.X
[U,S,V] = svd(flobsBasis.X,"econ",'vector');
[U,S,V] = svd(flobsBasis.X(2:end,:) - mean(flobsBasis.X(2:end,:),1),"econ",'vector');
figure;
plot(t,U(:,1:5));
all(flobsBasis.X(:,i*j*k*l*m*n)==0)

flobsBasis.param = doubleGamma;
flobsBasis.amp1 = 1;
flobsBasis.peak1 = 1:0.1:10;
flobsBasis.fwhm1 = 1:0.1:15;
flobsBasis.amp2 = 0:0.1:2;
flobsBasis.peak2 = flobsBasis.peak1 + 5;
flobsBasis.fwhm2 = flobsBasis.fwhm1 + 5;


    % Normalize range of y for better behaved gradient descent
sFactor = 1./max(abs(y));
y = y.*sFactor;

%%% Stage 1: Fit only gamma1, keep gamma2 fixed at 0
disp('Stage 1: Fitting gamma1 only...');

% Parameter fixing mask for gamma1 only [AMP1 PEAK1 FWHM1 AMP2 PEAK2 FWHM2]
paramMask1 = logical([1 1 1 0 0 0]);  % Only fit gamma1 parameters

% Initial guess with no gamma2
%%% Find peak of the main lobe of the response
Idx = 1:length(t);
Idx = Idx(t <= 10);
[~, maxIdx] = max(abs(y(Idx)));
b = Idx(maxIdx);
a = y(b);
%%% Find fwhm based on the point at half the peak amplitude on the left of the peak
[~,c] = sort(abs(y(b)/2 - y),'ascend');
c = c(ismember(c,Idx));
c = c(1);
[~,c(2)] = min(abs(t - (t(b) + (t(b) - t(c)))));
c = sort(c);

figure('Position', [1200 500 400 400]);
hData = plot(t, y, 'k', 'LineWidth', 2); hold on
if y(b)<0
    hGuess1 = plot(t([b b]), [y(b) 0], 'b-', 'LineWidth', 2);
else
    hGuess1 = plot(t([b b]), [0 y(b)], 'b-', 'LineWidth', 2);
end
hGuess2 = plot(t(c), [1 1].*mean([y(c)]), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Response');
title('Fit gamma1 only');
grid on;
axis tight

peak1guess = t(b);
fwhm1guess = diff(t(c([1 2])));
amp1guess = y(b);
% amp1guess = fwhm1guess*y(b)*sqrt(2);
p0      = [  amp1guess     peak1guess     fwhm1guess 0 0 0];  % AMP2=0 to eliminate gamma2
hFit0 = plot(t, doubleGamma(t,p0), 'm--', 'LineWidth', 1);

% Set bounds for parameters [AMP1 PEAK1 FWHM1 AMP2 PEAK2 FWHM2]
lb_full = [2*amp1guess 0.5*peak1guess 0.5*fwhm1guess 0  5  0];  % Lower bounds
ub_full = [          0 1.5*peak1guess 1.5*fwhm1guess 0 20 40];  % Upper bounds (note AMP2=0)
[ub_full; p0; lb_full]

% Perform the first fit
pFit1 = fitDoubleGamma(t,y,paramMask1,p0,lb_full,ub_full,hFit0);
[ub_full; pFit1; lb_full]
if min([abs([pFit1(paramMask1)-ub_full(paramMask1) pFit1(paramMask1)-lb_full(paramMask1)])])<1e-6
    warning([newline '!!!!!!!!!!!!!!!!!!!!!!' newline 'pFit1 is at the bounds' newline '!!!!!!!!!!!!!!!!!!!!!!']);
end

% Plot first stage fit results
hFit = hFit0;
hFit0 = plot(t, doubleGamma(t,p0), 'r--', 'LineWidth', 1);
hFit.Color = 'g'; hFit.LineWidth = 2;
hResid = plot(t, y - doubleGamma(t,pFit1), 'k:', 'LineWidth', 1);
legend([hData, hGuess1, hFit, hFit0, hResid],'Data', 'Initial guess', 'Fitted HRF', 'Initial HRF', 'Residual');




if noUndershootFlag
    pFit2 = nan(size(pFit1));
    pFit3 = nan(size(pFit1));
    pFit4 = pFit1;
else


%%% Stage 2: Fit gamma2, keep gamma1 fixed at the value from stage 1
disp('Stage 2: Fitting gamma2 only...');
% Parameter fixing mask for gamma2 only [AMP1 PEAK1 FWHM1 AMP2 PEAK2 FWHM2]
paramMask2 = logical([0 0 0 1 1 1]);  % Only fit gamma2 parameters

% Initial guess
p0 = pFit1;
amp2guess = -0.3*pFit1(1);
peak2guess = pFit1(2)+7.5;
fwhm2guess = pFit1(3)*1.5;
p0(4:6) = [amp2guess peak2guess fwhm2guess];

% Set bounds for parameters [AMP1 PEAK1 FWHM1 AMP2 PEAK2 FWHM2]
ub_full = [+inf +inf +inf pFit1(1)*-2 peak2guess*2 fwhm2guess*3];  % Upper bounds (note AMP2=0)
lb_full = [-inf -inf -inf           0 peak2guess*0.5 fwhm2guess*0.5];  % Lower bounds
[ub_full; p0; lb_full]

figure('Position', [1200 500 400 400]);
hData = plot(t, y, 'k', 'LineWidth', 2); hold on
hFit0 = plot(t, doubleGamma(t,p0), 'm--', 'LineWidth', 1);
title('Fit gamma2 with fixed gamma1');
grid on;
axis tight


% pX = p0;
% pX(4:6) = 0;
% plot(t, doubleGamma(t,pX), 'm');
% pX = p0;
% pX(1:3) = 0;
% plot(t, doubleGamma(t,pX), 'm');



% Perform the second fit
pFit2 = fitDoubleGamma(t,y,paramMask2,p0,lb_full,ub_full,hFit0);
[ub_full; pFit2; lb_full]
if min([abs([pFit2(paramMask2)-ub_full(paramMask2) pFit2(paramMask2)-lb_full(paramMask2)])])<1e-6
    warning([newline '!!!!!!!!!!!!!!!!!!!!!!' newline 'pFit2 is at the bounds' newline '!!!!!!!!!!!!!!!!!!!!!!']);
end


hFit = hFit0; hFit.Color = 'g'; hFit.LineWidth = 2;
hFit0 = plot(t, doubleGamma(t,p0), 'r--', 'LineWidth', 1);
hResid = plot(t, y - doubleGamma(t,pFit2), 'k:', 'LineWidth', 1);
uistack(hFit, 'top');
legend([hData, hFit, hFit0, hResid],'Data', 'Fitted HRF', 'Initial HRF', 'Residual');


% pX = pFit2;
% pX(4:6) = 0;
% plot(t, doubleGamma(t,pX), 'm');
% pX = pFit2;
% pX(1:3) = 0;
% plot(t, doubleGamma(t,pX), 'm');

%%% Stage 3: Fit gamma1 again with fixed gamma2
disp('Stage 3: Fitting gamma1 again with fixed gamma2...');
paramMask3 = logical([1 1 1 0 0 0]);  % Only fit gamma1 parameters

% Initial guess
p0 = pFit2;

% Set bounds for parameters [AMP1 PEAK1 FWHM1 AMP2 PEAK2 FWHM2]
lb_full = [2.0*p0(1) p0(2)-5 p0(3)-5 0.5*p0(4) p0(5)-7 p0(6)-10];  % Lower bounds
ub_full = [0.5*p0(1) p0(2)+5 p0(3)+5 5.0*p0(4) p0(5)+7 p0(6)+10];  % Upper bounds
[ub_full; p0; lb_full]

figure('Position', [1200 500 400 400]);
hData = plot(t, y, 'k', 'LineWidth', 2); hold on
hFit0 = plot(t, doubleGamma(t,p0), 'm--', 'LineWidth', 1);
title('Fit gamma1 with fixed gamma2');
grid on;
axis tight

% Perform the third fit
pFit3 = fitDoubleGamma(t,y,paramMask3,p0,lb_full,ub_full,hFit0);
[ub_full; pFit3; lb_full]
if min([abs([pFit3(paramMask3)-ub_full(paramMask3) pFit3(paramMask3)-lb_full(paramMask3)])])<1e-6
    warning([newline '!!!!!!!!!!!!!!!!!!!!!!' newline 'pFit3 is at the bounds' newline '!!!!!!!!!!!!!!!!!!!!!!']);
end

hFit = hFit0; hFit.Color = 'g'; hFit.LineWidth = 2;
hFit0 = plot(t, doubleGamma(t,p0), 'r--', 'LineWidth', 1);
hResid = plot(t, y - doubleGamma(t,pFit3), 'k:', 'LineWidth', 1);
uistack(hFit, 'top');


pX = pFit3;
pX(4:6) = 0;
plot(t, doubleGamma(t,pX), 'm');
pX = pFit3;
pX(1:3) = 0;
plot(t, doubleGamma(t,pX), 'm');


%%% Stage 4: Fit both gamma1 and gamma2
disp('Stage 4: Fitting both gamma1 and gamma2...');
% Parameter fixing mask for both gamma1 and gamma2 [AMP1 PEAK1 FWHM1 AMP2 PEAK2 FWHM2]
paramMask4 = logical([1 1 1 1 1 1]);  % Fit all parameters

% Initial guess
p0 = pFit3;

% Set bounds for parameters [AMP1 PEAK1 FWHM1 AMP2 PEAK2 FWHM2]
ub_full = [0.5*p0(1) p0(2)+5 p0(3)+2 5.0*p0(4) p0(5)+7 p0(6)+10];  % Upper bounds
lb_full = [2.0*p0(1) p0(2)-5 p0(3)-5 0.5*p0(4) p0(5)-7 p0(6)-10];  % Lower bounds
[ub_full; p0; lb_full]

figure('Position', [1200 500 400 400]);
hData = plot(t, y, 'k', 'LineWidth', 2); hold on
hFit0 = plot(t, doubleGamma(t,p0), 'm--', 'LineWidth', 1);
title('Fit both gamma1 and gamma2');
grid on;
axis tight

% Perform the fourth fit
pFit4 = fitDoubleGamma(t,y,paramMask4,p0,lb_full,ub_full,hFit0);
[ub_full; pFit4; lb_full]
if min([abs([pFit4(paramMask4)-ub_full(paramMask4) pFit4(paramMask4)-lb_full(paramMask4)])])<1e-6
    warning([newline '!!!!!!!!!!!!!!!!!!!!!!' newline 'pFit4 is at the bounds' newline '!!!!!!!!!!!!!!!!!!!!!!']);
end


hFit = hFit0; hFit.Color = 'g'; hFit.LineWidth = 2;
hFit0 = plot(t, doubleGamma(t,p0), 'r--', 'LineWidth', 1);
hResid = plot(t, y - doubleGamma(t,pFit4), 'k:', 'LineWidth', 1);
uistack(hFit, 'top');

pX = pFit4;
pX(4:6) = 0;
plot(t, doubleGamma(t,pX), 'm');
pX = pFit4;
pX(1:3) = 0;
plot(t, doubleGamma(t,pX), 'm');



end






%% Find fmristats parameters
hrf1 = doubleGamma(t,pFit4);
s1 = abs(sum(hrf1));
% hrf1 = hrf1./s1;
pFit4_fmristats = [pFit4([2 3 5 6]) abs(pFit4(4)/pFit4(1))];
X_cache =fmridesign(t,0,[1 0],[],pFit4_fmristats);
tr = mean(diff(t));
hrf2 = tr.*pFit4(1).*X_cache.X(:,1,1,1);
s2 = abs(sum(hrf2));
% hrf2 = hrf2./s2;

% figure('Position', [1200 500 400 400]);
% plot(t,hrf1,'k'); hold on
% plot(t,hrf2,'r--');

pFit4_fmristats(end+1) = s1/s2*tr.*pFit4(1);
% figure('Position', [1200 500 400 400]);
% plot(t,doubleGamma(t,pFit4),'k','LineWidth',2); hold on
% plot(t,pFit4_fmristats(end)*X_cache.X(:,1,1,1),'r--','LineWidth',1)

%% Remove scaling
pFit1([1 4]) = pFit1([1 4])./sFactor;
pFit2([1 4]) = pFit2([1 4])./sFactor;
pFit3([1 4]) = pFit3([1 4])./sFactor;
pFit4([1 4]) = pFit4([1 4])./sFactor;
pFit4_fmristats(end)   = pFit4_fmristats(end)./sFactor;

end 