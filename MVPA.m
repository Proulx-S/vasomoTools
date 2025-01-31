function [volTs,volResp]= MVPA(volTs,info,mask,volResp)


%% Load time series data if not already loaded
for r = 1:length(volTs)
    if isempty(volTs(r).mri.vol) && ~isfield(volTs(r).mri,'vec') || isempty(volTs(r).mri.vec)
        disp(['loadind run' num2str(r) '/' num2str(length(volTs))])
        volTs(r).mri = MRIload2(volTs(r).mri,mask);
    end
end

%% Run univariate response extraction
if ~exist('volResp','var') || isempty(volResp)
    forceThis   = 1;
    verboseThis = 1;
    if isfield(volTs,'dsgn') && ~isempty(volTs(1).dsgn.onsetList)
        [volResp, ~, info] = volTsGetResp3([],info,volTs,[],mask,forceThis,verboseThis);
    else
        dbstack; error('need dsgn as a subfield of volTs')
    end
end

%% Load response time series  if not already loaded
for r = 1:length(volResp)
    if isempty(volResp(r).ts.vol) && ~isfield(volResp(r).ts,'vec') || isempty(volResp(r).ts.vec)
        disp(['loadind run' num2str(r) '/' num2str(length(volTs))])
        volResp(r).ts = MRIload2(volResp(r).ts,mask);
    end
end


%% 
t = volResp.ts.t(2:end-1)';
imMask = MRIload2(MRIload2(mask)); imMask = imMask.vol2vec;
imxLim = [find(any(imMask,1),1,'first')-1.5 find(any(imMask,1),1,'last')+1.5];
imyLim = [find(any(imMask,2),1,'first')-1.5 find(any(imMask,2),1,'last')+1.5];

switch info.method
    case 'uniSVD'
        Y = volResp.ts.vec(2:end-1,:)';
        
        [U,S,V] = svd(Y,'vector','econ');
        
        Cvar = S';
        Ctime = permute(V.*Cvar,[3 1 2]);
        Cspace = permute(U.*Cvar,[1 3 2]);
        Cvar = permute(Cvar,[1 3 2]);

    case 'multiSVD'
        Y = [volTs.mri]; Y = cat(1,Y.vec);
        X = volResp.dsgn.dsgnMat;
        
        [beta,Sigma,E,CovB,logL] = mvregress(X,Y);
        
        regInd = ismember(volResp.dsgn.dsgnMatLabel(1,:),'stim');
        Y = beta(regInd,:)';
        [U,S,V] = svd(Y,'vector','econ');

        Cvar = S';
        Ctime = permute(V.*Cvar,[3 1 2]);
        Cspace = permute(U.*Cvar,[1 3 2]);
        Cvar = permute(Cvar,[1 3 2]);

        % Actually gives the same response time courses as uniSVD. The
        % difference may lie only in the error covariance and associated
        % statistics.
        % See
        % https://www.mathworks.com/matlabcentral/answers/108929-after-using-mvregress-how-can-i-find-the-rsquared-value-t-values-p-values-f-statistic-and-stand
        % for statistic on multivariate models

    case 'canon'
        Y = [volTs.mri]; Y = cat(1,Y.vec);
        X = volResp.dsgn.dsgnMat;
        [beta,Sigma,E,CovB,logL] = mvregress(X,Y);
        regInd = ~ismember(volResp.dsgn.dsgnMatLabel(1,:),'stim');
        Y = Y-X(:,regInd)*beta(regInd,:);
        
        regInd = ismember(volResp.dsgn.dsgnMatLabel(1,:),'stim');
        [A,B,r,U,V,stats] = canoncorr(X(:,regInd),Y);
        

        A = ( U' / (   X(:,regInd)-mean(X(:,regInd),1)   )'  )';
        B = ( V' / (   Y          -mean(Y          ,1)   )'  )';

        Cvar = permute(stats.F,[ 1 3 2]);
        Ctime = permute(A,[3 1 2]);
        Cspace = permute(B,[1 3 2]);

        regInd = ismember(volResp.dsgn.dsgnMatLabel(1,:),'stim');
        Y = beta(regInd,:)';

        % Since regressor estimates are the same for univariate and
        % multivariate regression, using multivariate regression for
        % detrending is a waste.

    case 'pls'
        Y = [volTs.mri]; Y = cat(1,Y.vec);
        X = volResp.dsgn.dsgnMat;
        [beta,Sigma,E,CovB,logL] = mvregress(X,Y);
        regInd = ~ismember(volResp.dsgn.dsgnMatLabel(1,:),'stim');
        Y = Y-X(:,regInd)*beta(regInd,:);

        regInd = ismember(volResp.dsgn.dsgnMatLabel(1,:),'stim');
        [XL,YL,XS,YS,~,PCTVAR,MSE,stats] = plsregress(X(:,regInd),Y);

        Cvar = PCTVAR(2,:);
        Ctime = permute(XL,[3 1 2]);
        Cspace = permute(YL,[1 3 2]);

        regInd = ismember(volResp.dsgn.dsgnMatLabel(1,:),'stim');
        Y = beta(regInd,:)';

        % Since regressor estimates are the same for univariate and
        % multivariate regression, using multivariate regression for
        % detrending is a waste.
end


figure('WindowStyle','docked');
ht = tiledlayout(3,3); ht.TileSpacing = 'tight'; ht.Padding = "tight";

%%% Plot original data
nexttile;
plot(t,Y); hold on
plot(t,mean(Y,1),'k','LineWidth',3);
xlabel('time (sec)')
title('voxel response timecourses')

nexttile;
plot(squeeze(Cvar),'-ok')
ylabel('variance explained')
xlabel('component number')

nexttile;
imagesc(imMask); colormap gray
xlim(imxLim); ylim(imyLim);
ax = gca; ax.PlotBoxAspectRatio = [1 1 1];
ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off';
title('voxel mask')

%%% Plot spactial compoentns
axSpace = {};
for c = 1:3
    nexttile;
    im = zeros(size(imMask));
    im(imMask) = Cspace(:,:,c);
    imagesc(im)
    xlim(imxLim); ylim(imyLim);
    ax = gca; ax.PlotBoxAspectRatio = [1 1 1];
    ax.YAxis.Visible = 'off'; ax.XAxis.Visible = 'off';
    ax.Colormap = parula; colorbar
    title(['component ' num2str(c)])
    axSpace{end+1} = ax;
end
cLim = get([axSpace{:}],'CLim');
cLim = max(abs([cLim{:}])).*[-1 1];
set([axSpace{:}],'CLim',cLim);

%%% Plot temporal compoentns
axTime = {};
for c = 1:3
    axTime{end+1} = nexttile;
    % switch info.method
    %     case {'uniSVD' 'multiSVD'}
            plot(t,Ctime(:,:,c),'k'); hold on
    %     otherwise
    %         plot(t,Ctime(:,:,c),'k'); hold on
    % end

    if c==1
        [~,iVox] = max(abs(Cspace(:,:,c)),[],1);
        y = Y(iVox,:).*sign(Cspace(iVox,:,c));
        y = y ./max(abs(y)) .*max(abs(Ctime(:,:,c)),[],2);
        plot(t,y,'--k')
    end

    grid on; grid minor;
    axis tight
    title(['component ' num2str(c)])
    xlabel('time (sec)')
    drawnow
    xLim = xlim; xLim(1) = 0; xlim(xLim);

    if c==1
        legend({'compenent' 'max vox'})
    end
end
hold on
plot(t,sum(Ctime(:,:,1:3),3),'r')
yLim = get([axTime{:}],'YLim'); yLim = max(abs([yLim{:}])).*[-1 1];
set([axTime{:}],'YLim',yLim);

ttlStr = {info.method};
if isfield(info,'label') && ~isempty(info.label)
    ttlStr{end+1} = info.label;
end
title(ht,strjoin(ttlStr,'; '));


drawnow

