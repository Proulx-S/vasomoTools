function viewMTsvd(anaType,funTs,fpass,id,f0,mask,kInd)

if ~exist('f0','var') || isempty(f0)
    f0 = 0.1;
end
if ~exist('id','var') || isempty(id)
    id = '';
elseif isnumeric(id)
    id = num2str(id);
end
if ~exist('mask','var') || isempty(mask)
    mask = funTs.(anaType).mask;
end
if ~exist('slc','var') || isempty(slc)
    slc = round(funTs.depth*0.5);
end
f = funTs.(anaType).f;
coh = funTs.(anaType).c;
w = funTs.(anaType).w;

%% Prepare titles
if exist('id','var') && ischar(id) && ~isempty(id)
    titleStr     = {['subj' id]};
else
    titleStr = {''};
end
% titleStr{end+1} = ['tw=' num2str(funPsd.psd.tw)];
titleStr{end+1} = ['w=' num2str(w)];
[~,b,c] = fileparts(funTs.fspec);
titleStr{end+1} = [b c];
titleStr = [strjoin(titleStr(1:2),'; ') titleStr(3)];


switch anaType
    case 'svdMitra'
        figure('WindowStyle','docked');
        tl = tiledlayout(2,1);
        tl.Padding = 'tight';
        tl.TileSpacing = 'tight';
        tlTitle = title(tl,titleStr,'interpreter','none');

        %% PSD
        nexttile
        psd = funTs.psd.psd(mask(funTs.psd.mask),:);
        normFact = exp( mean(log(psd),2) - mean(log(psd(:))));
        psd = exp(mean(log(psd./normFact),1));
        hPsd = plot(f,psd,'k'); hold on
        ax = gca; ax.YScale = 'log';
        ax.XScale = 'log';
%         ax.YAxisLocation = 'right';
        xlabel('Hz'); ylabel({'psd averaged within mask' 'same as for svd'});
        grid on
        if exist('fpass','var') && ~isempty(fpass)
            xlim(fpass)
        else
            xlim tight
        end
        ylim tight

        f0w = exp(linspace(log(f(2)),log(f(end)),6)); f0w([1 end]) = [];
        yLim = ylim;
        y = yLim(1).*[1 1];
        for ii = 1:length(f0w)
            x = f0w(ii)+[-1 1].*funTs.psd.w;
            if x(1)>f(1) && x(2)<f(end)
                hW = plot(x,y,'g','LineWidth',3)
            end
        end

        fpassSVD = [-1 1].*w + funTs.(anaType).param.fpass(1);
        yLim = ylim;
        hfpass = patch(fpassSVD([1 2 2 1 1]),yLim([1 1 2 2 1]),'r','LineStyle','none','FaceAlpha',0.1);
        uistack(hfpass,'bottom')
        
        legend([hPsd hW hfpass],{'psd averaged over mask' '2*w' 'frequency range for svd'},'box','off')


        %% Scree plot
        nexttile
        plot(squeeze(funTs.svdMitra.sv),'k')
        axis tight
        hold on
        ylabel('singular values');
        xlabel('sorted mode index');
        grid on
        if ~exist('kInd','var') || isempty(kInd)
            return
        end
        plot([1 1].*kInd+0.5,ylim,'-r')
        text(kInd+1,funTs.svdMitra.sv(kInd),['keeping ' num2str(kInd) ' modes'])
        
        

        %% Launch 3D viewer
        funTs = vol2vec(funTs,mask);
        funTs.vec = abs(permute(funTs.(anaType).sp(:,:,1:kInd),[3 1 2]));
        funTs.nframes = kInd;
        funTs.tr = 1000;
        funTs = vec2vol(funTs);
        funTs.vol(isnan(funTs.vol)) = 0;
        funTmpName = [tempname '.nii.gz'];
        MRIwrite(funTs,funTmpName);
        display('******************')
        display(['To view 3D volume of the first ' num2str(kInd) ' spatial modes :'])
        display(['freeview ' funTmpName ':colormap=turbo'])
        display('******************')


    case 'svdKlein'
        funTs = vol2vec(funTs,mask);

        [~,f0Ind] = min(abs(f-f0));
        funTs.vec = funTs.svdKlein.sp(:,f0Ind,1)';
        funTs.nframes = 1;
        funTs = vec2vol(funTs);

        [a,b] = max(abs(funTs.vol(:)));
        disp(a)
        [x,y,slc] = ind2sub(size(funTs.vol),b);
        slc = [-1 0 1]+slc;
        spTmp = funTs.vol(:,:,slc);
        cLim = [min(min(min(abs(spTmp(~isnan(spTmp)))))) max(max(max(abs(spTmp(~isnan(spTmp))))))];
        clear spTmp

        figure('WindowStyle','docked');
        tl = tiledlayout(3,3); tl.TileSpacing = "tight"; tl.Padding = "tight";
        tl.TileIndexing = 'rowmajor';
        % Title
        if exist('id','var') && ischar(id) && ~isempty(id)
            titleStr     = {['subj' id]};
        else
            titleStr = {''};
        end
        titleStr{end+1} = ['tw=' num2str(funTs.psd.tw)];
        titleStr{end+1} = ['w=' num2str(w)];
        [~,b,c] = fileparts(funTs.fspec);
        titleStr{end+1} = [b c];
        titleStr = [strjoin(titleStr(1:3),'; ') titleStr(4)];
        tlTitle = title(tl,titleStr,'interpreter','none');


        for slcInd = 1:length(slc)
            maskC = getMaskOutline(mask(:,:,slc(slcInd)),5);
            titleStr2 = {['f0=' num2str(f(f0Ind),'%0.4f') 'Hz']};
            titleStr2{end+1} = ['slc' num2str(slc(slcInd))];

            nexttile
            imagesc(abs(funTs.vol(:,:,slc(slcInd))),cLim)
            ylabel(colorbar,['singular vector weight @ f0=' num2str(f(f0Ind),'%0.4f') 'Hz'])
            ax = gca;
            %             ax.ColorScale = 'log';
            ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1]; ax.XTick = []; ax.YTick = []; ax.YDir = 'normal';
            xlabel(strjoin(titleStr2,'; '))
            hold on; plot(maskC,'FaceColor','none')

            plot([1 1].*y,ylim,':w','LineWidth',0.1)
            plot(xlim,[1 1].*x,':w','LineWidth',0.1)
        end




        nexttile(4,[2 3])
        plot(f,coh,'-k')
        xlabel('Hz'); ylabel('spatial coherence');
        grid on
        xlim tight
        if exist('fpass','var') && ~isempty(fpass)
            xlim(fpass)
        end
        hold on
        K = size(funTs.svdKlein.sv,3);
        plot(xlim,[1 1].*1/K,'--r')

        % plot([1 1].*f(f0Ind),ylim,':r');

        f0w = exp(linspace(log(f(2)),log(f(end)),6)); f0w([1 end]) = [];
        yLim = ylim;
        y = yLim(1).*[1 1];
        for ii = 1:length(f0w)
            x = f0w(ii)+[-1 1].*w;
            if x(1)>f(1) && x(2)<f(end)
                plot(x,y,'-g','LineWidth',8)
            end
        end
        ylim tight

        % x = f(f0Ind)+[-1 1].*w;
        % yLim = ylim;
        % y = yLim(2).*[1 1];
        % plot(x,y,'g','LineWidth',3)

        if isfield(funTs,'psd') && isfield(funTs.psd,'psd')
            yyaxis right
            psd = exp(mean(log(funTs.psd.psd),1));
            plot(funTs.psd.f,psd,'-')
            ylim tight
            %     if exist('fpass','var') && ~isempty(fpass)
            %         [~,bMin] = min(abs(funTs.psd.f-fpass(1)));
            %         [~,bMax] = min(abs(funTs.psd.f-fpass(2)));
            %         ylim([min(psd(bMin:bMax)) max(psd(bMin:bMax))])
            %     else
            %         ylim([min(psd) max(psd)])
            %     end
        end

        ax = gca;
        ax.YScale = 'log';
        ax.XScale = 'log';
        drawnow
        set(ax.YAxis,'LimitsMode','manual')

        ylabel('psd averaged within mask');

        plot([1 1].*f(f0Ind),ylim,':r')
        x = f(f0Ind)+[-1 1].*w;
        yLim = ylim;
        y = yLim(2).*[1 1];
        plot(x,y,'-g','LineWidth',8)
        legend({'coherence' 'f0' '2*w'},'box','off','Location','northeast')

        
        
        figure('WindowStyle','docked');
        rho = corr(funTs.svdKlein.sp(:,:,1));
        imagesc(funTs.svdKlein.f(2:end),funTs.svdKlein.f(2:end),abs(rho(2:end,2:end)))
        ax = gca;
        ax.PlotBoxAspectRatio = [1 1 1];
        ax.DataAspectRatio = [1 1 1];
        ax.YDir = 'normal';
        ylabel(colorbar,'Pearson''s rho')
        xlabel('Hz')
        ylabel('Hz')
        title(titleStr,'Interpreter','none')


        drawnow
        %% Launch 3D viewer
        funTs = vol2vec(funTs);
        funTs.vec = abs(funTs.svdKlein.sp(:,f0Ind,1)');
        funTs.nframes = size(funTs.vec,1);
        funTs.tr = 1000;
        funTs = vec2vol(funTs);
        funTs.vol(isnan(funTs.vol)) = 0;
        funTmpName = [tempname '.nii.gz'];
        MRIwrite(funTs,funTmpName);
        display('******************')
        display(['To view 3D volume of the first spatial singular vector @ ' num2str(f(f0Ind)) 'Hz :'])
        display(['freeview ' funTmpName ':colormap=turbo'])
        display('******************')
end
