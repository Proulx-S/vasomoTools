function viewMTsg(funTs,anaType,coord,f0,id)
wS = 0.005;
wSG = 0.01;
winSz = 200;
SGscalePSD = 'log';
SscalePSD = 'log';
scaleF = 'log';
fLow = 0.01;

% if ~exist('f0','var') || isempty(f0)
%     f0 = 0.1;
% end
if ~exist('f0','var')% || isempty(f0)
    f0 = [];
end
if ~exist('id','var') || isempty(id)
    id = '';
elseif isnumeric(id)
    id = num2str(id);
end
% if ~exist('mask','var') || isempty(mask)
%     mask = funTs.(anaType).mask;
% end
if ~exist('slc','var') || isempty(slc)
    slc = round(funTs.depth*0.5);
end
% f = funTs.(anaType).f;
% coh = funTs.(anaType).c;
% w = funTs.(anaType).w;
tr = funTs.tr/1000;
T = funTs.nframes*tr;
t = 0:tr:(T-tr);


%% Prepare titles
if exist('id','var') && ischar(id) && ~isempty(id)
    titleStr     = {['subj' id]};
else
    titleStr = {''};
end
% titleStr{end+1} = ['tw=' num2str(funPsd.psd.tw)];
% titleStr{end+1} = ['w=' num2str(w)];
[~,b,c] = fileparts(funTs.fspec);
titleStr{end+1} = [b c];
titleStr = [strjoin(titleStr(1:end-1),'; ') titleStr(end)];

switch anaType
    case 'singleVox'
        %% Make sure coord are matching
        if ~isempty(coord)
            figure('WindowStyle','docked');
            if isfield(funTs,'funBrain')
                funMean = funTs.funBrain.vol;
            else
                funMean = mean(funTs.vol,4);
            end
            imagesc(funMean(:,:,coord(3)));
            ax = gca;
            ax.Colormap = gray;
            ax.YDir = 'normal';
            ax.PlotBoxAspectRatio = [1 1 1]; ax.DataAspectRatio = [1 1 1];
            ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
            hold on
            plot([1 1].*coord(1),ylim,':r')
            plot(xlim,[1 1].*coord(2),':r')
            title(titleStr,'interpreter','none');
        end

        %% Initiate figure
        figure('WindowStyle','docked');
        tl = tiledlayout(4,4);
        tl.Padding = 'none';
        tl.TileSpacing = 'none';
        tlTitle = title(tl,titleStr,'interpreter','none');
        %% Define timeseries and plot
        mask = false(funTs.volsize);
        mask(coord(2),coord(1),coord(3)) = true; % flip x and y to match matlab imagesc
%         mask(coord(1),coord(2),coord(3)) = true; % flip x and y to match matlab imagesc
        funTs = vol2vec(funTs,mask);
        
        ax = {};
        nexttile(1,[3,1])
%         subplot(4,4,[1 5 9])
        plot(funTs.vec,t,'k'); ylim([0 t(end)])
        ylabel('time (s)')
        xlabel('0-mean BOLD')
        ax{end+1} = gca;
        ax{end}.YDir = 'reverse';
        box off
        grid on
        yLim = ylim;
        


        %% Compute spectrum and plot
        disp('Spectrum')
        [paramS,w] = setTaperParam(wS,[],funTs);
        paramS.err = [1 0.05];
        [S, F, Serr] = mtspectrumc(funTs.vec, paramS);
        nexttile(14,[1,3])
        hEr = shadedErrorBar(F(2:end),S(2:end)',[Serr(1,2:end)-S(2:end)'; S(2:end)'-Serr(2,2:end)]);
        delete(hEr.edge)
        if isempty(fLow)
            xlim([F(2) F(end)])
        else
            xlim([fLow F(end)])
        end
        hold on
        ax{end+1} = gca;
        ax{end}.YAxis.Scale = SscalePSD;
        ax{end}.XAxis.Scale = scaleF;
        ax{end}.YAxisLocation = 'right';
        grid on
        box off
        xLim = xlim;
        xlabel('Hz')
        ylabel('PSD')
        %%% Add w references
        addRef(w)
        

        %% Plot spectrogram
        disp('Spectrogram')
        [paramSG,w] = setTaperParam(wSG,[],funTs,winSz);
        [Sg,tg,fg] = mtspecgramc(funTs.vec, [winSz winSz/10], paramSG);
        tg = tg + t(1);
        nexttile(2,[3,3])
%         imagesc(fg,tg,Sg)
        surf(fg(2:end),tg,Sg(:,2:end),'EdgeColor','none','FaceColor','interp'); view(0,90)
        ax{end+1} = gca;
        ax{end}.YDir = 'reverse';
        ax{end}.ColorScale = SGscalePSD;
        ax{end}.XScale = scaleF;
        ax{end}.Color = 'none';
        ax{end}.YAxis.Visible = 'off';
        ax{end}.XAxis.Visible = 'off';
        grid on
        grid minor
        box off
        xlim(xLim)
        ylim(yLim)
        ylabel(colorbar,'psd')
        %%% Plot bw
        addRef(w,[],winSz)
        
    otherwise
        error('Allowed anaType: ''singleVox''')
end
