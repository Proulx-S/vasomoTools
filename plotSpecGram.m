function [ax,F] = plotSpecGram(volPsd,timeLabel,metricLabel,H)
if ~exist('H','var'); H = [];                             end
if isempty(H);        H = figure('WindowStyle','docked'); end

if isfield(volPsd,'dsgn') && isfield(volPsd.dsgn,'f0'); f0 = volPsd.dsgn.f0; else; f0 = []; end

if ~exist('timeLabel','var');     timeLabel = [];                             end
if ~exist('metricLabel','var'); metricLabel = [];                             end
if isempty(timeLabel);            timeLabel = 'gram'; end % 'gram' 'trialGram' 'trialGramMD'
if isempty(metricLabel);        metricLabel = 'psd'; end % 'psd' 'psdEPC' 'coh' 'cohEPC' 'cohEK'


switch class(H)
    case 'matlab.graphics.layout.TiledChartLayout'
        F = H.Parent;
    case 'matlab.ui.Figure'
        F = H;
    otherwise
end

figure(F);
ax = {};
ax{end+1} = nexttile;

%% Select approrpiate data
switch timeLabel
    case 'gram'
        switch metricLabel
            case 'psd'
                mt    = volPsd.psdGram;
                vec   = mt.PSD;
                label = 'spectrogram';
            case 'coh'
                mt    = volPsd.svdGram;
                vec   = mt.COH;
                label = 'coherogram';
            otherwise
                dbstack; error('code trhat')
        end
    case 'trialGram'
    %     case 'av'
    %     title(['evoked spectrogram ' paramStr])
    % case 'pc'
    %     title(['phase-locked evoked spectrogram ' paramStr])
    % case 'MD'
    %     title(['discontinuous-taper evoked spectrogram ' paramStr])
    % case 'MDpc'
    %     title(['discontinuous-taper phase-locked evoked spectrogram ' paramStr])
        switch metricLabel
            case 'psd'
                mt    = volPsd.psdTrialGram;
                vec   = mt.vec.psd;
                label = 'evoked spectrogram';
            case 'psdEPC'
                mt    = volPsd.psdTrialGram;
                vec   = mt.vec.psdPC;
                label = 'phase-locked evoked spectrogram';
            case 'coh'
                mt    = volPsd.svdTrialGram;
                vec   = mt.vec.coh;
                label = 'evoked coherogram';
            case 'cohEPC'
                mt    = volPsd.svdTrialGram;
                vec   = mt.vec.cohEPC;
                label = 'phase-locked evoked coherogram';
            case 'cohEK'
                mt    = volPsd.svdTrialGram;
                vec   = mt.vec.cohEK;
                label = 'evoked cross-trial coherogram';
            otherwise
                dbstack; error('code trhat')
        end
    case 'trialGramMD'
        switch metricLabel
            case 'psd'
                mt    = volPsd.psdTrialGramMD;
                vec   = mt.vec.psd;
                label = 'discontinuous-taper evoked spectrogram';
            case 'psdEPC'
                mt    = volPsd.psdTrialGramMD;
                vec   = mt.vec.psdPC;
                label = 'discontinuous-taper phase-locked evoked spectrogram';
            case 'coh'
                mt    = volPsd.svdTrialGramMD;
                vec   = mt.vec.coh;
                label = 'discontinuous-taper evoked coherogram';
            case 'cohEPC'
                mt    = volPsd.svdTrialGramMD;
                vec   = mt.vec.cohEPC;
                label = 'discontinuous-taper phase-locked evoked coherogram';                
            otherwise
                dbstack; error('code trhat')
        end
    otherwise
        dbstack; error('code trhat')
end



%% spatial average
vec = permute(mean(vec(:,:,:,:,:,:,:,1),6),[5 7 2 8 1 3 4 6]);
f   = permute(mt.f       ,[5 7 2 8 1 3 4 6]);
Fs  = mt.param.Fs;
t   = permute(mean(mt.t(:,1,:,:,:,:,:,1),1),[5 7 2 8 1 3 4 6]);
K   = mt.K;
T   = mt.T;
[TW,W] = K2W(T,K,0);

imagesc(t,f,vec)

xlabel('time (s)')
ylabel('freq (Hz)')
switch metricLabel
    case {'psd' 'psdEPC'}
        ylabel(colorbar,'psd');
        ax{end}.ColorScale = 'log';
    case {'coh' 'cohEPC' 'cohEK'}
        ylabel(colorbar,'coherence');
    otherwise
        dbstack; error('code trhat')
end




% if contains(timeLabel,'trial')
%     mt.t(:,1,:,:,:,:,:)
% else
runDur = volPsd.nframes/mt.param.Fs;
% end
xlim([0 runDur])

paramStr = ['(K=' num2str(K) '; 2W=' num2str(W*2,'%0.4f') 'Hz; T=' num2str(T,'%0.2f') 'sec; TW=' num2str(TW) ')'];
title([label ' ' paramStr])


cLim = vec(f>0.01,:);
cLim = [min(cLim(:)) max(cLim(:))];
clim(cLim)




addWin([],mt)
addW([],mt)




if isfield(volPsd,'psdTrialGram') && ~isempty(volPsd.psdTrialGram)
    addFreq([],volPsd.psdTrialGram.onsetList,volPsd.psdTrialGram.durList)
end

if ~isempty(f0)
    yline(f0,'--b')
end






ax = [ax{:}];


%% %%%%%%%%%%%%%%%%%%
% Add to timeseries %
%%%%%%%%%%%%%%%%%% %%

axTs = findobj(allchild(F.Children),'type','axes'); ttl = get(axTs,'Title'); ttl = get([ttl{:}],'String');
axTs = axTs(contains(ttl,'timeseries'));

%%% delete previous window size visual elements (magentat lines)
hLine = findobj(axTs.Children,'type','Line'); %hLine = {hLine(:)};
mLine = get(hLine,'Color'); if iscell(mLine); mLine = cat(1,mLine{:}); end
mLine = all(cat(1,mLine)==[1 0 1],2);
delete(hLine(mLine));

% % if length(hLine)==1
% %     hLine = {hLine};
% % end
% % 
% %     mLine = get(hLine,'Color');
% % else
% %     mLine = {get([axTs.Children(:)],'Color')};
% % end
% mLine = get(hLine,'Color'); if ~iscell(mLine); mLine = {mLine}; end
% ind = all(cat(1,mLine{:})==[1 0 1],2);
% mLine(ind)
% 
% mLine = hLine(all(cat(1,mLine{:})==[1 0 1],2));
% delete(mLine);
% % if length(axTs.Children)>1
% %     mLine = get(hLine,'Color');
% % else
% %     mLine = {get([axTs.Children(:)],'Color')};
% % end
% % mLine = axTs.Children(all(cat(1,mLine{:})==[1 0 1],2));
% % delete(mLine);

%%% add window size
addWin(axTs,mt)
