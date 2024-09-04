function cmd = fsCommand2(runCond,stimCondLabel)
global srcFs

if ~exist('stimCondLabel','var'); stimCondLabel = []; end


runCondStimList = fields(runCond);
if isempty(stimCondLabel)
    runFlag = 0;
else
    runFlag = 1;
    runCondStimList = runCondStimList(ismember(runCondStimList,stimCondLabel));
end




cmd = {srcFs};
if runFlag
    cmd{end+1} = ['freeview -subtitle sub-' runCond.(runCondStimList{1}).sub '_task-' replace(char(runCondStimList),'_','-') ' \'];
else
    cmd{end+1} = ['freeview -subtitle sub-' runCond.(runCondStimList{1}).sub ' \'];
end



%% sub underlay
volAnatSub = fields(runCond); volAnatSub = runCond.(volAnatSub{1}).volAnatSub;
fList = cellstr(volAnatSub.mask{1}.ulay.f);
label = strjoin({'subUlay'},'_');
opt = {...
    ['name=' label]...
    };
for f = 1:length(fList)
    cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
end


%% resp individual run
if runFlag
    for rc = 1:length(runCondStimList)
        stimCondLabel = runCondStimList{rc};
        if ~isfield(runCond,runCondStimList{rc}); continue; end
        volResp = cat(1,runCond.(stimCondLabel).volTs.volResp);
        
        ses = [volResp.ts]; ses = {ses.fspec}'; [~,ses] = fileparts(fileparts(ses)); ses = cellstr(ses); for i = 1:size(ses,1); ses{i} = strsplit(ses{i},'_'); ses{i} = ses{i}{contains(ses{i},'ses-')}; end
        run = [volResp.ts]; run = {run.fspec}'; [~,run] = fileparts(fileparts(run)); run = cellstr(run); for i = 1:size(run,1); run{i} = strsplit(run{i},'_'); run{i} = run{i}{contains(run{i},'run-')}; end
    
        fList = [volResp.ts];
        fList = {fList.fspec};
        fBaseList = [volResp.base];
        fBaseList = {fBaseList.fspec};
        for f = 1:length(fList)
            label = strjoin({...
                replace(runCondStimList{rc},'_','-')...
                ses{f}...
                run{f}...
                'resp'...
                },'_');
            opt = {...
                ['name=' label]...
                'visible=0'};
            cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];

            label = strjoin({...
                replace(runCondStimList{rc},'_','-')...
                ses{f}...
                run{f}...
                'base'...
                },'_');
            opt = {...
                ['name=' label]...
                'visible=0'};
            cmd{end+1} = [strjoin([fBaseList(f) opt],':') ' \'];
            clear volResp
        end
    end
end


%% resp individual ses (cross-run average)
for rc = 1:length(runCondStimList)
    stimCondLabel = runCondStimList{rc};
    if ~isfield(runCond,runCondStimList{rc}); continue; end
    volResp = cat(1,runCond.(stimCondLabel).volTsSes.volResp);
    
    ses = [volResp.ts]; ses = {ses.fspec}'; [~,ses] = fileparts(fileparts(ses)); ses = cellstr(ses); for i = 1:size(ses,1); ses{i} = strsplit(ses{i},'_'); ses{i} = ses{i}{contains(ses{i},'ses-')}; end
    run = [volResp.ts]; run = {run.fspec}'; [~,run] = fileparts(fileparts(run)); run = cellstr(run); for i = 1:size(run,1); run{i} = strsplit(run{i},'_'); run{i} = run{i}{contains(run{i},'run-')}; end

    fList = [volResp.ts];
    fList = {fList.fspec};
    fBaseList = [volResp.base];
    fBaseList = {fBaseList.fspec};
    for f = 1:length(fList)
        label = strjoin({...
            replace(runCondStimList{rc},'_','-')...
            ses{f}...
            run{f}...
            'resp'...
            },'_');
        opt = {...
            ['name=' label]...
            'visible=0'};
        cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];

        label = strjoin({...
            replace(runCondStimList{rc},'_','-')...
            ses{f}...
            run{f}...
            'base'...
            },'_');
        opt = {...
            ['name=' label]...
            'visible=0'};
        cmd{end+1} = [strjoin([fBaseList(f) opt],':') ' \'];
    end
    clear volResp
end

%% resp individual sub (cross-ses cross-run average)
if ~runFlag
    for rc = 1:length(runCondStimList)
        stimCondLabel = runCondStimList{rc};
        if ~isfield(runCond,runCondStimList{rc}); continue; end
        volResp = cat(1,runCond.(stimCondLabel).volTsSub.volResp);

        ses = [volResp.ts]; ses = {ses.fspec}'; [~,ses] = fileparts(fileparts(ses)); ses = cellstr(ses); for i = 1:size(ses,1); ses{i} = strsplit(ses{i},'_'); ses{i} = ses{i}{contains(ses{i},'ses-')}; end
        run = [volResp.ts]; run = {run.fspec}'; [~,run] = fileparts(fileparts(run)); run = cellstr(run); for i = 1:size(run,1); run{i} = strsplit(run{i},'_'); run{i} = run{i}{contains(run{i},'run-')}; end

        fList = [volResp.ts];
        fList = {fList.fspec};
        for f = 1:length(fList)
            label = strjoin({...
                replace(runCondStimList{rc},'_','-')...
                ses{f}...
                run{f}...
                'resp'...
                },'_');
            opt = {...
                ['name=' label]...
                'visible=0'};
            cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
        end
        clear volResp
    end
end


%% respF individual run
if runFlag
    for rc = 1:length(runCondStimList)
        stimCondLabel = runCondStimList{rc};
        if ~isfield(runCond,runCondStimList{rc}); continue; end
        volResp = cat(1,runCond.(stimCondLabel).volTs.volResp);

        ses = [volResp.ts]; ses = {ses.fspec}'; [~,ses] = fileparts(fileparts(ses)); ses = cellstr(ses); for i = 1:size(ses,1); ses{i} = strsplit(ses{i},'_'); ses{i} = ses{i}{contains(ses{i},'ses-')}; end
        run = [volResp.ts]; run = {run.fspec}'; [~,run] = fileparts(fileparts(run)); run = cellstr(run); for i = 1:size(run,1); run{i} = strsplit(run{i},'_'); run{i} = run{i}{contains(run{i},'run-')}; end
    
        fList = [volResp.F];
        fList = {fList.fspec};
        for f = 1:length(fList)
            label = strjoin({...
                replace(runCondStimList{rc},'_','-')...
                ses{f}...
                run{f}...
                'respF'...
                },'_');
            opt = {...
                ['name=' label]...
                'visible=1'...
                'colormap=heat'...
                'heatscale=2.5,5'};
            cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
        end
        clear volResp
    end
end


%% respF individual ses (cross-run average)
for rc = 1:length(runCondStimList)
    stimCondLabel = runCondStimList{rc};
    if ~isfield(runCond,runCondStimList{rc}); continue; end
    volResp = cat(1,runCond.(stimCondLabel).volTsSes.volResp);
    
    ses = [volResp.ts]; ses = {ses.fspec}'; [~,ses] = fileparts(fileparts(ses)); ses = cellstr(ses); for i = 1:size(ses,1); ses{i} = strsplit(ses{i},'_'); ses{i} = ses{i}{contains(ses{i},'ses-')}; end
    run = [volResp.ts]; run = {run.fspec}'; [~,run] = fileparts(fileparts(run)); run = cellstr(run); for i = 1:size(run,1); run{i} = strsplit(run{i},'_'); run{i} = run{i}{contains(run{i},'run-')}; end
    
    fList = [volResp.F];
    fList = {fList.fspec};
    for f = 1:length(fList)
        label = strjoin({...
            replace(runCondStimList{rc},'_','-')...
            ses{f}...
            run{f}...
            'respF'...
            },'_');
        opt = {...
            ['name=' label]...
            'visible=1'...
            'colormap=heat'...
            'heatscale=2.5,5'};
        cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
    end
    clear volResp
end


%% respF individual sub (cross-ses cross-run average)
if ~runFlag
    for rc = 1:length(runCondStimList)
        stimCondLabel = runCondStimList{rc};
        if ~isfield(runCond,runCondStimList{rc}); continue; end
        volResp = cat(1,runCond.(stimCondLabel).volTsSub.volResp);
        
        ses = [volResp.ts]; ses = {ses.fspec}'; [~,ses] = fileparts(fileparts(ses)); ses = cellstr(ses); for i = 1:size(ses,1); ses{i} = strsplit(ses{i},'_'); ses{i} = ses{i}{contains(ses{i},'ses-')}; end
        run = [volResp.ts]; run = {run.fspec}'; [~,run] = fileparts(fileparts(run)); run = cellstr(run); for i = 1:size(run,1); run{i} = strsplit(run{i},'_'); run{i} = run{i}{contains(run{i},'run-')}; end
    
        fList = [volResp.F];
        fList = {fList.fspec};
        for f = 1:length(fList)
            label = strjoin({...
                replace(runCondStimList{rc},'_','-')...
                ses{f}...
                run{f}...
                'respF'...
                },'_');
            opt = {...
                ['name=' label]...
                'visible=1'...
                'colormap=heat'...
                'heatscale=2.5,5'};
            cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
        end
        clear volResp
    end
end





warning('PSD and COH is broken')

% %% PSD individual ses (cross-run average)
% for rc = 1:length(runCondStimList)
%     stimCondLabel = runCondStimList{rc};
%     if ~isfield(runCond,runCondStimList{rc}); continue; end
%     volPsd = cat(1,runCond.(stimCondLabel).volTsSes.volPsd);
% 
%     if isMRI(volPsd)
%         fList = volPsd;
% 
%         ses = volPsd;
%         run = volPsd;
%     else
%         fList = [volPsd.psd];
% 
%         ses = [volPsd.psd]; 
%         run = [volPsd.psd];
%     end
%     fList = {fList.fspec};
% 
%     ses = {ses.fspec}'; [~,ses] = fileparts(fileparts(ses)); ses = cellstr(ses); for i = 1:size(ses,1); ses{i} = strsplit(ses{i},'_'); ses{i} = ses{i}{contains(ses{i},'ses-')}; end
%     run = {run.fspec}'; [~,run] = fileparts(fileparts(run)); run = cellstr(run); for i = 1:size(run,1); run{i} = strsplit(run{i},'_'); run{i} = run{i}{contains(run{i},'run-')}; end
% 
%     % fBaseList = [volResp.base];
%     % fBaseList = {fBaseList.fspec};
%     for f = 1:length(fList)
%         label = strjoin({...
%             replace(runCondStimList{rc},'_','-')...
%             ses{f}...
%             run{f}...
%             'psd'...
%             },'_');
%         opt = {...
%             ['name=' label]...
%             'visible=0'...
%             'colormap=turbo'};
%         cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
% 
%         % label = strjoin({...
%         %     replace(runCondStimList{rc},'_','-')...
%         %     ses{f}...
%         %     run{f}...
%         %     'base'...
%         %     },'_');
%         % opt = {...
%         %     ['name=' label]...
%         %     'visible=0'};
%         % cmd{end+1} = [strjoin([fBaseList(f) opt],':') ' \'];
%     end
%     clear volPsd
% end
% 
% 
% 
% %% COH individual ses (cross-run average)
% for rc = 1:length(runCondStimList)
%     stimCondLabel = runCondStimList{rc};
%     if ~isfield(runCond,runCondStimList{rc}); continue; end
%     volPsd = cat(1,runCond.(stimCondLabel).volTsSes.volPsd);
% 
%     if isMRI(volPsd)
%         fList = volPsd;
% 
%         ses = volPsd;
%         run = volPsd;
%     else
%         fList = [volPsd.svd];
% 
%         ses = [volPsd.svd]; 
%         run = [volPsd.svd];
%     end
%     fList = {fList.fspec};
% 
%     ses = {ses.fspec}'; [~,ses] = fileparts(fileparts(ses)); ses = cellstr(ses); for i = 1:size(ses,1); ses{i} = strsplit(ses{i},'_'); ses{i} = ses{i}{contains(ses{i},'ses-')}; end
%     run = {run.fspec}'; [~,run] = fileparts(fileparts(run)); run = cellstr(run); for i = 1:size(run,1); run{i} = strsplit(run{i},'_'); run{i} = run{i}{contains(run{i},'run-')}; end
% 
% 
%     % ses = [volPsd.svd]; ses = [ses.fspec]'; ses = {ses.spSVmag}'; [~,ses] = fileparts(fileparts(ses)); ses = cellstr(ses); for i = 1:size(ses,1); ses{i} = strsplit(ses{i},'_'); ses{i} = ses{i}{contains(ses{i},'ses-')}; end
%     % run = [volPsd.svd]; run = [run.fspec]'; run = {run.spSVmag}'; [~,run] = fileparts(fileparts(run)); run = cellstr(run); for i = 1:size(run,1); run{i} = strsplit(run{i},'_'); run{i} = run{i}{contains(run{i},'run-')}; end
%     % 
%     % fList = [volPsd.svd];
%     % fList = [fList.fspec];
%     fListMag   = {fList.spSVmag}';
%     fListPhase = {fList.spSVphase}';
%     % fBaseList = [volResp.base];
%     % fBaseList = {fBaseList.fspec};
%     for f = 1:length(fList)
%         label = strjoin({...
%             replace(runCondStimList{rc},'_','-')...
%             ses{f}...
%             run{f}...
%             'cohMag'...
%             },'_');
%         opt = {...
%             ['name=' label]...
%             'visible=0'...
%             'colormap=turbo'};
%         cmd{end+1} = [strjoin([fListMag(f) opt],':') ' \'];
% 
%         label = strjoin({...
%             replace(runCondStimList{rc},'_','-')...
%             ses{f}...
%             run{f}...
%             'cohPhase'},'_');
%         opt = {...
%             ['name=' label]...
%             'visible=0'...
%             'colormap=turbo'};
%         cmd{end+1} = [strjoin([fListPhase(f) opt],':') ' \'];
% 
%         % label = strjoin({...
%         %     replace(runCondStimList{rc},'_','-')...
%         %     ses{f}...
%         %     run{f}...
%         %     'base'...
%         %     },'_');
%         % opt = {...
%         %     ['name=' label]...
%         %     'visible=0'};
%         % cmd{end+1} = [strjoin([fBaseList(f) opt],':') ' \'];
%     end
%     clear volPsd
% end





cmd{end}   = replace(cmd{end},' \',' &');
cmd{end+1} = newline;
clipboard('copy',strjoin(cmd,newline))
disp(strjoin(cmd,[' \\' newline]))







