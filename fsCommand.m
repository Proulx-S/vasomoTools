function cmd = fsCommand(runCond,acqCondLabel,stimCondLabel)
global srcFs

if ~exist('acqCondLabel','var');   acqCondLabel = []; end
if ~exist('stimCondLabel','var'); stimCondLabel = []; end

if isempty(acqCondLabel); acqCondLabel = 'vfMRI'; end


runCondAcqList = fields(runCond);
runCondAcq = runCondAcqList{ismember(runCondAcqList,acqCondLabel)};

runCondStimList = fields(runCond.(runCondAcq));
if isempty(stimCondLabel)
    runFlag = 0;
else
    runFlag = 1;
    runCondStimList = runCondStimList(ismember(runCondStimList,stimCondLabel));
end




cmd = {srcFs};
cmd{end+1} = ['freeview -subtitle ' runCond.(runCondAcq).(runCondStimList{1}).sub ' \'];


%%% sub underlay
% if ~runFlag
    volAnatSub = fields(runCond.(runCondAcq)); volAnatSub = runCond.(runCondAcq).(volAnatSub{1}).volAnatSub;
    fList = cellstr(volAnatSub.mask.ulay.f);
    label = strjoin({'subUlay'},'_');
    opt = {...
        ['name=' label]...
        };
    for f = 1:length(fList)
        cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
    end
% else
    % volAnat = fields(runCond.(runCondAcq)); volAnat = runCond.(runCondAcq).(volAnat{1}).volTs(1).volResp.base;
    % runCond.(runCondAcq).(volAnat{1}).volTsSes.mri
    % runCond.(runCondAcq).(volAnat{1}).fPreprocUnderSesCatRunCatAvList
    
% end



for rc = 1:length(runCondStimList)
    acqCondLabel = runCondStimList{rc};
    if ~isfield(runCond.(runCondAcq),runCondStimList{rc}); continue; end
    volTs = runCond.(runCondAcq).(acqCondLabel).volTs;
    volTsSes = runCond.(runCondAcq).(acqCondLabel).volTsSes;
    volResp = [volTs.volResp];
    volRespCat = [volTsSes.volRespCat];
    volRespSes = [volTsSes.volRespSes];
    ses = cellstr(strcat('ses-',runCond.(runCondAcq).(acqCondLabel).ses));
    run = runCond.(runCondAcq).(acqCondLabel).bidsList(:,contains(runCond.(runCondAcq).(acqCondLabel).bidsList(1,:),'run-'));

    %%% resp
    if runFlag
        %%%% individual run
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
                'respBase'...
                },'_');
            opt = {...
                ['name=' label]...
                'visible=0'};
            cmd{end+1} = [strjoin([fBaseList(f) opt],':') ' \'];
        end
    end

end



for rc = 1:length(runCondStimList)
    acqCondLabel = runCondStimList{rc};
    if ~isfield(runCond.(runCondAcq),runCondStimList{rc}); continue; end
    volTs = runCond.(runCondAcq).(acqCondLabel).volTs;
    volTsSes = runCond.(runCondAcq).(acqCondLabel).volTsSes;
    volResp = [volTs.volResp];
    volRespCat = [volTsSes.volRespCat];
    volRespSes = [volTsSes.volRespSes];
    ses = cellstr(strcat('ses-',runCond.(runCondAcq).(acqCondLabel).ses));
    run = runCond.(runCondAcq).(acqCondLabel).bidsList(:,contains(runCond.(runCondAcq).(acqCondLabel).bidsList(1,:),'run-'));

    %%% resp
    %%%% run catenated
    fList = [volRespCat.ts];
    fList = {fList.fspec};
    for f = 1:length(fList)
        label = strjoin({...
            replace(runCondStimList{rc},'_','-')...
            % ses{f}...
            'respCat'...
            },'_');
        opt = {...
            ['name=' label]...
            'visible=0'};
        cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
    end

end



for rc = 1:length(runCondStimList)
    acqCondLabel = runCondStimList{rc};
    if ~isfield(runCond.(runCondAcq),runCondStimList{rc}); continue; end
    volTs = runCond.(runCondAcq).(acqCondLabel).volTs;
    volTsSes = runCond.(runCondAcq).(acqCondLabel).volTsSes;
    volResp = [volTs.volResp];
    volRespCat = [volTsSes.volRespCat];
    volRespSes = [volTsSes.volRespSes];
    ses = cellstr(strcat('ses-',runCond.(runCondAcq).(acqCondLabel).ses));
    run = runCond.(runCondAcq).(acqCondLabel).bidsList(:,contains(runCond.(runCondAcq).(acqCondLabel).bidsList(1,:),'run-'));

    % if ~runFlag
    %     %%% resp
    %     %%%% run averaged
    %     fList = [volRespSes.ts];
    %     fList = {fList.fspec};
    %     for f = 1:length(fList)
    %         label = strjoin({...
    %             replace(runCondStimList{rc},'_','-')...
    %             % ses{f}...
    %             'respAv'...
    %             },'_');
    %         opt = {...
    %             ['name=' label]...
    %             'visible=0'};
    %         cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
    %     end
    % end

end



for rc = 1:length(runCondStimList)
    acqCondLabel = runCondStimList{rc};
    if ~isfield(runCond.(runCondAcq),runCondStimList{rc}); continue; end
    volTs = runCond.(runCondAcq).(acqCondLabel).volTs;
    volTsSes = runCond.(runCondAcq).(acqCondLabel).volTsSes;
    volResp = [volTs.volResp];
    volRespCat = [volTsSes.volRespCat];
    volRespSes = [volTsSes.volRespSes];
    ses = cellstr(strcat('ses-',runCond.(runCondAcq).(acqCondLabel).ses));
    run = runCond.(runCondAcq).(acqCondLabel).bidsList(:,contains(runCond.(runCondAcq).(acqCondLabel).bidsList(1,:),'run-'));

    if runFlag
        %%% respF
        %%%% individual run
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
    end

end



for rc = 1:length(runCondStimList)
    acqCondLabel = runCondStimList{rc};
    if ~isfield(runCond.(runCondAcq),runCondStimList{rc}); continue; end
    volTs = runCond.(runCondAcq).(acqCondLabel).volTs;
    volTsSes = runCond.(runCondAcq).(acqCondLabel).volTsSes;
    volResp = [volTs.volResp];
    volRespCat = [volTsSes.volRespCat];
    volRespSes = [volTsSes.volRespSes];
    ses = cellstr(strcat('ses-',runCond.(runCondAcq).(acqCondLabel).ses));
    run = runCond.(runCondAcq).(acqCondLabel).bidsList(:,contains(runCond.(runCondAcq).(acqCondLabel).bidsList(1,:),'run-'));

    %%% respF
    %%%% run catenated
    fList = [volRespCat.F];
    fList = {fList.fspec};
    for f = 1:length(fList)
        label = strjoin({...
            replace(runCondStimList{rc},'_','-')...
            % ses{f}...
            'respFcat'...
            },'_');
        opt = {...
            ['name=' label]...
            'visible=1'...
            'colormap=heat'...
            'heatscale=2.5,5'};
        cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
    end

end



for rc = 1:length(runCondStimList)
    acqCondLabel = runCondStimList{rc};
    if ~isfield(runCond.(runCondAcq),runCondStimList{rc}); continue; end
    volTs = runCond.(runCondAcq).(acqCondLabel).volTs;
    volTsSes = runCond.(runCondAcq).(acqCondLabel).volTsSes;
    volResp = [volTs.volResp];
    volRespCat = [volTsSes.volRespCat];
    volRespSes = [volTsSes.volRespSes];
    ses = cellstr(strcat('ses-',runCond.(runCondAcq).(acqCondLabel).ses));
    run = runCond.(runCondAcq).(acqCondLabel).bidsList(:,contains(runCond.(runCondAcq).(acqCondLabel).bidsList(1,:),'run-'));


    % if ~runFlag
    %     %%% respF
    %     %%%% run averaged
    %     fList = [volRespSes.F];
    %     fList = {fList.fspec};
    %     for f = 1:length(fList)
    %         label = strjoin({...
    %             replace(runCondStimList{rc},'_','-')...
    %             % ses{f}...
    %             'respFav'...
    %             },'_');
    %         opt = {...
    %             ['name=' label]...
    %             'visible=1'...
    %             'colormap=heat'...
    %             'heatscale=2.5,5'};
    %         cmd{end+1} = [strjoin([fList(f) opt],':') ' \'];
    %     end
    % end



end
cmd{end} = replace(cmd{end},' \','');
clipboard('copy',strjoin(cmd,newline))

disp(strjoin(cmd,[' \\' newline]))







