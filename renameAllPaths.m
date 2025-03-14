function var = renameAllPaths(var, oldName, newName)
varInfo = whos('var');

switch varInfo.class
    case {'char'}
        if size(var,1) == 1
            var = replace(var,oldName,newName);
        else
            for i = 1:size(var,1)
                var(i,:) = renameAllPaths(var(i,:), oldName, newName);
            end
        end
    case {'cell'}
        for i = 1:numel(var)
                var{i} = renameAllPaths(var{i}, oldName, newName);
        end
    case {'struct', 'runCond', 'runDsgn'}
        if numel(var) == 1
            fields = fieldnames(var);
            for i = 1:length(fields)
                    var.(fields{i}) = renameAllPaths(var.(fields{i}), oldName, newName);
            end
        else
            for i = 1:numel(var)
                var(i) = renameAllPaths(var(i), oldName, newName);
            end
        end
    case {'double', 'single', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64',...
        'datetime', 'logical','duration',...
        'matlab.ui.Figure'}
    otherwise
        % keyboard
        dbstack; error('code that');
end