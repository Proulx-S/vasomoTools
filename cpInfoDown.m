function [down,up] = cpInfoDown(up,down,fieldList)

if ~exist('fieldList','var')
    fieldList = []                                                                              ; end
if isempty(fieldList)
    fieldList = {'sub' 'ses' 'acq' 'task' 'dsgn' 'acqTime' 'nDummy' 'nFrame' 'nFrameOrig' 'fOrigList'}; end


for i = 1:length(fieldList)


    %% Reconcile old and new convention
    switch fieldList{i}
        case 'acq'
            if ~isfield(up,fieldList{i})  ||  isempty(up.(fieldList{i})) || strcmp(up.(fieldList{i}),'?')
                if isfield(up,'labelAcq')  &&  ~isempty(up.labelAcq)
                    up.(fieldList{i}) = up.labelAcq;
                end
            end
        case 'task'
            if ~isfield(up,fieldList{i})  ||  isempty(up.(fieldList{i})) || strcmp(up.(fieldList{i}),'?')
                if isfield(up,'label')  &&  ~isempty(up.label)
                    up.(fieldList{i}) = up.label;
                end
            end
    end

    
    %% Skip otherwise non-existent fields
    % if ~isfield(up,fieldList{i})  ||  isempty(up.(fieldList{i}))  ; continue; end
    if isempty(up.(fieldList{i}))                                 ; continue; end
    if isfield(down,fieldList{i}) && ~isempty(down.(fieldList{i})); continue; end
    
    
    
    %% Copy fields

    % One up, many down
    if size(up.(fieldList{i}),1)==1 && size(down,1)~=1
        [down.(fieldList{i})] = deal(up.(fieldList{i}));
        continue
    end

    % Many up, as many down
    if all(  size(up.(fieldList{i})) == size(down)  )
        for ii = 1:numel(down)
            down(ii).(fieldList{i}) = up.(fieldList{i})(ii);
        end
        continue
    end

end
