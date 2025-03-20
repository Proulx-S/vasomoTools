function out = anyDontExist(fList)
    out = false;
    for i = 1:length(fList)
        if ~isempty(fList{i}) && ~exist(fList{i},'file')
            out = true;
            break;
        end
    end
end