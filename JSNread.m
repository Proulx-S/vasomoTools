function out = JSNread(f,fld)
if ~exist('fld','var'); fld = []; end


if iscell(f)

    if isempty(fld)
        JSNread(f{1});
        return;
    end
    for r = 1:length(f)
        out{r} = JSNread(f{r},fld);
    end
    out = cat(2,out{:})';


elseif ischar(f)

    fid = fopen(replace(f,'.nii.gz','.json'), 'r');
    if fid == -1; error('Could not open JSON file: %s', j); end
    j = jsondecode(fread(fid, '*char')'); fclose(fid);
    
    if isempty(fld)
        disp('--------------------');
        disp('Fields in json file:');
        disp(char(fields(j)));
        disp('--------------------');
        return;
    end

    out = cell(length(fld),1);
    for i = 1:length(fld)
        if isfield(j,fld{i})
            out{i} = j.(fld{i});
        end
    end
else
    dbstack; error('X');
end

