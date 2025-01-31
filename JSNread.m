function out = JSNread(f,fld)

if ~exist('fld','var'); fld = []; end


if iscell(f)
    for r = 1:length(f)
        out{r} = JSNread(f{r},fld);
    end
elseif ischar(f)
    dbstack; error('code that')
    dir(replace(f,'.nii.gz','.jsn'))
    
else
    dbstack; error('X');
end

