function cmdFs = viewResp(volResp,volAnat,volFile,volFileOpt)
global srcAfni srcFs
if ~exist('volFile','var');       volFile = {}; end
if ~iscell(volFile);              volFile = {volFile}; end
if ~exist('volFileOpt','var'); volFileOpt = {}; end
if ~iscell(volFileOpt);        volFileOpt = {volFileOpt}; end


cmdFs = {srcFs};
cmdFs{end+1} = 'freeview \';
cmdFs{end+1} = [volResp.base.fspec ':visible=1 \'];
cmdFs{end+1} = [volResp.tsOnBase.fspec ':visible=0 \'];
cmdFs{end+1} = [volResp.ts.fspec ':visible=0 \'];
cmdFs{end+1} = [volResp.F.fspec ':visible=1:colormap=turbo \'];
cmdFs{end+1} = [volResp.Fq.fspec ':visible=0'];

%% Add roi
if ~isempty(volAnat)
    label = 'vesselCenters';
    volRoi = volAnat.fun.roi.vesselCenter;
    if isfield(volRoi,'mri') && isfield(volRoi.mri,'fspec')
        dbstack; error('code that')
    else
        mri = volAnat.fun.mask.crop;
        mri.fspec = [tempname '.nii.gz'];
        mri.vol = any(volRoi.mri.vol,4);
    end
    MRIwrite(mri,mri.fspec);
    cmdFs{end} = [cmdFs{end} ' \'];
    cmdFs{end+1} = [mri.fspec ':visible=0:name=' label];
end

%% Add extra
if ~isempty(volFile)
    cmdFs{end} = [cmdFs{end} ' \'];
    for i = 1:length(volFile)
        if isempty(volFileOpt)
            cmdFs{end+1} = [volFile{i} ':visible=0'];
        else
            cmdFs{end+1} = [volFile{i} ':visible=0' volFileOpt{i}];
        end
        if i ~= length(volFile)
            cmdFs{end} = [cmdFs{end} ' \'];
        end
    end
end



disp(strjoin(cmdFs,newline));
clipboard('copy',strjoin(cmdFs,newline))




% 
% 
% cmdFs = {srcAfni};
% cmdFs{end+1} = 'afni \';
% cmdFs{end+1} = [volResp.base.fspec ':visible=1 \'];
% cmdFs{end+1} = [volResp.tsOnBase.fspec ':visible=0 \'];
% cmdFs{end+1} = [volResp.ts.fspec ':visible=0 \'];
% cmdFs{end+1} = [volResp.F.fspec ':visible=1:colormap=turbo \'];
% cmdFs{end+1} = [volResp.Fq.fspec ':visible=0'];
% 
% %% Add roi
% if ~isempty(volAnat)
%     label = 'vesselCenters';
%     volRoi = volAnat.fun.roi.vesselCenter;
%     if isfield(volRoi,'mri') && isfield(volRoi.mri,'fspec')
%         dbstack; error('code that')
%     else
%         mri = volAnat.fun.mask.crop;
%         mri.fspec = [tempname '.nii.gz'];
%         mri.vol = any(volRoi.mri.vol,4);
%     end
%     MRIwrite(mri,mri.fspec);
%     cmdFs{end} = [cmdFs{end} ' \'];
%     cmdFs{end+1} = [mri.fspec ':visible=0:name=' label];
% end
% 
% disp(strjoin(cmdFs,newline));


