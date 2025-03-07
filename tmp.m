cd /Users/sebastienproulx/Library/CloudStorage/OneDrive-Stanford/vsm
addpath('vasomoTools')

saveFlag = 1;

% Get list of files to process
files = dir('rCond_*.mat');
nFiles = length(files);

% Loop through each file
iFile = 1%:nFiles
%% Load data
% Get current filename
currentFile = files(iFile).name;
fprintf('Processing file %d/%d: %s\n', iFile, nFiles, currentFile);

% Load data
rCond = load(currentFile);


oldName = '/autofs/space/takoyaki_001/users';
newName = '/cluster/meso/users/seb';
rCond = renameAllPaths(rCond, oldName, newName)

rCond.rCond_vsmDrivenP1.vfMRI.task_05sPrd1sDur.fPreprocList


