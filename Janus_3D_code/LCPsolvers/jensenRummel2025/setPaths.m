function [dirname, basedir] = setPaths()
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname, ~,~] = fileparts(mfilePath);
basedir = getAbsPath(fullfile(dirname, '..','..'));
addpath(basedir);
addpath(dirname);
addpath(fullfile(basedir,'support'));
addpath(genpath(fullfile(basedir, 'LCPsolvers/solvers')))
addpath(genpath(fullfile(basedir, 'LCPsolvers/test')))
addpath(fullfile(basedir, 'FMMLIB/fmmlib3d-1.2/matlab'));
addpath(fullfile(basedir,'FMMLIB/stfmmlib3d-1.2/matlab'));

function absPath = getAbsPath(obj)
getAbsFolderPath = @(y) string(unique(arrayfun(@(x) x.folder, dir(y), 'UniformOutput', false)));
getAbsFilePath = @(y) string(arrayfun(@(x) fullfile(x.folder, x.name), dir(y), 'UniformOutput', false));
if isfolder(obj)
    absPath = getAbsFolderPath(obj);
elseif isfile(obj)
    absPath = getAbsFilePath(obj);
else
    error('The specified object does not exist.');
end