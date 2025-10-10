function saveDenseMat(srcFile, dstDir, ix)
%% Params (perhaps these should be inputs?)
ps = 8:-1:2; 
gmresTols = 10 .^(-(8:-1:4));
% ps = [2]; % 2:8;
% gmresTols = [10 .^(-4)];
%% set path
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname, ~,~] = fileparts(mfilePath);
basedir = fullfile(dirname, '..','..');
addpath(basedir);
addpath(fullfile(basedir,'support'));
addpath(genpath(fullfile(basedir, 'LCPsolvers/solvers')))
addpath(fullfile(basedir, 'FMMLIB/fmmlib3d-1.2/matlab'));
addpath(fullfile(basedir,'FMMLIB/stfmmlib3d-1.2/matlab'));
%% set defaults
if ~exist('srcFile', 'var') || isempty(srcFile)
    srcFile = fullfile(dirname, '../data/amphiLCPs.n_2.p_8.cDist_2.3.mat');
end
if ~exist('dstDir', 'var') || isempty(dstDir)
    dstDir = fullfile(dirname, '../data/amphiLCPs.n_2.p_8.cDist_2.3');
end
if ~exist('ix', 'var') || isempty(ix)
    ix = 1;
end
%% Load from file
disp(['Loading ' srcFile])
load(srcFile, 'Fparams', 'lcp_list');
%%
disp(['Creating dense matrices for ix = ' num2str(ix)])
dstFile = fullfile(dstDir, [num2str(ix) '.mat']);
if ~exist(dstDir, 'dir')
    mkdir(dstDir)
end
F = lcp_list(ix).F;
nc = size(F,2);
C = lcp_list(ix).C;
numPs = numel(ps);
numTols = numel(gmresTols);
out = cell(numPs,numTols);
disp(['nc = ' num2str(nc)])
for l = 1:numPs
    p = ps(l);
    disp(['p = ' num2str(p)])
    for k = 1:numTols
        gmresTol = gmresTols(k);
        Amatvec = getMatVec(Fparams, F, C, p, gmresTol);
        A = zeros(nc,nc);
        for i = 1:nc
            disp(['  i = ' num2str(i)])
            ei = zeros(nc,1);
            ei(i) = 1;
            tic
            A(:,i) = Amatvec(ei);
            dt = toc;
            disp(['dt = ' num2str(dt)])
        end
        out{l,k} = A;
        disp(['Saving to ' dstFile])
        save(dstFile, 'out', 'ps', 'gmresTols')
    end
end

