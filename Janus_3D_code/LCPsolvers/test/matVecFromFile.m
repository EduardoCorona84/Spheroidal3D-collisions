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

load(fullfile(dirname, '../../data/amphiLCPs.n_2.p_8.cDist_2.3.prt_0.mat'), 'Fparams', 'lcp_list');
%%
mc = 2;
F = lcp_list(mc).F;
nc = size(F,2);
C = lcp_list(mc).C;
p = 8;
A = getMatVec(Fparams, F, C, p);
disp(['Mat Vec for ix ' num2str(mc)])
% disp('run time ')

% rng('default')
% rng(1)
np = 1;
out = cell(np,1);
for n = 1:np
    disp(['n = ' num2str(n)])
    v = rand(nc, 1);
    tic
    u = A(v);
    toc
    % disp(' seconds')
    out{n} = u;
end
% save('OldOut.mat', 'out')
save('FastOut.mat', 'out')