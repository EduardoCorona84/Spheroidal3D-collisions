%% Get top level directory
[dirname, basedir] = setPaths();
%% Set paths to 
% mono-disperse 
root = fullfile(basedir,'offlineResults');
prefix = 'amphi.lcp.lattice.n_5.p_8.cDist_2.5';
% Hyper parameters 
ps = [8,6,4,3];
tols = [1e-6, 1e-6, 1e-6, 1e-5];
%% plots
finalComparison(root, prefix, ps, tols, false)
plotDenseMats(root, prefix)