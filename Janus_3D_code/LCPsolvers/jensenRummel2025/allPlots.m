%% Get top level directory
[dirname, basedir] = setPaths();
%% Set paths to 
% poly-disperse 
root = '/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord.01.31.2026/';
prefix = 'amphi.lcp.lattice.n_5.p_8.cDist_3.lcpSlvr_proxquasinewton.polyDisperseRatio_0.2';
% mono-disperse 
% root = fullfile(basedir,'goodData');
% prefix = 'amphi.lcp.lattice.n_5.p_8.cDist_2.5';
% Hyper parameters 
ps = [8,6,4,3];
tols = [1e-6, 1e-6, 1e-6, 1e-5];
%% plots
finalComparison(root, prefix, ps, tols, false)
% plotDenseMats(root, prefix)