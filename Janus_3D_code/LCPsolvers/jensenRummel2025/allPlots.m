%% Set paths to 
% poly disperse 
% root = '/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord/';
% prefix = 'amphi.lcp.lattice.n_5.p_8.cDist_3.lcpSlvr_proxquasinewton.polyDisperseRatio_0.2';
% mono disperse (old)
root = '/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/goodData';
prefix = 'amphi.lcp.lattice.n_5.p_8.cDist_2.5';
%% 
ps = [8,6,4,3];
tols = [1e-6, 1e-6, 1e-6, 1e-5];
%% accuracy plots
monoFidelityComparison(root, prefix, ps, tols)
warmStartComparison(root, prefix, ps, tols)
bifidelityComparison(root, prefix, ps, tols, false) % no warm start
bifidelityComparison(root, prefix, ps, tols, true) % with warm start
%% Matrix Plots
% For now use these as a place holder
root = '/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord/';
prefix = 'amphi.lcp.lattice.n_5.p_8.cDist_3.lcpSlvr_proxquasinewton.polyDisperseRatio_0.2';
plotDenseMats(root, prefix)