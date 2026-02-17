function saveDenseMat(srcFile, dstDir, ix, p, gmresTol)
global DATA_DIR
%% set path
[dirname, ~] = setPaths();
%% set defaults
if ~exist('srcFile', 'var') || isempty(srcFile)
    srcFile = '/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/test/../resultsForRecord/amphi.lcp.lattice.n_5.p_8.cDist_3.lcpSlvr_proxquasinewton.polyDisperseRatio_0.2.mat';
end
if ~exist('dstDir', 'var') || isempty(dstDir)
    dstDir = '/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/test/../resultsForRecord/amphi.lcp.lattice.n_5.p_8.cDist_3.lcpSlvr_proxquasinewton.polyDisperseRatio_0.2.denseMats';
end
if ~exist('ix', 'var') || isempty(ix)
    ix = 50;
end
if ~exist('p', 'var') || isempty(p)
    p = 4;
end
if ~exist('gmresTol', 'var') || isempty(gmresTol)
    gmresTol = 1e-6;
end
%% Set DATA_DIR 
DATA_DIR = fullfile(getenv('SLURM_SCRATCH'), ...
    ['ix_' num2str(ix) '.p_' num2str(p) '.tol_' num2str(gmresTol)]);
if ~exist(DATA_DIR, 'dir')
    mkdir(DATA_DIR)
end
%% Load from file
disp(['Loading ' srcFile])
load(srcFile, 'Fparams', 'lcp_list');
%%
disp(['Creating dense matrices for ix = ' num2str(ix) ...
    ', p = ' num2str(p) ', tol = ' num2str(gmresTol)])
dstFile = fullfile(dstDir, ['ix_' num2str(ix) '.p_' num2str(p) '.tol_' num2str(gmresTol) '.mat']);
if ~exist(dstDir, 'dir')
    mkdir(dstDir)
end
F = lcp_list(ix).F;
nc = size(F,2);
C = lcp_list(ix).C;
disp(['nc = ' num2str(nc)])

if nc ==0
    disp('Empty LCP problem')
    A = [];
    save(dstFile, 'A', 'p', 'gmresTol')
    return
end
IX = 1:nc;
A = zeros(nc,nc);
dt = zeros(nc,1);
if exist(dstFile, 'file')
    res_ = load(dstFile);
    disp(['Loaded precomputed result from ' dstFile])
    if ~isempty(res_.A)
        IX = [];
        for i = 1:nc
            if all(res_.A(:,i) == 0) || norm(res_.A(:,i)) > 4
                IX(end+1) = i; %#ok<AGROW> 
            else
                A(:,i) = res_.A(:,i);
                dt(i) = res_.dt(i);
            end
        end
    end
    % Grep... is so slow
    % if flag 
        % dt = getRunTimesFromLog(dstDir, ix, p, gmresTol);
    % end
    disp([' Need to compute ' num2str(numel(IX)) '/' num2str(nc)]);
    fprintf('  IX=%d\n' ,IX')
end

Amatvec = getLCPMatVec(Fparams, F, [], [], p, gmresTol, C, false, false);
for ii = 1:numel(IX)
    i = IX(ii);
    disp(['    i = ' num2str(i)])
    ei = zeros(nc,1);
    ei(i) = 1;
    tic
    A(:,i) = Amatvec(ei);
    dt(i) = toc;
    if norm(A(:,ii)) > 4 
        disp('Numerical Error in GMRES, so perturbing input slightly.')
        ei = ei + rand(nc,1)*eps;
        A(:,i) = Amatvec(ei);
    end
    disp(['dt = ' num2str(dt(i))])
    disp(['Saving to ' dstFile])
    save(dstFile, 'A', 'p', 'gmresTol','dt')
end
disp('final save')
save(dstFile, 'A', 'p', 'gmresTol','dt')



