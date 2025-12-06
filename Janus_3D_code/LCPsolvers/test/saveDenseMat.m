function saveDenseMat(srcFile, dstDir, ix, p, gmresTol)
global DATA_DIR
%% set path
[dirname, ~] = setPaths();
%% set defaults
if ~exist('srcFile', 'var') || isempty(srcFile)
    srcFile = fullfile(dirname, '../data/amphi.lattice.n_5.p_8.cDist_2.5.mat');
end
if ~exist('dstDir', 'var') || isempty(dstDir)
    dstDir = fullfile(dirname, '../data/amphi.lattice.n_5.p_8.cDist_2.5');
end
if ~exist('ix', 'var') || isempty(ix)
    ix = 598;
end
if ~exist('p', 'var') || isempty(p)
    p = 8;
end
if ~exist('gmresTol', 'var') || isempty(gmresTol)
    gmresTol = 1e-8;
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
A = zeros(nc,nc);
dt = zeros(nc,1);
i = 1;
if exist(dstFile, 'file')
    res_ = load(dstFile);
    disp(['Loaded precomputed result from ' dstFile])
    if isempty(res_.A)
        i =1;
    else
        i = nc+1;
        for ii = 1:nc
            if all(res_.A(:,ii) == 0)
                disp(['---- intermediate results found up to column'  num2str(ii) '/' num2str(nc)]);
                i = ii;
                break
            end
            A(:,ii) = res_.A(:,ii);
            try
                dt(ii) = res_.dt(ii);
            catch
                i = 1;
                break;
            end
        end
    end
    % Grep... is so slow
    % if flag 
        % dt = getRunTimesFromLog(dstDir, ix, p, gmresTol);
    % end
    disp([' Found results up to ' num2str(i-1) '/' num2str(nc)]);
end

Amatvec = getMatVec(Fparams, F, C, p, gmresTol);
for ii = i:nc
    disp(['    ii = ' num2str(ii)])
    ei = zeros(nc,1);
    ei(ii) = 1;
    tic
    A(:,ii) = Amatvec(ei);
    dt(ii) = toc;
    if norm(A(:,ii)) > 1e4 
        disp('Numerical Error in GMRES, so perturbing input slightly.')
        ei = ei + rand(nc,1)*eps;
        A(:,ii) = Amatvec(ei);
    end
    disp(['dt = ' num2str(dt(ii))])
    disp(['Saving to ' dstFile])
    save(dstFile, 'A', 'p', 'gmresTol','dt')
end
disp('final save')
save(dstFile, 'A', 'p', 'gmresTol','dt')



