function saveDenseMat(srcFile, dstDir, ix)
global DATA_DIR
%% Params (perhaps these should be inputs?)
ps = 8:-1:2; 
gmresTols = 10 .^(-(8:-1:4));
% ps = [2]; % 2:8;
% gmresTols = [10 .^(-4)];
%% set path
[dirname, ~] = setPaths();
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
%% Set DATA_DIR 
DATA_DIR = fullfile(getenv('SLURM_SCRATCH'), num2str(ix));
if ~exist(DATA_DIR, 'dir')
    mkdir(DATA_DIR)
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

disp(['nc = ' num2str(nc)])
try 
    load(dstFile, 'out');
    disp(['Loaded precomputed result from ' dstFile])
    ll = 1;
    kk = 1;
    for l = 1:numPs
        for k = 1:numTols
            if isempty(out{l,k})
                ll = l; kk=k;
                break;
            end
        end 
        if isempty(out{ll,kk})
            break;
        end
    end
    disp([' Found results up to ' num2str([ll,kk])])
catch 
    disp('No intermediate result found. Initializing with empty')
    out = cell(numPs,numTols);
    ll = 1;
    kk = 1;
end
for l = ll:numPs
    p = ps(l);
    disp(['l = ' num2str(l) ' <= ' num2str(numPs) ', p = ' num2str(p)])
    for k = kk:numTols
        gmresTol = gmresTols(k);
        disp(['  k = ' num2str(l) ' <= ' num2str(numTols) ', gmresTol = ' num2str(gmresTol)])
        Amatvec = getMatVec(Fparams, F, C, p, gmresTol);
        A = zeros(nc,nc);
        for i = 1:nc
            disp(['    i = ' num2str(i)])
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

