function results = compareLofiPrec()
%% Hyper parameters
prefix = 'amphi.lattice.n_5.p_8.cDist_2.5.multiFidelity';
plotDebug = true;
minSz = 100;
maxSz = inf;
max_iter = 1000;
tol = 1e-8;
%% Set paths
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
%% Try to initialize cvx if on the cluster
try %#ok<TRYNC>
    run(fullfile('/projects', getenv('USER'), 'cvx', 'cvx_startup.m'))
end
%% Preallocate structures
opts = struct( ...
    'max_iter',max_iter, ...
    'kkt_rel',tol, ...
    'kkt_abs',tol, ...
    'storeIts', true, ...
    'subspaceMin', struct( ...
    'innerStepSelection', 'pqn', ...
    'innerSolver','cvx', ...
    'orthoMethod','qr', ...
    'kThresh', 0, ...
    'useOneStepIter', false ...
    ) ...
    );
%% Switch load LCP's saved in a file
res_ = load('/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/data.11.13.2025/amphi.lattice.n_5.p_8.cDist_2.5.allMats.mat');
[~,jHi] = max(res_.ps);
[~,kHi]= min(res_.tols);
Nt = size(res_.A,1);
ttlStr = [sprintf('A($p=%d', res_.ps(jHi)) ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g)', res_.tols(kHi)) '$'];
A_list = res_.A(:,jHi,kHi); % arrayfun(@(i) (res_.A{i,jHi,kHi} + res_.A{i,jHi,kHi}')/2, 1:Nt, 'UniformOutput',false);
b_list = res_.b;
dt_list = arrayfun(@(i) res_.dt{i,jHi,kHi} , 1:Nt, 'UniformOutput',false);
[pLo,jLo] = min(res_.ps);
[tolLo,kLo]= max(res_.tols);
IX = find(cellfun(@(A) ~isempty(A), res_.A(:,jHi,kHi)) & ...
    cellfun(@(A) ~isempty(A), res_.A(:,jLo,kLo))); 
% cellfun(@(A) ~isempty(A), res_.A(:,5,k)) &...
% cellfun(@(A) ~isempty(A), res_.A(:,4,end)) & ...
% cellfun(@(A) ~isempty(A), res_.A(:,3,end)) & ...
% cellfun(@(A) ~isempty(A), res_.A(:,2,end)) & ...

I = numel(IX);
%% Set up all algorithm callers
algoNames = {
    'PGD ($\kappa = \tau_{bb_1}, \eta = 1$)';
    'PQN (BFGS)';
    };
algoSlvr = {
    'pgd';
    'proxquasinewton';
    };
algoHndls = {
    @projectedGradientDescent;
    @proxQuasiNewton;
    };
% TODO CHECK MORE
mRatio = .5;
for p = pLo
    for tol = tolLo %, 1e-6, 1e-7, 1e-8]
        jMid = find(res_.ps == p,1,'first');
        kMid = find(abs(res_.tols - tol) < eps*10, 1,'first');
        % if any(arrayfun(@(i) any(eig((res_.A{i,jMid,kMid} +res_.A{i,jMid,kMid}) /2) <0) ,IX))
        %     continue 
        % end
        % Use the full Ahat (basically cheat)
        algoNames{end+1} = ['PQN (BFGS, $B_0 = \hat{A}(' sprintf('p=%d', p) ...
            ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g', tol) '$)']; %#ok<*AGROW>
        algoSlvr{end+1} = 'proxquasinewton';
        algoHndls{end+1} = @(fg, x0, opts) precProxQuasiNewton(fg, x0, opts,...
            res_.A{opts.i,jMid,kMid}, res_.b{opts.i}, ...
            'fullAhat',p,tol,mRatio);
        % Use the truncated svd to approximate Ahat 
        % mPercent = round(mRatio * 100);
        % algoNames{end+1} = ['PQN (BFGS, $B_0 = \hat{A}_{' sprintf('%d',mPercent) ...
        %     '\%}(' sprintf('p=%d', p) ...
        %     ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g', tol) '$)']; %#ok<*AGROW>
        % algoSlvr{end+1} = 'proxquasinewton';
        % algoHndls{end+1} = @(fg, x0, opts) precProxQuasiNewton(fg, x0, opts,...
        %      res_.A{opts.i,jMid,kMid}, res_.b{opts.i}, ...
        %     'truncAhat',p,tol,mRatio);
        % Use the BFGS to get a low rank approx of Ahat
        algoNames{end+1} = ['PQN (BFGS, $B_0 = B(' sprintf('p=%d', p) ...
            ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g', tol) '$)'];
        algoSlvr{end+1} = 'proxquasinewton';
        algoHndls{end+1} = @(fg, x0, opts) precProxQuasiNewton(fg, x0, opts,...
             res_.A{opts.i,jMid,kMid}, res_.b{opts.i}, ...
            'solveLoFi',p,tol,mRatio);
    end
end
%%
numAlgo = numel(algoNames);
assert(numel(algoNames) == numel(algoHndls));
results = repmat(...
    struct( ...
    'algo', '',...
    'time',   [], ...
    'estimTime',   [], ...
    'iters', [], ...
    'kkt', [], ...
    'matVecs', [], ...
    'AhatMVP', [], ...
    'x', [], ...
    'errHist', [], ...
    'iterHist', [] ...
    ), [I,numAlgo] ...
    );
mcGood = false(I,1);
for ii = 1:I
    i = IX(ii);
    disp(['i = ' num2str(ii) '/' num2str(I)])
    %% Load problem from list
    A = A_list{i};
    dt = dt_list{i};
    b = b_list{i};
    n = size(A,2);
    if isempty(A) || n < minSz || maxSz <= n
        continue;
    end
    %% Build cost function
    Acnt = @(x) Acounter(x,A, false);
    x0 = zeros(n,1);
    fg = @(x, Ax) quadraticLoss(x, Acnt, b, Ax);
    %% Fill the opts with problem specific information
    opts.i = i;
    opts.A = Acnt;
    opts.b = b;
    xlb = A \ -b;
    opts.fstar = 1/2*dot(xlb, A*xlb) + dot(xlb, b);
    evals = eig(A);
    opts.L = max(evals);
    opts.mu = min(evals);
    opts.AA = A;
    try %#ok<TRYNC>
        [xstar,~] = callCVX([], x0, opts);
        opts.metricNames = {'abs_kkt', 'rel_kkt', 'MVP', 'rel_iter', 'abs_iter', 'obj'};
        opts.errFcn = {
            @(x) abs_kkt(x, A, b);
            @(x) rel_kkt(x, A, b);
            @(x) Acnt('cnt');
            @(x) rel_iter(x,xstar);
            @(x) abs_iter(x,xstar);
            @(x) dot(x,0.5*A*x+b);
            };
        mcGood(i) = true;
    end
    %% Run all the algorithms
    for ixAlgo = 1:numAlgo
        if ~mcGood(i)
            break
        end
        name = algoNames{ixAlgo};
        this_opts= opts;
        this_opts.name = name;
        this_opts.solver = algoSlvr{ixAlgo};
        [this_opts,~] = defaultLCPOpts(this_opts, x0);
        results(i,ixAlgo).estimTime = 0;
        Acnt('reset');
        algo = algoHndls{ixAlgo};
        tic
        [x, info] = algo(fg, x0, this_opts);
        results(i,ixAlgo).name = name;
        results(i,ixAlgo).x = x;
        results(i,ixAlgo).time = toc();
        results(i,ixAlgo).matVecs = Acnt('reset');
        results(i,ixAlgo).estimTime = results(i,ixAlgo).estimTime + results(i,ixAlgo).matVecs*mean(dt);
        results(i,ixAlgo).iters = info.iter;
        results(i,ixAlgo).kkt = info.kkt;
        results(i,ixAlgo).errHist = info.errHist;
        try %#ok<TRYNC>
            results(i,ixAlgo).iterHist{i} = info.iterHist;
        end
        for ixMetric = 1:numel(opts.errFcn)
            hndl = opts.errFcn{ixMetric};
            try %#ok<TRYNC>
                hndl('reset');
            end
        end
    end
    %% Debug Plot
    plotSubProblem(results(i,:), A, mcGood(i) && plotDebug, [ttlStr sprintf(', i = %d', i)]);
end
%% Print debug information
fprintf('Algo            | Time      | matVec    | Iter | kkt\n')
for ixAlgo = 1:numAlgo
    name = algoNames{ixAlgo};
    while length(name) < length('Projected QuasiNewton (BFGS)')
        name = [name ' ']; %#ok<AGROW>
    end
    time = [results(mcGood,ixAlgo).time];
    iters = [results(mcGood,ixAlgo).iters];
    matVec = [results(mcGood,ixAlgo).matVecs];
    kkt = [results(mcGood,ixAlgo).kkt];
    fprintf('%s \t| %.1e s | %.3g\t| %.3g | %.2g\n', ...
        name(1:10), mean(time), mean(matVec), mean(iters), mean(kkt));
end
tmp = mcGood;
mcGood = find(tmp);
%% Plot Overall Statistics
iterationBarChart
%% Save results to File
disp(['Saving to ' prefix]);
save(['/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/' prefix '.mat'], ...
    'results', 'mcGood', 'minSz', 'maxSz', 'opts');
end


function [x, info] = precProxQuasiNewton(fg, x0, opts, AHat, b, mode, p,tol, mRatio)
AHatCnt = @(x) Acounter(x, AHat, false);
ACnt = opts.A;
opts.A = AHatCnt;
fgMid =  @(x, Ax) quadraticLoss(x, AHatCnt, b, Ax);
[x0, ~, opts] = proxQuasiNewton(fgMid, x0, opts);
opts.A = ACnt;
for ixMetric = 1:numel(opts.errFcn)
    hndl = opts.errFcn{ixMetric};
    try %#ok<TRYNC>
        hndl('reset');
    end
end
switch mode
    case 'fullAhat'
        n = opts.n;
        S = eye(n);
        Y = AHat;
        opts.qn.rho = arrayfun(@(i) 1/dot(S(:,i), Y(:,i)), 1:n);
        opts.qn.m = 2*n;
        opts.qn.S = [S zeros(n)];
        opts.qn.Y = [Y zeros(n)];
    case 'truncAhat'
        n = size(AHat,1); 
        m = round(mRatio*n);
        [S, L] = eig(AHat);
        Y = S*L;
        Y = Y(:,end-m+1:end);
        S = S(:,end-m+1:end);
        lambdas = diag(L);
        opts.qn.rho = 1.0./lambdas(end-m+1:end);
        opts.qn.m = m+n;
        opts.qn.S = [S zeros(n)];
        opts.qn.Y = [Y zeros(n)];
    case 'solveLoFi'
       n = size(AHat,1);
       r = find(opts.qn.rho, 1,'last');
       opts.qn.m = r +n;
       opts.qn.rho = [opts.qn.rho(1:r); zeros(n,1)];
       opts.qn.S = [opts.qn.S(:,1:r); zeros(n)];
       opts.qn.Y = [opts.qn.Y(:,1:r); zeros(n)];
    otherwise
        error(['Not implemented mode ' mode])
end
opts.A('reset');
[x, info] = proxQuasiNewton(fg, x0, opts);
end