[dirname, basedir] = setPaths();
load(fullfile(basedir, 'data.10.31.2025/amphi.lattice.n_5.p_8.cDist_2.5.allMats.mat'))
%%
MC = size(A, 1);
numP = size(A, 2);
numTols = size(A,3);
jHi = numP -1;
kHi = numTols;
jLo = numP -2;
kLo = numTols;
ii = find(arrayfun(@(i) ~isempty(A{i,jHi, kHi}) && ~isempty(A{i,jLo, kLo}), 1:MC), 1, 'last');
%%
bb = b{ii};
pHi = ps(jHi);
tolHi = tols(kHi);
AHi = A{ii,jHi, kHi};
AHi = (AHi + AHi') /2;
dtHi = dt{ii, jHi, kHi};
pLo = ps(jLo);
tolLo = tols(kLo);
ALo = A{ii,jLo, kLo};
ALo = (ALo + ALo') /2;
sqrtALoinv = inv(sqrtm(ALo));
dtLo = dt{ii, jLo, kLo};
%
fprintf('RelErr = %.3g\n', norm(AHi - ALo) / norm(AHi))
fprintf('Cond of A_{hi} = %.3g\n', cond(AHi))
fprintf('Cond of A_{lo}^{-1} A_{hi} = %.3g\n', cond(ALo \ AHi))
fprintf('Cond of A_{lo}^{-1/2} A_{hi} A_{lo}^{-1/2} = %.3g\n', cond(sqrtALoinv*AHi*sqrtALoinv))
%%
n = size(AHi,1);
max_iter = 1000;
tol = 1e-8;
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
algoNames = {
    'Min-Map Newton';
    'PGD (\kappa = \tau_{bb_1}, \eta = 1)'; 'PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    'zeroSR1 (\kappa = 1, \eta = \eta^*)';
    'L-BFGS-B';
    'PQN (BFGS, \kappa = 1, \eta = \eta^*)';
    'Subspace Minimization (PQN)';
    };
algoSlvr = {
    'semismoothnewton';
    'pgd'; 'pgd';
    'zerosr1';
    'l-bfgs-b';
    'proxquasinewton';
    'subspacemin';
    };
algoHndls = {
    @minmap_newton;
    @projectedGradientDescent; @projectedGradientDescent;
    @zeroSr1_nic;
    @L_BFGS_B;
    @proxQuasiNewton;
    @subspaceMin;
    };
numAlgo = numel(algoNames);
assert(numel(algoNames) == numel(algoHndls));
results = repmat(...
    struct( ...
    'algo', '',...
    'time',   [], ...
    'iters', [], ...
    'kkt', [], ...
    'matVecs', [], ...
    'x', [], ...
    'errHist', [], ...
    'iterHist', [] ...
    ), [1,numAlgo] ...
    );

Acnt = @(x) Acounter(x,AHi, false);
x0 = zeros(n,1);
fg = @(x, Ax) quadraticLoss(x, Acnt, bb, Ax);
%% Fill the opts with problem specific information
opts.A = Acnt;
opts.b = bb;
xlb = AHi \ -bb;
opts.fstar = 1/2*dot(xlb, AHi*xlb) + dot(xlb, bb);
evals = eig(AHi);
opts.L = max(evals);
opts.mu = min(evals);
opts.AA = AHi;
mcGoodFlag = true;
try %#ok<TRYNC>
    [xstar,~] = callCVX([], x0, opts);
    opts.metricNames = {'abs_kkt', 'rel_kkt', 'MVP', 'rel_iter', 'abs_iter', 'obj'};
    opts.errFcn = {
        @(x) abs_kkt(x, AHi, bb);
        @(x) rel_kkt(x, AHi, bb);
        @(x) Acnt('cnt');
        @(x) rel_iter(x,xstar);
        @(x) abs_iter(x,xstar);
        @(x) dot(x,0.5*AHi*x+bb);
        };
    mcGood(mc) = true;
end
%% Run all the algorithms
Acnt('reset')
for ixAlgo = 1:numAlgo
    this_opts = opts;
    name = algoNames{ixAlgo};
    this_opts.name = name; 
    this_opts.solver = algoSlvr{ixAlgo};
    [this_opts,~] = defaultLCPOpts(this_opts, x0);
    algo = algoHndls{ixAlgo};
    tic
    [x, info] = algo(fg, x0, this_opts);
    results(ixAlgo).name = name;
    results(ixAlgo).x = x;
    results(ixAlgo).time = toc();
    results(ixAlgo).iters = info.iter;
    results(ixAlgo).kkt = info.kkt;
    results(ixAlgo).matVecs = Acnt('reset');
    results(ixAlgo).errHist = info.errHist;
    try %#ok<TRYNC>
        results(ixAlgo).iterHist = info.iterHist;
    end
    % Reset all the metric handles so that we don't mess up any persistent
    % var silly-ness
    for ixMetric = 1:numel(opts.errFcn)
        hndl = opts.errFcn{ixMetric};
        try %#ok<TRYNC>
            hndl('reset');
        end
    end
end
plotSubProblem(results, AHi)