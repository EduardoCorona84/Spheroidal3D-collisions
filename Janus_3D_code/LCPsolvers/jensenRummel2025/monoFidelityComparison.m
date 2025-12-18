function monoFidelityComparison(root, prefix, ps, tols)
%% Set paths
[dirname, basedir] = setPaths();
%% Try to initialize cvx if on the cluster
try %#ok<TRYNC>
    run(fullfile('/projects', getenv('USER'), 'cvx', 'cvx_startup.m'))
end
%% Load Data from file
resFile = fullfile(root, [prefix '.denseMats.allMats.mat']);
res = load(resFile);
pHi = ps(1);
tolHi = tols(1);
jHi = find(abs(res.ps - pHi) < 1e-8);
kHi = find(abs(res.tols - tolHi) < 1e-8);
Nt = size(res.A,1);
%% Hyper parameters
plotDebug = false;
condNum = 1e2;
percentLarge = 0.25;
max_iter = 1000;
tol = 1e-8;
% initialize opts structure
opts = struct( ...
    'max_iter',max_iter, ...
    'kkt_rel',tol, ...
    'kkt_abs',tol, ...
    'storeIts', true);
%% Specify all algorithms to be compared
algoNames = {
    'PGD';
    'Accelerated PGD';
    'zeroSR1';
    'L-BFGS-B';
    'PQN';
    'Min-Map Newton';
    };
algoSlvr = {
    'bbpgd';
    'fista'
    'zerosr1';
    'l-bfgs-b';
    'proxquasinewton';
    'semismoothnewton';
    };
algoHndls = {
    @projectedGradientDescent; 
    @fista;
    @zeroSr1_nic;
    @L_BFGS_B;
    @proxQuasiNewton;
    @minmap_newton;
    };
numAlgo = numel(algoNames);
assert(numel(algoNames) == numel(algoHndls));
%% Initialize results struct
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
    ), [Nt,numAlgo] ...
    );
%% Run all solvers on all problems 
mcGood = getMCGood(resFile, ps, tols);
badII = [];
for ii = 1:numel(mcGood)
    i = mcGood(ii);
    disp(['- i = ' num2str(i) '/' num2str(Nt)])
    A = res.A{i,jHi,kHi};
    dt = res.dt{i,jHi,kHi};
    b = res.b{i};
    n = size(A,2);
    if isempty(A) || norm(A) > 10 
        % the first check is if this problem doesn't exist (LCP is empty)
        % the second check is from numerical instabilities when GMRES ran
        % with such a high tolerence
        badII(end+1) = ii;
        continue;
    end
    %% Build cost function
    Acnt = @(x) Acounter(x,A, false);
    x0 = zeros(n,1);
    fg = @(x, Ax) quadraticLoss(x, Acnt, b, Ax);
    %% Fill the opts with problem specific information
    opts.A = Acnt;
    opts.b = b;
    xlb = A \ -b;
    try 
        [xstar,~] = callCVX(x0, A, b);
        opts.metricNames = {'abs_kkt', 'rel_kkt', 'MVP', 'rel_iter', 'abs_iter', 'obj'};
        opts.errFcn = {
            @(x) abs_kkt(x, A, b);
            @(x) rel_kkt(x, A, b);
            @(x) Acnt('cnt');
            @(x) rel_iter(x,xstar);
            @(x) abs_iter(x,xstar);
            @(x) dot(x,0.5*A*x+b);
            };
    catch 
        badII(end+1) = ii;
        continue
    end
    %% Run all the algorithms
    for ixAlgo = 1:numAlgo
        name = algoNames{ixAlgo};
        this_opts = opts;
        this_opts.name = name;
        this_opts.solver = algoSlvr{ixAlgo};
        [this_opts,~] = defaultLCPOpts(this_opts, x0);
        Acnt('reset');
        algo = algoHndls{ixAlgo};
        % disp(['--- ' name])
        tic
        [x, info] = algo(fg, x0, this_opts);
        results(i,ixAlgo).name = name;
        results(i,ixAlgo).x = x;
        results(i,ixAlgo).time = toc();
        results(i,ixAlgo).matVecs = Acnt('reset');
        results(i,ixAlgo).estimTime = results(i,ixAlgo).time + results(i,ixAlgo).matVecs*mean(dt);
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
    plotSubProblem(results(i,:), A, plotDebug, sprintf('i = %d', i) );
end
if ~isempty(badII)
    warning('did not identify all the good indices')
    mcGood(badII) = [];
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
        name(1:10), median(time), median(matVec), median(iters), median(kkt));
end
%% Plot Overall Statistics
iterationBarChart
%% Save results to File
saveFile = fullfile(root, [prefix '.mono.mat']);
disp(['Saving to ' saveFile]);
save(saveFile, ...
    'results', 'mcGood', 'opts');


