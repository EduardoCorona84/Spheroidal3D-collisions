function results = testLCPWrappers()
%% Hyper parameters
plotDebug = false;
fromFile = false;
diagDom = false;
minSz = 100;
maxSz = 150;
condNum = 1e3;
percentLarge = 0.1;
MC = 100;
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
algoNames = {
    'PGD (\kappa = \tau_{bb_1}, \eta = 1)'; 'PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    % 'Optimal Diagonal Preconditioned PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    % 'Online Preconditioned PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    'zeroSR1 (\kappa = 1, \eta = \eta^*)';
    'L-BFGS-B';
    % 'Projected QuasiNewton (BFGS, \kappa = 1, \eta = \eta^*)';
    'PQN (BFGS, \kappa = 1, \eta = \eta^*)';
    'Subspace Minimization (PQN)';% 'Subspace Minimization (PGD)';
    % 'Binding Proximal QuasiNewton (BFGS, \kappa = 1, \eta = \eta^*)'
    };
algoHndls = {
    @projectedGradientDescent; @projectedGradientDescent;
    % @optDiagPrecond;
    % @onlineScaledGradient;
    @zeroSr1_nic;
    @L_BFGS_B;
    % @projectQuasiNewton_nic;
    @proxQuasiNewton;
    @subspaceMin; %@subspaceMin;
    % @bindingProxQuasiNewton
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
    ), [MC,numAlgo] ...
    );
%% Switch between random problems and those save to a file
if fromFile
    fname = 'all_data'; %#ok<*UNRCH>
    load([fname '.mat'], ...
        'A_list', 'b_list');
    MC = length(A_list);
else
    try %#ok<TRYNC>
        rng('default')
    end
    rng(1);
end
mcGood = false(MC,1);
for mc = 1:MC
    disp(['mc = ' num2str(mc) '/' num2str(MC)])
    if fromFile
        %% Load problem from list
        A = A_list{mc};
        A = 1/2*(A +A');
        b = b_list{mc};
        n = size(A,2);
        if n < minSz || maxSz <= n
            continue;
        end
    else
        %% Generate random problem
        n = minSz + int64(round((maxSz-minSz)*rand(1)));
        B = randn(n,n);
        numLarge = round(percentLarge*n);
        vv = [1; 1 + rand(n-numLarge-1,1); 1 + condNum * rand(numLarge-1,1); condNum];
        if diagDom
            A = (B+B') + diag(vv);
        else
            [Q, ~] = qr(B);
            A = Q*diag(vv)*Q';
        end
        [vecs,vals] = eig(A);
        vals = max(1,diag(vals));
        A = vecs*diag(vals)*vecs';
        x_unconstrained = randn(n,1);
        while ~any(x_unconstrained<0)
            x_unconstrained = randn(n,1);
        end
        b = -A*x_unconstrained;
    end
    %% Build cost function
    Acnt = @(x) Acounter(x,A, false);
    x0 = zeros(n,1);
    fg = @(x, Ax) quadraticLoss(x, Acnt, b, Ax);
    %% Fill the opts with problem specific information
    opts.A = Acnt;
    opts.b = b;
    xlb = A \ -b;
    opts.fstar = 1/2*dot(xlb, A*xlb) + dot(xlb, b);
    evals = eig(A);
    opts.L = max(evals);
    opts.mu = min(evals);
    opts.AA = A;
    mcGoodFlag = true;
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
        mcGood(mc) = true;
    end
    %% Run all the algorithms
    for ixAlgo = 1:numAlgo
        if ~mcGoodFlag
            break
        end
        name = algoNames{ixAlgo};
        [this_opts,~] = defaultLCPOpts(opts, x0);
        if contains(name, 'kappa') && contains(name, 'bb')
            this_opts.stepSize.kappa = 'bb1';
        elseif contains(name, 'kappa = 1')
            this_opts.stepSize.kappa = 'uniform';
        elseif contains(name, 'kappa^*')
        elseif contains(lower(name), 'subspace')
            if contains(lower(name), 'pgd')
                this_opts.subspaceMin.innerStepSelection = 'pgd';
                this_opts.stepSize.kappa = 'opt';
            else
                this_opts.subspaceMin.innerStepSelection = 'pqn';
                this_opts.stepSize.kappa = 'uniform';
            end
        else
            this_opts.stepSize.kappa = 'uniform';
        end
        if contains(name, 'eta = 1')
            this_opts.stepSize.eta = 'uniform';
        elseif contains(name, 'eta^*')
            this_opts.stepSize.eta = 'opt';
        end
        algo = algoHndls{ixAlgo};
        tic
        [x, info] = algo(fg, x0, this_opts);
        results(mc,ixAlgo).name = name;
        results(mc,ixAlgo).x = x;
        results(mc,ixAlgo).time = toc();
        results(mc,ixAlgo).iters = info.iter;
        results(mc,ixAlgo).kkt = info.kkt;
        results(mc,ixAlgo).matVecs = Acnt('reset');
        results(mc,ixAlgo).errHist = info.errHist;
        try %#ok<TRYNC>
            results(mc,ixAlgo).iterHist{mc} = info.iterHist;
        end
        for ixMetric = 1:numel(opts.errFcn)
            hndl = opts.errFcn{ixMetric};
            try %#ok<TRYNC>
                hndl('reset');
            end
        end
    end
    %% Debug Plot
    plotSubProblem(results(mc,:), A, vals, mcGood(mc) && plotDebug);
end
%% Print debug information
fprintf('Algo                         | Time      | matVec | Iter | kkt\n')
for ixAlgo = 1:numAlgo
    name = algoNames{ixAlgo};
    while length(name) < length('Projected QuasiNewton (BFGS)')
        name = [name ' ']; %#ok<AGROW>
    end
    time = results(ixAlgo).time(mcGood);
    iters = results(ixAlgo).iters(mcGood);
    matVec = results(ixAlgo).matVecs(mcGood);
    kkt = results(ixAlgo).kkt(mcGood);
    fprintf('%s | %.1e s | %.3g\t| %.3g | %.2g\n', ...
        name, mean(time), mean(matVec), mean(iters), mean(kkt));
end
%% Plot Overall Statistics
iterationBarChart
%% Save results to File
if fromFile
    fname = ['results_n_' num2str(minSz) '_' num2str(maxSz)];
elseif diagDom
    fname = ['results_randProblems_diagDom_n_' num2str(minSz) '_' num2str(maxSz)];
elseif percentLarge ~= 1
    fname = ['10.30.2025.results_randProblems_percentLarge_' num2str(percentLarge) '_n_' num2str(minSz) '_' num2str(maxSz)];
else
    fname = ['results_randProblems_noStructure_updateStepSize_n_' num2str(minSz) '_' num2str(maxSz)];
end
tmp = mcGood;
mcGood = find(tmp);
disp(['Saving to ' fname]);
save(['/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/' fname '.mat'], ...
    'results', 'mcGood', 'minSz', 'maxSz', 'opts');
end


