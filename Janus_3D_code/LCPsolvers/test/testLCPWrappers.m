function results = testLCPWrappers()
%% Hyper parameters
prefix = 'amphi.lattice.n_5.p_8.cDist_2.5.warmStart';
plotDebug = true;
warmStart = false;
fromFile = false;
diagDom = false; %#ok<*NASGU>
minSz = 100;
maxSz = 200;
condNum = 1e2;
percentLarge = 0.25;
lid = 1;
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
    'storeIts', true);
%% Switch between random problems and those save to a file
if fromFile %#ok<*UNRCH>
    res_ = load('/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/data.11.03.2025/amphi.lattice.n_5.p_8.cDist_2.5.allMats.mat');
    [~,j] = max(res_.ps);
    j = j - 1; % p =8 is not all done yet
    [~,k]= min(res_.tols);
    MC = size(res_.A,1);
    A_list = arrayfun(@(i) (res_.A{i,j,k} + res_.A{i,j,k}')/2, 1:MC, 'UniformOutput',false);
    dt_list = arrayfun(@(i) res_.dt{i,j,k}, 1:MC, 'UniformOutput',false);
    jLo = 1;
    [~,kLo]= max(res_.tols);
    pLo = res_.ps(jLo);
    tolLo = res_.tols(kLo);
    ALo_list = arrayfun(@(i) (res_.A{i,jLo,kLo} + res_.A{i,jLo,kLo}') / 2, 1:MC, 'UniformOutput',false);
    dtLo_list = arrayfun(@(i) res_.dt{i,jLo,kLo}, 1:MC, 'UniformOutput',false);
    jMid = find(res_.ps == 5,1,'first');
    kMid = find(res_.tols == 1e-8, 1,'first');
    pMid = res_.ps(jMid);
    tolMid = res_.tols(kMid);
    AMid_list = arrayfun(@(i) (res_.A{i,jMid,kMid} + res_.A{i,jMid,kMid}') / 2, 1:MC, 'UniformOutput',false);
    dtMid_list = arrayfun(@(i) res_.dt, 1:MC, 'UniformOutput',false);
    b_list = arrayfun(@(i) res_.b{i}, 1:MC, 'UniformOutput',false);
    res_ = load('/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/data.11.03.2025/lcp.amphi.lattice.n_5.p_8.cDist_2.5.mat');
    F_list = arrayfun(@(i) res_.lcp_list(i).F, 1:MC, 'UniformOutput',false);
    lid = find(arrayfun(@(i) ~isempty(A_list{i}) && ~isempty(F_list{i}) , 1:MC), 1, 'first');
    MC = find(arrayfun(@(i) ~isempty(A_list{i}) && ~isempty(F_list{i}) , 1:MC), 1, 'last');
    pairs = cell(MC,1);
    for i = lid:MC
        A = A_list{i};
        F = F_list{i};
        if isempty(A) || isempty(F)
            continue
        end
        n = size(A,1);
        N = size(F,1) / 6;
        pairs{i} = zeros(n,2);
        % For spheres torque is 0, so F will only be nonzero at the
        % positional places
        pairs{i}(:,1) = arrayfun(@(ii) (find(F(:,ii), 1,'first')-1) / 6 + 1, 1:n);
        pairs{i}(:,2) = arrayfun(@(ii) (find(F(:,ii), 1,'last')-3) / 6 + 1, 1:n);
    end
else
    try %#ok<TRYNC>
        rng('default')
    end
    rng(1);
end
algoNames = {
    'PGD (\tau = \tau_{bb_1}, \eta = 1)'; 'PGD (\tau = \tau_{bb_1}, \eta = \eta^*)';
    % 'Optimal Diagonal Preconditioned PGD (\tau = \tau_{bb_1}, \eta = \eta^*)';
    % 'Online Preconditioned PGD (\tau = \tau_{bb_1}, \eta = \eta^*)';
    'zeroSR1 (\tau = 1, \eta = \eta^*)';
    'L-BFGS-B';
    % 'Projected QuasiNewton (BFGS, \tau = 1, \eta = \eta^*)';
    'PQN (BFGS, \tau = 1, \eta = \eta^*)';
    'Subspace Minimization (PQN)';% 'Subspace Minimization (PGD)';
    % 'Binding Proximal QuasiNewton (BFGS, \tau = 1, \eta = \eta^*)'
    'Min-Map Newton';
    % ['PQN (BFGS, $B_0 = \hat{A}(' sprintf('p=%d', pLo) ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g', tolLo) '$)'];
    % ['PQN (BFGS, $B_0 = \hat{A}(' sprintf('p=%d', pMid) ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g', tolMid) '$)'];
    % ['PQN (BFGS, $B_0 = B(' sprintf('p=%d', pMid) ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g', tolMid) '$)'];
    };
algoSlvr = {
    'pgd'; 'pgd';
    'zerosr1';
    'l-bfgs-b';
    'proxquasinewton';
    'subspacemin';
    'semismoothnewton';
    % 'proxquasinewton';
    % 'proxquasinewton';
    % 'proxquasinewton';
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
    @minmap_newton;
    % @proxQuasiNewton;
    % @proxQuasiNewton;
    % @proxQuasiNewton;
    % @bindingProxQuasiNewton
    };
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
    ), [MC,numAlgo] ...
    );
mcGood = false(MC,1);
for mc = lid:MC
    disp(['mc = ' num2str(mc) '/' num2str(MC)])
    if fromFile
        %% Load problem from list
        A = A_list{mc};
        dt = dt_list{mc};
        ALo = ALo_list{mc};
        AMid = AMid_list{mc};
        dtMid = dt_list{i};
        if isempty(ALo) || isempty(AMid)
            continue
        end
        A = 1/2*(A +A');
        b = b_list{mc};
        n = size(A,2);
        if isempty(A) || n < minSz || maxSz <= n
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
    if warmStart && mc > 1 && ~isempty(results(mc-1,1).x)
        ix_im1 = pairs{mc-1};
        ix_i = pairs{mc};
        for ii = 1:n
            n_im1 = size(ix_im1,1);
            jj = find(arrayfun(@(jj) all(ix_i(ii,:) == ix_im1(jj,:)), 1:min(n, n_im1)));
            if ~isempty(jj)
                assert(isscalar(jj), 'We should never find two perfect matches...');
                x0(ii) = results(mc-1,1).x(jj);
            end
        end
    end
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
        if ~mcGood(mc)
            break
        end
        name = algoNames{ixAlgo};
        this_opts= opts;
        this_opts.name = name;
        this_opts.solver = algoSlvr{ixAlgo};
        [this_opts,~] = defaultLCPOpts(this_opts, x0);
        results(mc,ixAlgo).estimTime = 0;
        if contains(name, 'B_0') 
            if contains(name, num2str(pLo))
                this_opts.B0 = @(x) ALo*x;
                L = chol(ALo);
                this_opts.H0 = @(x) L \ (L' \ x);
                this_opts.prox_B0 = @(y) solveGeneralProx(y, ALo);
            elseif contains(name, num2str(pMid)) && contains(name, '\hat{A}')
                L = chol(AMid);
                [U, lambdas] = eig(AMid);
                b0 = min(diag(lambdas)) / 2;
                U = U*sqrtm(diag(diag(lambdas)-b0));
                relErr = norm(b0*eye(n) + U*U' - AMid) / norm(AMid);
                this_opts.B0 = @(x) b0*x;
                this_opts.H00 = @(x) x ./ b0;
                this_opts.U0 = U;
                this_opts.H0 = @(x) L \ (L' \ x);
            else
                AMidCnt = @(x) Acounter(x, AMid, false);
                mid_opts = this_opts;
                mid_opts.A = AMidCnt; 
                fgMid =  @(x, Ax) quadraticLoss(x, AMidCnt, b, Ax);
                [x0, ~] = proxQuasiNewton(fgMid, x0, mid_opts);
                this_opts.useCurrentB = true;
                results(mc,ixAlgo).AhatMVP = AMidCnt('reset');
                results(mc,ixAlgo).estimTime = results(mc,ixAlgo).AhatMVP*mean(dtMid);
                for ixMetric = 1:numel(opts.errFcn)
                    hndl = opts.errFcn{ixMetric};
                    try %#ok<TRYNC>
                        hndl('reset');
                    end
                end
            end
        end
        Acnt('reset');
        algo = algoHndls{ixAlgo};
        tic
        [x, info] = algo(fg, x0, this_opts);
        results(mc,ixAlgo).name = name;
        results(mc,ixAlgo).x = x;
        results(mc,ixAlgo).time = toc();
        results(mc,ixAlgo).matVecs = Acnt('reset');
        % results(mc,ixAlgo).estimTime = results(mc,ixAlgo).estimTime + results(mc,ixAlgo).matVecs*mean(dt);
        results(mc,ixAlgo).iters = info.iter;
        results(mc,ixAlgo).kkt = info.kkt; 
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
    plotSubProblem(results(mc,:), A, mcGood(mc) && plotDebug, sprintf('mc = %d', mc) );
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


