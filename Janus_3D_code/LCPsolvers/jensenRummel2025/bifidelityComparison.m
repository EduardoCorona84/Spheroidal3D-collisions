function bifidelityComparison(root, prefix, ps, tols, warmStart)
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
%% Get the indecies of the contact pairs
res_ = load(fullfile(root, [prefix '.mat']));
contactPairIX = cell(Nt,1);
for i = 1:Nt
    F = res_.lcp_list(i).F;
    if isempty(F)
        continue
    end
    disp("Processing contact pairs for problem " + i + " of " + Nt);
    n = size(F,2); % number of contact pairs
    N = size(F,1) / 6; % number of particles
    % For spheres torque is 0, so F will only be nonzero at the
    % positional places
    thesePairs = zeros(n,1);
    for ii = 1:n
        l0 = (find(F(:,ii), 1,'first')-1) / 6;
        l1 = (find(F(:,ii), 1,'last')-3) / 6;
        % linear indexing from 0
        % pair 0,1 -> 1, N,0 -> N(N-1), and so on
        thesePairs(ii) = l0*N + l1; 
    end
    contactPairIX{i} = thesePairs;
end
%% Hyper parameters
symmetrize = false;
plotDebug = false;
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
    'PQN';
    'B-PQN$\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)$';
    'B-PQN$\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)$';
    'B-PQN$\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)$';
    '(fixed) B-PQN$\bigl(\hat{\mathbf{A}}(p=3, \epsilon_{\mathrm{gmres}}=10^{-5})\bigr)$';
    '(fixed) B-PQN$\bigl(\hat{\mathbf{A}}(p=4, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)$';
    '(fixed) B-PQN$\bigl(\hat{\mathbf{A}}(p=6, \epsilon_{\mathrm{gmres}}=10^{-6})\bigr)$';
    };
algoSlvr = {
    'bbpgd';
    'proxquasinewton';
    'bifi';
    'bifi';
    'bifi';
    'bifi';
    'bifi';
    'bifi';
    };
algoHndls = {
    @projectedGradientDescent; 
    @proxQuasiNewton;
    @bifidelityProxQuasiNewton;
    @bifidelityProxQuasiNewton;
    @bifidelityProxQuasiNewton
    @bifidelityProxQuasiNewton;
    @bifidelityProxQuasiNewton;
    @bifidelityProxQuasiNewton
    };
biFi = {
    struct(); 
    struct();
    struct('p', 3, 'tol', 1e-5, 'fixed', false);
    struct('p', 4, 'tol', 1e-6, 'fixed', false);
    struct('p', 6, 'tol', 1e-6, 'fixed', false);
    struct('p', 3, 'tol', 1e-5, 'fixed', true);
    struct('p', 4, 'tol', 1e-6, 'fixed', true);
    struct('p', 6, 'tol', 1e-6, 'fixed', true);
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
    'matVecsLo', [], ...
    'eMatVecs', [], ...
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
    if symmetrize 
        A = (A +A')/2;
    end
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
    Acnt = @(x) Acounter(x,A,0);
    fg = @(x, Ax) quadraticLoss(x, Acnt, b, Ax);
    %% Fill the opts with problem specific information
    opts.A = Acnt;
    opts.b = b;
    xlb = A \ -b;
    try 
        [xstar,~] = callCVX(zeros(n,1), A, b);
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
        this_opts = opts;
        %% Warm start if possible
        x0 = zeros(n,1);
        if warmStart && i > 1 && ~isempty(results(i-1,ixAlgo).x)
            x_im1 = results(i-1,ixAlgo).x; 
            ix_im1 = contactPairIX{i-1};
            ix_i = contactPairIX{i};
            for ii = 1:n
                if any(ix_i(ii) == ix_im1)
                    jj = ix_i(ii) == ix_im1;
                    assert(sum(jj) == 1, 'We should never find two perfect matches...');
                    x0(ii) = x_im1(jj);
                end
            end
        end
        %% Set up bifi method 
        if isfield(biFi{ixAlgo},'p')
            global fixBifi 
            fixBifi = biFi{ixAlgo}.fixed;
            pLo = biFi{ixAlgo}.p;
            tolLo = biFi{ixAlgo}.tol;
            jLo = find(res.ps == pLo);
            kLo = find(abs(res.tols - tolLo)/abs(tolLo) < 1e-8);
            ALo = res.A{i,jLo,kLo};
            if symmetrize 
                ALo = (ALo +ALo')/2;
            end
            if isempty(ALo)
                badII(end+1) = ii;
                continue
            end
            ALoCnt = @(x) Acounter(x, ALo, 1);
            dtLo = res.dt{i,jLo,kLo};
            fctr = mean(dtLo) / mean(dt);
            this_opts.errFcn{3} = @(x) Acnt('cnt') + fctr * ALoCnt('cnt');
            this_opts.sub.solver = 'proxquasinewton';
            this_opts.low.A = @(x) ALoCnt(x);
            this_opts.low.kkt_abs = 1e-8;
            this_opts.low.kkt_rel = 1e-8;
            this_opts.low.b = b;
            this_opts.low.initWithLofi = ~warmStart;
        end
        %% Run Algo
        name = algoNames{ixAlgo};
        this_opts.name = name;
        this_opts.solver = algoSlvr{ixAlgo};
        [this_opts,~] = defaultLCPOpts(this_opts, x0);
        Acnt('reset');
        algo = algoHndls{ixAlgo};
        tic
        [x, info] = algo(fg, x0, this_opts);
        results(i,ixAlgo).name = name;
        results(i,ixAlgo).x = x;
        results(i,ixAlgo).time = toc();
        results(i,ixAlgo).matVecs = Acnt('cnt');
        results(i,ixAlgo).matVecsLo = 0;
        results(i,ixAlgo).eMatVecs = Acnt('cnt');
        results(i,ixAlgo).estimTime = results(i,ixAlgo).time + results(i,ixAlgo).eMatVecs*mean(dt);
        if isfield(biFi{ixAlgo},'p')
            matVecsLo = this_opts.low.A('cnt');
            lofiTime = matVecsLo * mean(dtLo);
            eMatVecs = mean(dtLo) / mean(dt) * matVecsLo;
            results(i,ixAlgo).matVecsLo = matVecsLo;
            results(i,ixAlgo).eMatVecs = results(i,ixAlgo).matVecs + eMatVecs;
        end
        results(i,ixAlgo).iters = info.iter;
        results(i,ixAlgo).kkt = info.kkt; 
        results(i,ixAlgo).errHist = [info.errHist results(i,ixAlgo).time + info.errHist(:,3)*mean(dt)] ;
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
fprintf('Algo            | Time      | eMatVec   | Iter | kkt\n')
for ixAlgo = 1:numAlgo
    name = algoNames{ixAlgo};
    while length(name) < length('Projected QuasiNewton (BFGS)')
        name = [name ' ']; %#ok<AGROW>
    end
    time = [results(mcGood,ixAlgo).time];
    iters = [results(mcGood,ixAlgo).iters];
    matVec = [results(mcGood,ixAlgo).eMatVecs];
    kkt = [results(mcGood,ixAlgo).kkt];
    fprintf('%s \t| %.1e s | %.3g\t| %.3g | %.2g\n', ...
        name(1:10), median(time), median(matVec), median(iters), median(kkt));
end
%% Plot Overall Statistics
iterationBarChart
%% Save results to File
if warmStart 
    saveFile = fullfile(root, [prefix '.warmstart.bifi.mat']);
else 
    saveFile = fullfile(root, [prefix '.bifi.mat']);
end
disp(['Saving to ' saveFile]);
save(saveFile, ...
    'results', 'mcGood',  'opts');


