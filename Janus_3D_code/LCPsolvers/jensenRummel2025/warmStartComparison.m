%% Set paths
[dirname, basedir] = setPaths();
%% Try to initialize cvx if on the cluster
try %#ok<TRYNC>
    run(fullfile('/projects', getenv('USER'), 'cvx', 'cvx_startup.m'))
end
%% Load Data from file
resFile = '/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/goodData/amphi.lattice.n_5.p_8.cDist_2.5.allMats.mat';
res = load(resFile);
[~,jHi] = max(res.ps);
[~,kHi]= min(res.tols);
Nt = size(res.A,1);
%% Get the indecies of the contact pairs
res_ = load('/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/data.11.03.2025/lcp.amphi.lattice.n_5.p_8.cDist_2.5.mat');
contactPairIX = cell(Nt,1);
for i = 1:Nt
    F = res_.lcp_list(i).F;
    if isempty(F)
        continue
    end
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
prefix = 'amphi.lattice.n_5.warmStart';
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
    'Accelerated PGD'
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
mcGood = getMCGood(resFile);
badII = [];
for ii = 1:numel(mcGood)
    i = mcGood(ii);
    disp(['i = ' num2str(i) '/' num2str(i)])
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
        %% Warm start if possible
        x0 = zeros(n,1);
        if i > 1 && ~isempty(results(i-1,ixAlgo).x)
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
        name = algoNames{ixAlgo};
        this_opts = opts;
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
        name(1:10), mean(time), mean(matVec), mean(iters), mean(kkt));
end
%% Plot Overall Statistics
iterationBarChart
%% Save results to File
disp(['Saving to ' prefix]);
save(['/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/' prefix '.mat'], ...
    'results', 'mcGood',  'opts');


