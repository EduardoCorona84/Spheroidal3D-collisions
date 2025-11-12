[dirname, basedir] = setPaths();
srcFile = 'ampHi.special.n_3.p_8.cDist_2.3.allMats.mat';
res_ = load(fullfile(basedir, 'data.11.03.2025', srcFile));
A = res_.A;
dt = res_.dt;
ps = res_.ps;
tols = res_.tols;
%%
MC = size(A, 1);
numP = size(A, 2);
numTol = size(A,3);
[~,jHi] = max(res_.ps);
[~,kHi]= min(res_.tols);
i = find(arrayfun(@(i) ~isempty(A{i,jHi, kHi}) && ~isempty(A{i,jLo, kLo}), 1:MC), 1, 'last');
b = res_.b{i};
pHi = ps(jHi);
tolHi = tols(kHi);
AHi = A{i,jHi, kHi};
AHi = (AHi + AHi') /2;
dtHi = dt{i, jHi, kHi};
%% Timings pLot
% figure()
% f1 = figure;
% for j = 1:numP
%     for k = 1:numTol
%         subpLot(numP, numTol, numTol*(j-1) + k);
%         edges = 10.^(-6:.5:3);
%         [counts,edges] = Histcounts(dt{i,j,k},edges);
%         g = Histogram('BinEdges',edges,'BinCounts',counts);
%         set(gca, "Xscale", "Log")
%         xticks(edges)
%         title({ ...
%             sprintf('p = %d, tol = %.0e', ps(j), tols(k)), ...
%             sprintf('mean = %.3f',mean(dt{i,j,k})),...
%             sprintf('std = %.3f', std(dt{i,j,k}))}...
%         );
%     end
% end
% name = split(srcFile,'_');
% name = join(name,'\_');
% name = name{1};
% sgtitle({sprintf('Run Time for mc=%d', i),name})
% set(f1, 'Position',  [0, 0, 1000, 1200])
% timingFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_runTime.png']);
% disp(['Saving to ' timingFile])
% saveas(f1, timingFile);
%% Relative Error 
% relErr = zeros(numP, numTol);
% preCond = zeros(numP, numTol);
% timePerHi = zeros(numP, numTol);
% 
% for jLo = 1:numP
%     for kLo = 1:numTol
%         ALo = A{i,jLo, kLo};
%         ALo = (ALo + ALo') /2;
%         sqrtALoinv = inv(sqrtm(ALo));
%         relErr(jLo, kLo) = norm(AHi-ALo) / norm(AHi);
%         preCond(jLo, kLo) = cond(sqrtALoinv*AHi*sqrtALoinv);
%         timePerHi(jLo, kLo) = mean(dt{i,jHi,kHi}) / mean(dt{i,jLo,kLo});
%     end
% end
% name = split(srcFile,'_');
% name = join(name,'\_');
% name = name{1};
% 
% f2 = figure;
% subpLot(1,1,1)
% heatmap(ps, tols, relErr', 'CoLorLimits', [0 1e2])
% sgtitle(sprintf('Relative Error for mc=%d', i) )
% title(name)
% set(f2, 'Position',  [0, 0, 1000, 1200])
% relErrFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_relErr.png']);
% disp(['Saving to ' relErrFile])
% saveas(f2, relErrFile);
% 
% f3 = figure;
% subpLot(1,1,1)
% heatmap(ps, tols, preCond', 'CoLorLimits', [1 1e2])
% sgtitle(['Condition  Number of $\hat{A}^{-1/2}A\hat{A}^{-1/2}$ for ' sprintf('mc=%d', i)], 'Interpreter', 'latex');
% title(name)
% set(f3, 'Position',  [0, 0, 1000, 1200])
% timePerHiFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_preCond.png']);
% disp(['Saving to ' timePerHiFile])
% saveas(f3, timePerHiFile);
% 
% f4 = figure;
% subpLot(1,1,1)
% heatmap(ps, tols, timePerHi', 'CoLorLimits', [1 1e2])
% sgtitle(['Time to apply $\hat{A}$ vs $A$ for ' sprintf('mc=%d', i)], 'Interpreter', 'latex');
% title(name)
% set(f4, 'Position',  [0, 0, 1000, 1200])
% timePerHiFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_timePerHi.png']);
% disp(['Saving to ' timePerHiFile])
% saveas(f4, timePerHiFile);
%% Lo-Fidelity
jLo = 1;
[~,kLo]= max(res_.tols);
pLo = ps(jLo);
tolLo = tols(kLo);
ALo = A{i,jLo, kLo};
ALo = (ALo + ALo') /2;
sqrtALoinv = inv(sqrtm(ALo));
dtLo = dt{i, jLo, kLo};
%% Mid-Fidelity
jMid = find(ps == 5,1,'first');
kMid = find(tols ==1e-6, 1,'first');
pMid = ps(jMid);
tolMid = tols(kMid);
AMid = A{i,jMid, kMid};
AMid = (AMid + AMid') /2;
sqrtAMidinv = inv(sqrtm(AMid));
dtMid = dt{i, jMid, kMid};
%% Quick metrics
fprintf('Hi vs Lo RelErr = %.3g\n', norm(AHi - ALo) / norm(AHi))
fprintf('Cond of A_{Hi} = %.3g\n', cond(AHi))
fprintf('Cond of A_{Lo}^{-1} A_{Hi} = %.3g\n', cond(ALo \ AHi))
fprintf('Cond of A_{Lo}^{-1/2} A_{Hi} A_{Lo}^{-1/2} = %.3g\n', cond(sqrtALoinv*AHi*sqrtALoinv)) %#ok<*MINV>
fprintf('Hi vs Mid RelErr = %.3g\n', norm(AHi - AMid) / norm(AHi))
fprintf('Cond of A_{Mid}^{-1} A_{Hi} = %.3g\n', cond(AMid \ AHi))
fprintf('Cond of A_{Mid}^{-1/2} A_{Hi} A_{Mid}^{-1/2} = %.3g\n', cond(sqrtAMidinv*AHi*sqrtAMidinv))
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
    'PGD ($\kappa = \tau_{bb_1}, \eta = 1$)'; 'PGD ($\kappa = \tau_{bb_1}, \eta = \eta^*$)';
    'zeroSR1 ($\kappa = 1, \eta = \eta^*$)';
    'L-BFGS-B';
    'PQN (BFGS, $\kappa = 1, \eta = \eta^*$)';
    'Subspace Minimization (PQN)';
    'Min-Map Newton';
    ['PQN (BFGS, $B_0 = \hat{A}(' sprintf('p=%d', pLo) ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g', tolLo) '$)'];
    ['PQN (BFGS, $B_0 = \hat{A}(' sprintf('p=%d', pMid) ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g', tolMid) '$)']
    };
algoSlvr = {
    'pgd'; 'pgd';
    'zerosr1';
    'l-bfgs-b';
    'proxquasinewton';
    'subspacemin';
    'semismoothnewton';
    'proxquasinewton';
    'proxquasinewton';
    };
algoHndls = {
    @projectedGradientDescent; @projectedGradientDescent;
    @zeroSr1_nic;
    @L_BFGS_B;
    @proxQuasiNewton;
    @subspaceMin;
    @minmap_newton;
    @proxQuasiNewton;
    @proxQuasiNewton;
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
fg = @(x, Ax) quadraticLoss(x, Acnt, b, Ax);
%% Fill the opts with problem specific information
opts.A = Acnt;
opts.b = b;
xlb = AHi \ -b;
opts.fstar = 1/2*dot(xlb, AHi*xlb) + dot(xlb, b);
evals = eig(AHi);
opts.L = max(evals);
opts.mu = min(evals);
opts.AA = AHi;
mcGoodFlag = true;
try %#ok<TRYNC>
    [xstar,~] = callCVX([], x0, opts);
    opts.metricNames = {'abs_kkt', 'rel_kkt', 'MVP', 'rel_iter', 'abs_iter', 'obj'};
    opts.errFcn = {
        @(x) abs_kkt(x, AHi, b);
        @(x) rel_kkt(x, AHi, b);
        @(x) Acnt('cnt');
        @(x) rel_iter(x,xstar);
        @(x) abs_iter(x,xstar);
        @(x) dot(x,0.5*AHi*x+b);
        };
    mcGood(mc) = true;
end
%% Run all the algorithms
Acnt('reset')
for ixAlgo = 1:numAlgo
    tHis_opts = opts;
    name = algoNames{ixAlgo};
    if contains(name, 'B_0') 
        if contains(name, num2str(pLo))
            tHis_opts.B0 = @(x) ALo*x;
            L = chol(ALo);
            tHis_opts.H0 = @(x) L \ (L' \ x);
            tHis_opts.prox_B0 = @(y) solveGeneralProx(y, ALo);
        elseif contains(name, num2str(pMid))
            tHis_opts.B0 = @(x) AMid*x;
            L = chol(AMid);
            tHis_opts.H0 = @(x) L \ (L' \ x);
            tHis_opts.prox_B0 = @(y) solveGeneralProx(y, AMid);
        end

    end
    tHis_opts.name = name; 
    tHis_opts.solver = algoSlvr{ixAlgo};
    [tHis_opts,~] = defaultLCPOpts(tHis_opts, x0);
    algo = algoHndls{ixAlgo};
    tic
    [x, info] = algo(fg, x0, tHis_opts);
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
plotSubProblem(results, AHi,true, 'Special Config')