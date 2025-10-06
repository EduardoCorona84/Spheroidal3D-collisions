function results = testLCPWrappers()
try %#ok<TRYNC>
    rng('default')
end
rng(2);
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname, ~,~] = fileparts(mfilePath);
addpath(genpath(fileparts(dirname)))
fname = 'all_data';
load([fname '.mat'], ...
    'A_list', 'b_list');
minSz = 100;
maxSz = 150;

opts = struct( ...
    'max_iter',50, ...
    'kkt_rel',1e-12, ...
    'kkt_abs',1e-12, ...
    'arg_rel',[], ...
    'arg_abs',[], ...
    'step_abs',[], ...
    'stepSize',struct('init','uniform',...
        'kappa','uniform',...
        'eta','opt'),...
    'prox',struct( ...
        'maxiter', 1000, ...
        'res_abstol', 0, ...
        'res_reltol', 0, ...
        'alp_abstol', 0, ...
        'alp_reltol', 1e-12, ...
        'verbose', false, ...
        'runCVX', false ...
    ), ...
    'r', 20, ...
    'qnUpdate', 'bfgs', ...
    'storeIts', true, ...
    'subSpaceMin', struct( ...'
        'innerStepSelection', 'pqn', ...
        'innerSolver','cvx', ...
        'orthoMethod','qr', ...
        'kThresh', 1, ...
        'useOneStepIter', false ...
    ) ...
);
algoNames = {
    % 'PGD (\kappa = \tau_{bb_1}, \eta = 1)';
    'PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    % 'Optimal Diagonal Preconditioned PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    % 'Online Preconditioned PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    % 'zeroSR1 (\kappa = 1, \eta = \eta^*)';
    % 'L-BFGS-B';
    % 'Projected QuasiNewton (BFGS, \kappa = 1, \eta = \eta^*)';
    'PQN (BFGS, \kappa = 1, \eta = \eta^*)';
    'Subspace Minimization (PGD)';
    'Subspace Minimization (PQN)';
    % 'Binding Proximal QuasiNewton (BFGS, \kappa = 1, \eta = \eta^*)'
    };
algoHndls = {
    @projectedGradientDescent; %@projectedGradientDescent;
    % @optDiagPrecond;
    % @onlineScaledGradient;
    % @zeroSr1_nic;
    % @L_BFGS_B;
    % @projectQuasiNewton_nic;
    @proxQuasiNewton;%
    @subspaceMin;
    @subspaceMin;
    % @bindingProxQuasiNewton
    };
numAlgo = numel(algoNames);
assert(numel(algoNames) == numel(algoHndls));
fromFile = false;
diagDom = false;
if fromFile
    MC = length(A_list);
else
    MC = 10;
end
results = repmat(...
    struct( ...
    'algo', '',...
    'time',   zeros(MC,1), ...
    'iters', zeros(MC,1), ...
    'kkt', zeros(MC,1), ...
    'matVecs', zeros(MC,1), ...
    'errHist', {cell(MC,1)}, ...
    'iterHist', {cell(MC,1)} ...
    ), [1,numAlgo] ...
    );
mcGood = [];
for mc = 1:MC
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
        condNum = 1e2;
        n = minSz + int64(round((maxSz-minSz)*rand(1)));
        numLarge = 50;
        B = randn(n,n);
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
        while ~any(x_unconstrained)
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
    try
        [xstar,~] = callCVX([], x0, opts);
        opts.errFcn = {
            @(x) abs_kkt(x, A, b);
            @(x) rel_kkt(x, A, b);
            @(x) Acnt('cnt');
            @(x) rel_iter(x,xstar);
            @(x) abs_iter(x,xstar);
            };
    catch
        mcGoodFlag = false;
    end
    %% Run all the algorithms
    for ixAlgo = 1:numAlgo
        if ~mcGoodFlag
            break
        end
        name = algoNames{ixAlgo};
        this_opts = opts;
        if contains(name, 'kappa') && contains(name, 'bb')
            this_opts.stepSize.kappa = 'bb1';
        elseif contains(name, 'kappa = 1')
            this_opts.stepSize.kappa = 'uniform';
        elseif contains(name, 'kappa^*')
        elseif contains(lower(name), 'subspace') 
            if contains(lower(name), 'pgd')
                this_opts.subSpaceMin.innerStepSelection = 'pgd';
                this_opts.stepSize.kappa = 'opt';
            else
                this_opts.subSpaceMin.innerStepSelection = 'pqn';
                this_opts.stepSize.kappa = 'uniform';
            end
        else    
             this_opts.stepSize.kappa = 'uniform';
        end
        if contains(name, 'eta = 1')
            this_opts.stepSize.eta = 'uniform';
        elseif strcmpi(name, 'eta^*')
            this_opts.stepSize.eta = 'opt';
        end
        algo = algoHndls{ixAlgo};
        tic
        [x, info] = algo(fg, x0, this_opts);
        results(ixAlgo).name = name;
        results(ixAlgo).time(mc) = toc();
        results(ixAlgo).iters(mc) = info.iter;
        results(ixAlgo).kkt(mc) = info.kkt;
        results(ixAlgo).matVecs(mc) = Acnt('reset');
        results(ixAlgo).errHist{mc} = info.errHist;
        results(ixAlgo).iterHist{mc} = info.iterHist;
        for ixMetric = 1:numel(opts.errFcn)
            hndl = opts.errFcn{ixMetric};
            try %#ok<TRYNC>
                hndl('reset');
            end
        end
    end
    if mcGoodFlag
        mcGood = [mcGood mc]; %#ok<AGROW>
        disp(['mc ' num2str(mc)])
        linespec = {"-o", "--s",":*","-.diamond"};
        
        objVal = cell(numAlgo,1);
        for ixAlgo = 1:numAlgo
            iterHist = results(ixAlgo).iterHist{mc};
            errHist = results(ixAlgo).errHist{mc};
            numIter = size(errHist,1)-1;
            try
                objVal{ixAlgo} = arrayfun(@(i) fg(iterHist(i,:)',[]), 1:numIter+1);
            catch
                objVal{ixAlgo} = errHist(:,end);
            end
        end
        minObjVal = min(cellfun(@min, objVal));
        figure()
        for ixAlgo = 1:numAlgo
            name = algoNames{ixAlgo};
            subplot(3,1,1)
            errHist = results(ixAlgo).errHist{mc};
            numIter = size(errHist,1)-1;
            % matVecs = errHist(:,3);
            kk = mod(ixAlgo-1, numel(linespec)) + 1;
            semilogy(0:numIter, objVal{ixAlgo} - minObjVal + 1e-12, linespec{kk}, 'LineWidth',4,'MarkerSize',10) % abs(objVal - cvxObjVal) / abs(cvxObjVal))
            hold on
            xlabel('Iteration','FontSize', 20)
            ylabel('$f(x_k)$','interpreter', 'latex', 'FontSize', 25)
            subplot(3,2,3)
            semilogy(0:numIter, errHist(:,4), linespec{kk},'LineWidth',4,'MarkerSize',10);
            hold on
            ylabel('$\frac{\|x_k - x^*\|}{\|x^*\|}$', 'interpreter', 'latex', 'FontSize', 30)
            xlabel('Iteration','FontSize', 20)
            subplot(3,1,3)
            semilogy(0:numIter, errHist(:,1), linespec{kk}, 'LineWidth',4, 'MarkerSize',10);
            hold on
            xlabel('Iteration','FontSize', 20)
            ylabel('$\varphi(x_k)$','interpreter', 'latex','FontSize', 25)
        end
        subplot(3,2,4)
        plot(sort(vals,1,'descend'), '-o', 'Color', '#808080', ...
            'LineWidth',4, 'MarkerSize',10)
        ylabel('Eigen Values of $A$', 'interpreter', 'latex','FontSize', 25)
        subplot(3,1,1)
        legend(algoNames,'FontSize', 20,'Location','northeastoutside')
        sgtitle({['Deep Dive for MC = ' num2str(mc)], ...
            sprintf('n = %d, n^\\prime = %d, \\kappa(A) = %.2g', ...
            n, numLarge, cond(A))},'FontSize', 30)
    end
end
%%
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
iterationBarChart
if fromFile
    fname = ['results_n_' num2str(minSz) '_' num2str(maxSz)];
elseif diagDom
    fname = ['results_randProblems_diagDom_n_' num2str(minSz) '_' num2str(maxSz)];
else
    fname = ['results_randProblems_noStructure_updateStepSize_n_' num2str(minSz) '_' num2str(maxSz)];
end
save(['/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/' fname '.mat'], ...
    'results', 'mcGood', 'minSz', 'maxSz', 'opts');
end

function Ax = Acounter(x, A, transpose)
persistent matVecCnt

if isempty(matVecCnt)
    matVecCnt = 0;
end
if ~exist('transpose','var') || isempty(transpose)
    transpose = false;
end

if ischar(x)
    if strcmpi(x, 'cnt')
        Ax = matVecCnt;
        return
    elseif strcmpi(x, 'reset')
        Ax = matVecCnt;
        matVecCnt = 0;
        return
    else
        assert(false, ['Option: ' x ' not recognized'])
    end
end
if transpose
    Ax = A'*x;
else
    Ax = A*x;
end
matVecCnt = matVecCnt + 1;
end

function [x, info] = callCVX(~, x0, opts) %#ok<STOUT>
n = length(x0); %#ok<NASGU>
% f = @(x) 1/2*sum_square(opts.A(x)+opts.b);
f = @(x) 1/2*dot(x,opts.A(x)) + dot(opts.b,x); %#ok<NASGU>
cvx_begin quiet
variable x(n)
minimize f(x)
subject to
0 <= x %#ok<NODEF,NOPRT>
cvx_end

info.iter = NaN;
info.kkt = NaN;
info.errHist = NaN;
info.iterHist = NaN;
end

function e = abs_kkt(x, A, b)
phi = min(x,A*x + b);
e = 1/2*dot(phi, phi);
end % abs_kkt

function diff = rel_kkt(x, A, b)
persistent olde
if ischar(x) && strcmpi(x,'reset')
    olde = [];
    diff = NaN;
    return
end
phi = min(x,A*x + b);
e = 1/2*dot(phi, phi);
if isempty(olde)
    diff = NaN;
else
    diff = abs(e - olde) / abs(e);
end
olde = e;
end % rel_kkt

function e = abs_iter(x, xstar)
e = norm(x - xstar);
end % abs_iter

function e = rel_iter(x, xstar)
e = norm(x-xstar) / norm(xstar);
end % rel_iter