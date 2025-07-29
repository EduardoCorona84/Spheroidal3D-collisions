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
minSz = 50; 
maxSz = 100;

opts = struct( ...
    'max_iter',1000, ...
    'tol_rel',1e-12, ...
    'tol_abs',1e-12, ... 
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
        'runCVX', false), ...
    'r', 20, ...
    'qnUpdate', 'bfgs', ...
    'storeIts', true...
);
MC = length(A_list);

algoNames = {
    % 'CVX';
    % 'PGD (\kappa = \tau_{bb_1}, \eta = 1)'; 
    % 'PGD (\kappa = \kappa^*, \eta = \eta^*)'; 
    'PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)'; 
    % 'zeroSR1 (\kappa = 1, \eta = \eta^*)';
    % 'Proximal QuasiNewton (BFGS, \kappa = \kappa^*, \eta = \eta^*)';
    'Proximal QuasiNewton (BFGS, \kappa = 1, \eta = \eta^*)';
    'Optimal Preconditioned PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    'Online Preconditioned PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    % 'L-BFGS-B';
    % 'Projected QuasiNewton (BFGS)';
    % 'Binding Proximal QuasiNewton (BFGS)'
};
algoHndls = {
    % @callCVX, 
    @projectedGradientDescent;% @projectedGradientDescent; @projectedGradientDescent;
    % @zeroSr1_nic; 
    @proxQuasiNewton;% @proxQuasiNewton; 
    @optDiagPrecond;
    @onlineScaledGradient;
    % @L_BFGS_B; 
    % @projectQuasiNewton_nic;  
    % @bindingProxQuasiNewton
};
numAlgo = numel(algoNames);
assert(numel(algoNames) == numel(algoHndls));
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
MC = numel(A_list);
fromFile =false ;
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
        condNum = 1e4;
        n = minSz + int64(round((maxSz-minSz)*rand(1)));
        B = randn(n,n);
        % [Q, ~] = qr(B);
        vv = 1 + condNum * rand(n,1);
        % A = Q*diag(vv)*Q';
        A = (B+B') + diag(vv);
        [vecs,vals] = eig(A);
        vals = max(1,diag(vals));
        A = vecs*diag(vals)*vecs';
        x_unconstrained = randn(n,1);
        while ~any(x_unconstrained < 0)
            x_unconstrained = randn(n,1);
        end
        b = -A*x_unconstrained;
    end
    %% Build cost function 
    Acnt = @(x) Acounter(x,A, false);
    x0 = zeros(n,1);
    fg = @(x, Ax, Aq, eta) quadraticLoss(x, Acnt,b, Ax, Aq, eta);
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
            this_opts.stepSize.kappa = 'opt';
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
        figure() 
        for ixAlgo = 1:numAlgo
            name = algoNames{ixAlgo};
            subplot(3,1,1)
            iterHist = results(ixAlgo).iterHist{mc};
            errHist = results(ixAlgo).errHist{mc};
            numIter = size(errHist,1);
            try
                objVal = arrayfun(@(i) fg(iterHist(i,:)',[],[],[]), 1:numIter);
            catch 
                objVal = errHist(:,end);
            end
            semilogy(1:numIter, objVal - min(objVal) + 1e-12, 'LineWidth',4) % abs(objVal - cvxObjVal) / abs(cvxObjVal))
            hold on 
            ylabel('$f(x_k)$','interpreter', 'latex', 'FontSize', 25)
            subplot(3,1,2)
            semilogy(1:numIter, errHist(:,4), 'LineWidth',4);
            hold on
            ylabel('$\frac{\|x_k - x^*\|}{\|x^*\|}$', 'interpreter', 'latex', 'FontSize', 30)
            subplot(3,1,3)
            semilogy(1:numIter, errHist(:,1), 'LineWidth',4);
            hold on 
            xlabel('Iteration','FontSize', 20)
            ylabel('$\varphi(x_k)$','interpreter', 'latex','FontSize', 25)
        end
        subplot(3,1,1)
        legend(algoNames,'FontSize', 20,'Location','northeastoutside')
        sgtitle(['Iteration Metrics for MC = ' num2str(mc)],'FontSize', 30)
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
fname = ['results_randProblems_n_' num2str(minSz) '_' num2str(maxSz)];
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