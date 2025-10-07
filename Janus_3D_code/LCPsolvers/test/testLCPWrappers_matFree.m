function results = testLCPWrappers()
%% Set Up Matlab Path
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
%% Set the load and save paths
root = fullfile(dirname, '../data');
loadFile = 'amphiLCPs.n_2.p_8.cDist_2.3.prt_0.mat';
saveFile = 'results.n_2.p_8.cDist_2.3.mat';
%% Load scenario data from file 
load(fullfile(root, loadFile), ...
    'Fparams', 'lcp_list');
%% Set Hyperparameters
lcpOpts = struct( ...
    'max_iter',50, ...
    'kkt_rel',1e-12, ...
    'kkt_abs',1e-12, ...
    'storeIts', true, ...
    'subSpaceMin', struct( ...
        'innerStepSelection', 'pqn', ...
        'innerSolver','cvx', ...
        'orthoMethod','qr', ...
        'kThresh', 1, ...
        'useOneStepIter', false ...
    ) ...
);
%% Specify LCP solvers
algoNames = {
    % 'PGD (\kappa = \tau_{bb_1}, \eta = 1)';
%     'PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    % 'Optimal Diagonal Preconditioned PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    % 'Online Preconditioned PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)';
    % 'zeroSR1 (\kappa = 1, \eta = \eta^*)';
    % 'L-BFGS-B';
    % 'Projected QuasiNewton (BFGS, \kappa = 1, \eta = \eta^*)';
    'Subspace Minimization (PQN)';
    'PQN (BFGS, \kappa = 1, \eta = \eta^*)';
%     'Subspace Minimization (PGD)';
    % 'Binding Proximal QuasiNewton (BFGS, \kappa = 1, \eta = \eta^*)'
    };
algoHndls = {
%     @projectedGradientDescent; 
    % @zeroSr1_nic;
    % @L_BFGS_B;
    @subspaceMin;
    @proxQuasiNewton;%
%     @subspaceMin;

};
numAlgo = numel(algoNames);
% idiot check
assert(numel(algoNames) == numel(algoHndls), 'Must specify the same number of algo names and algos...');
MC = numel(lcp_list);
%% Preallocate results struct
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
p_list = 2:8;
numP = numel(p_list);
for mc = 1:MC 
    %% Build Mat-Vec
    b = lcp_list(mc).b;
    F = lcp_list(mc).F;
    C = lcp_list(mc).C;
    n = numel(b);
    p2A = cell(numP,1);
    for ix_p = 1:numP
        p = p_list(ix_p);
        p2A{ix_p} = getMatVec(Fparams, F, C, p);
    end
    A = @(x) AWithMem(x, p2A{end}); % HiFi A is considered 'true' A
    %% Build cost function 
    Acnt = @(x) Acounter(x,A); % HiFi p matVec
    x0 = zeros(n,1);
    fg = @(x, Ax) quadraticLoss(x, Acnt, b, Ax);
    %% Fill the lcpOpts with problem specific information
    lcpOpts.A = Acnt;
    lcpOpts.b = b;
    mcGoodFlag = true;
    %% Obtain 'true solution' by running PGD for a long time
%     try 
        this_lcpOpts = defaultLCPOpts(lcpOpts,x0);
        this_lcpOpts.stepSize.kappa = 'uniform';
        this_lcpOpts.stepSize.eta = 'opt';
        this_lcpOpts.kkt_rel = 1e-9;
        this_lcpOpts.kkt_abs = 1e-9;
        [xstar,~] = proxQuasiNewton(fg, x0, this_lcpOpts);
        lcpOpts.errFcn = {
            @(x) abs_kkt(x, A, b); 
            @(x) rel_kkt(x, A, b);
            @(x) Acnt('cnt');
            @(x) rel_iter(x,xstar);
            @(x) abs_iter(x,xstar); 
        };
      
%     catch
%     mcGoodFlag = false;
%     end
    Acnt('reset');
    %% Run all the algorithms
    for ixAlgo = 1:numAlgo
        if ~mcGoodFlag
            break
        end
        name = algoNames{ixAlgo};
        this_lcpOpts = defaultLCPOpts(lcpOpts,x0);
        if contains(name, 'kappa') && contains(name, 'bb')
            this_lcpOpts.stepSize.kappa = 'bb1';
        elseif contains(name, 'kappa = 1')
            this_lcpOpts.stepSize.kappa = 'uniform';
        elseif contains(lower(name), 'subspace') 
            if contains(lower(name), 'pgd')
                this_lcpOpts.subSpaceMin.innerStepSelection = 'pgd';
                this_lcpOpts.stepSize.kappa = 'opt';
            else
                this_lcpOpts.subSpaceMin.innerStepSelection = 'pqn';
                this_lcpOpts.stepSize.kappa = 'uniform';
            end
        else    
             this_lcpOpts.stepSize.kappa = 'uniform';
        end
        if contains(name, 'eta = 1')
            this_lcpOpts.stepSize.eta = 'uniform';
        elseif strcmpi(name, 'eta^*')
            this_lcpOpts.stepSize.eta = 'opt';
        end
        algo = algoHndls{ixAlgo};
        tic
        [x, info] = algo(fg, x0, this_lcpOpts);
        results(ixAlgo).name = name;
        results(ixAlgo).time(mc) = toc();
        results(ixAlgo).iters(mc) = info.iter;
        results(ixAlgo).kkt(mc) = info.kkt;
        results(ixAlgo).matVecs(mc) = Acnt('reset');
        results(ixAlgo).errHist{mc} = info.errHist;
        results(ixAlgo).iterHist{mc} = info.iterHist;
        for ixMetric = 1:numel(lcpOpts.errFcn)
            hndl = lcpOpts.errFcn{ixMetric};
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
save(fullfile(root, saveFile), ...
    'results', 'mcGood', 'lcpOpts');
end

function Ax = AWithMem(x, A)
persistent mem
if isempty(mem) || isempty(x)
    mem = cell(2,0);
end
ii = size(mem,2);
for i = ii:-1:1
    xi = mem{1,i};
    if norm(xi -x)/norm(x) < 1e-12 || norm(xi -x) < 1e-12
        Ax = mem{2,i};
        return 
    end
end
Ax = A(x);
mem{1,ii+1} = x;
mem{2,ii+1} = Ax;
end % AWithMem

function Ax = Acounter(x, A)
persistent matVecCnt 

if isempty(matVecCnt)
    matVecCnt = 0;
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

Ax = A(x);
matVecCnt = matVecCnt + 1;
end

function e = abs_kkt(x, A, b)
    phi = min(x,A(x) + b);
    e = 1/2*dot(phi, phi);
end % abs_kkt

function diff = rel_kkt(x, A, b)
    persistent olde 
    if ischar(x) && strcmpi(x,'reset')
        olde = [];
        diff = NaN;
        return
    end
    phi = min(x,A(x) + b);
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