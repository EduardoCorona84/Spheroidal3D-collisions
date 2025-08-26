function results = testLCPWrappers()
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname, ~,~] = fileparts(mfilePath);
root = '/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/';
loadFile = 'all_data';
saveFile = ['results' '.mat'];
%% Load scenario data from file 
load(fullfile(root, loadFile), ...
    'Fparams',
    'lcp_list');
%% Set Hyperparameters

lcpOpts = struct( ...
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
%% Specify LCP solvers
algoNames = {
    'PGD (\kappa = \tau_{bb_1}, \eta = 1)'; 
    'PGD (\kappa = \tau_{bb_1}, \eta = \eta^*)'; 
    'zeroSR1 (\kappa = 1, \eta = \eta^*)';
    'L-BFGS-B';
    'Projected QuasiNewton (BFGS, \kappa = 1, \eta = \eta^*)';
    'Proximal QuasiNewton (BFGS, \kappa = 1, \eta = \eta^*)';
};
algoHndls = {
    @projectedGradientDescent; @projectedGradientDescent; 
    @zeroSr1_nic; 
    @L_BFGS_B; 
    @projectQuasiNewton_nic;  
    @proxQuasiNewton;%
};
numAlgo = numel(algoNames);
% idiot check
assert(numel(algoNames) == numel(algoHndls), 'Must specify the same number of algo names and algos...');
MC = numel(b_list);
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
    %
    A_list = cell(numP,1);
    mc = ixTime; % TODO CHANGE
    p = 4; % TODO CHANGE
    b = lcp_list(mc).b;
    C = lcp_list(mc).C;
    F = lcp_list(mc).F;
    lofi_A = getMatVec(Fparams, F, C, p);% TODO CHANGE
    for ix_p = 1:numP
        p = p_list(ix_p);
        A_list{ix_p} = getMatVec(Fparams, F, C, p);
    end
    %% Build cost function 
    Acnt = @(x) Acounter(x,A, false);
    x0 = zeros(n,1);
    fg = @(x, Ax, Aq, eta) quadraticLoss(x, Acnt,b, Ax, Aq, eta);
    %% Fill the lcpOpts with problem specific information
    lcpOpts.A = Acnt;
    lcpOpts.b = b;
    mcGoodFlag = true;
    %% Obtain 'true solution' by running PGD for a long time
    try 
        this_lcpOpts = lcpOpts;
        this_lcpOpts.stepSize.kappa = 'uniform';
        this_lcpOpts.stepSize.eta = 'uniform';
        [xstar,~] = projectedGradientDescent(fg, x0, this_lcpOpts);
        lcpOpts.errFcn = {
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
        this_lcpOpts = lcpOpts;
        if contains(name, 'kappa') && contains(name, 'bb') 
            this_lcpOpts.stepSize.kappa = 'bb1';
        elseif contains(name, 'kappa = 1') 
            this_lcpOpts.stepSize.kappa = 'uniform';
        elseif contains(name, 'kappa^*') 
            this_lcpOpts.stepSize.kappa = 'opt';
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