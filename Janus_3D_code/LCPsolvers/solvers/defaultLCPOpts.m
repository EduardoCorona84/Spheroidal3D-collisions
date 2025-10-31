function [opts, info] = defaultLCPOpts(opts,x0)
n = numel(x0);
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'solver')
    opts.solver = 'bbpgd';
end

if ~isfield(opts, 'n')
    opts.n = n;
end

if ~isfield(opts, 'max_iter')
    opts.max_iter = 100;
end

if ~isfield(opts, 'kkt_rel')
    opts.kkt_rel = 1e-4;
end

if ~isfield(opts, 'kkt_abs')
    opts.kkt_abs = 1e-6;
end

if ~isfield(opts, 'arg_rel')
    opts.arg_rel = [];
end

if ~isfield(opts, 'arg_abs')
    opts.arg_abs = [];
end

if ~isfield(opts, 'step_abs')
    opts.step_abs = [];
end

if ~isfield(opts, 'gamma')
    opts.gamma = 0.8;
end

if ~isfield(opts, 'tau')
    % warning("Tau should be set as the 1/Lhat the the best estimate of the lipschitz constant of the the A mat");
    opts.tau = 1;
end

if ~isfield(opts, 'tau_min')
    opts.tau_min = 1e-14;
end

if ~isfield(opts, 'tau_max')
    opts.tau_max = Inf;
end

if ~isfield(opts, 'stepSize')
    opts.stepSize = struct('init', 'uniform',...
        'kappa', [],...
        'eta', []);
    switch lower(opts.solver)
        case 'bbpgd'
            opts.stepSize.kappa = 'bb1';
            opts.stepSize.eta = 'uniform';
        case 'proxquasinewton'
            opts.stepSize.kappa = 'uniform';
            opts.stepSize.eta = 'opt';
        case 'subspacemin'
            opts.stepSize.kappa = 'uniform';
            opts.stepSize.eta = 'opt';
        case 'semismoothnewton'
            opts.stepSize.init = 'uniform';
            opts.stepSize.kappa = 'uniform';
            opts.stepSize.eta = 'uniform';
        otherwise
            opts.stepSize.kappa = 'uniform';
            opts.stepSize.eta = 'opt';
    end
    % For testing scripts, all the name of the algo 
    if isfield(opts, 'name')
        if contains(opts.name, 'kappa') && contains(opts.name, 'bb')
            opts.stepSize.kappa = 'bb1';
        elseif contains(opts.name, 'kappa = 1')
            opts.stepSize.kappa = 'uniform';
        elseif contains(opts.name, 'kappa^*')
            opts.stepSize.kappa = 'opt';
        end
        if contains(opts.name, 'eta = 1')
            opts.stepSize.eta = 'uniform';
        elseif contains(opts.name, 'eta^*')
            opts.stepSize.eta = 'opt';
        end
    end
end

if ~isfield(opts, 'r')
    opts.r = min(20, n);
else 
    assert(opts.r > 0, 'Memory/effective-rank or hessian must by positive');
    opts.r = min(opts.r, n);
end

if ~isfield(opts, 'qnUpdate') || isempty(opts.qnUpdate)
    opts.qnUpdate = 'bfgs';
else 
    assert( strcmpi('bfgs', opts.qnUpdate) ...
        || strcmpi('sr1', opts.qnUpdate), 'Must choose SR1 or BFGS')
end

  
if ~isfield(opts, 'prox')|| isempty(opts.prox)
    opts.prox = struct( ...
        'maxiter', 1000, ...
        'res_abstol', 0, ...
        'res_reltol', 0, ...
        'alp_abstol', 0, ...
        'alp_reltol', 1e-12, ...
        'verbose', false, ...
        'runCVX', false);
end

%% initialize info
info = struct('kkt', [], ...
    'iter',[],...
    'flag', [],...
    'msg',[]);
if ~isfield(opts, 'errFcn') || isempty(opts.errFcn)
    opts.errFcn = [];
elseif isa(opts.errFcn,'function_handle')
    info.errHist = zeros(opts.max_iter+1,1);
    info.errHist(1) = opts.errFcn(x0);
else
    assert(iscell(opts.errFcn) && ...
            all(cellfun(@(f) isa(f,'function_handle'), opts.errFcn)));
    info.errHist = zeros(opts.max_iter+1,numel(opts.errFcn));
    for i = 1:numel(opts.errFcn)
        fcn = opts.errFcn{i};
        info.errHist(1,i) = fcn(x0);
    end
end

if ~isfield(opts, 'storeIts') || isempty(opts.storeIts)
    opts.storeIts = false;
elseif opts.storeIts
    info.iterHist = zeros(opts.max_iter+1, n);
end

if ~isfield(opts, 'subspaceMin') 
    opts.subSpaceMin = struct( ...
        'innerStepSelection', 'pqn', ...
        'innerSolver','cvx', ...
        'orthoMethod','qr', ...
        'kThresh', 0, ...
        'useOneStepIter', false ...
    );
end
end % defaultLCPOpts
