function [opts, info] = defaultLCPOpts(opts,x0,resetQN)
%% Set Defaults
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end
if ~exist('x0', 'var') || isempty(x0)
    x0 = [];
end
if ~exist('resetQN', 'var') || isempty(resetQN)
    resetQN = false;
end
%% High level params
n = numel(x0);
if ~isfield(opts, 'solver')
    opts.solver = 'proxquasinewton';
end
opts.n = n;
%% Convergence Parameters
if ~isfield(opts, 'max_iter')
    opts.max_iter = 100;
end

if ~isfield(opts, 'kkt_rel')
    opts.kkt_rel = 1e-8;
end

if ~isfield(opts, 'kkt_abs')
    opts.kkt_abs = 1e-8;
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
%% WarmStart (default is off)
if ~isfield(opts, 'warmStart')
    opts.warmStart = false;
end
%% Step Size Parameters
if ~isfield(opts, 'stepSize')
    opts.stepSize = struct('fwd', [], 'bwd', []);
    switch lower(opts.solver)
        case 'bbpgd'
            opts.stepSize.fwd = 'bb1';
            opts.stepSize.bwd = 'uniform';
        case 'proxquasinewton'
            opts.stepSize.fwd = 'uniform';
            opts.stepSize.bwd = 'opt';
        case 'fista'
            opts.stepSize.fwd = 'bb1';
            opts.stepSize.bwd = 'opt';
        case 'subspacemin'
            opts.stepSize.fwd = 'uniform';
            opts.stepSize.bwd = 'opt';
        case 'semismoothnewton'
            opts.stepSize.init = 'uniform';
            opts.stepSize.fwd = 'uniform';
            opts.stepSize.bwd = 'uniform';
        otherwise
            opts.stepSize.fwd = 'uniform';
            opts.stepSize.bwd = 'opt';
    end
    % For testing scripts, all the name of the algo 
    if isfield(opts, 'name')
        if contains(opts.name, 'tau') && contains(opts.name, 'bb')
            opts.stepSize.fwd = 'bb1';
        elseif contains(opts.name, 'tau = 1')
            opts.stepSize.fwd = 'uniform';
        elseif contains(opts.name, 'tau^*')
            opts.stepSize.fwd = 'opt';
        end
        if contains(opts.name, 'eta = 1')
            opts.stepSize.bwd = 'uniform';
        elseif contains(opts.name, 'eta^*')
            opts.stepSize.bwd = 'opt';
        end
    end
end
%% Linesearch Parameters
if ~isfield(opts, 'linesearch')
    % Hyper parameters from pg 62 of N&W
    opts.linesearch = struct('budget', 1,...
        'c1', 1e-4, ... 
        'c2', 0.9, ...
        'tol',  1e-8);
end
%% QN Parameters
if ~isfield(opts, 'qn')
    opts.qn= struct('resetMem',false);
end
if ~isfield(opts.qn, 'm') || isempty(opts.qn.m) 
    if n == 0 
        opts.qn.m = Inf;
    else
        opts.qn.m = n;
    end
end
if ~isfield(opts.qn, 'update') || isempty(opts.qn.update)
    opts.qn.update = 'bfgs';
end
if ~isfield(opts.qn, 'S') || isempty(opts.qn.S) || resetQN
    if isinf(opts.qn.m)
        % Corresponds to using full memory
        m = n;
    else
        m = opts.qn.m;
    end
    opts.qn.rho = zeros(m,1);
    opts.qn.S = zeros(n, m);
    opts.qn.Y = zeros(n, m);
end
%% Proximal operator parameters
if ~isfield(opts, 'prox')|| isempty(opts.prox)
    opts.prox = struct('prox_B0',@(xtilde) max(xtilde,0), ...
        'maxiter', 1000, ...
        'res_abstol', 0, ...
        'res_reltol', 0, ...
        'alp_abstol', 0, ...
        'alp_reltol', 1e-5, ...
        'verbose', false, ...
        'runCVX', false);
end
%% Parameters specific to subspaceMin
if ~isfield(opts, 'subspaceMin') 
    opts.subspaceMin = struct( ...
        'innerSolver','cvx', ...
        'orthoMethod','qr', ...
        'm', n ...
    );
end
%% Acceleration Parameters
if ~isfield(opts, 'acceleration')
    opts.acceleration = struct('alpha_k', 1, 'method', 'adaptive', 'cond', 1e-1, 'restart', true, 'curvature', 'feasible', 'direction', 'feasible');
end
%% bifi 
if contains(lower(opts.solver), 'bifi')
    if ~isfield(opts.qn, 'U') || ~isfield(opts.qn, 'V') || resetQN
        opts.qn.U = zeros(n,n);
        opts.qn.V = zeros(n,n);
    end
    if ~isfield(opts,'low')
        opts.low.initWithLofi = true;
        opts.low.p = 6;
        opts.low.gmresTol = 1e-6;
        opts.low.solver = 'proxquasinewton';
        opts.low.kkt_rel = opts.kkt_rel;
        opts.low.kkt_abs = opts.kkt_abs;
        opts.low.max_iter = int64(floor(opts.max_iter/10));
        opts.low = defaultLCPOpts(opts.low,x0,true);
    end
end

if nargout == 1
    return 
end
%% Initialize info struct
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
end % defaultLCPOpts
