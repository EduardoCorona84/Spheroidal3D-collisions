function [opts, info] = defaultLCPOpts(opts,x0)
n = numel(x0);
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'solver')
    opts.solver = 'proxquasinewton';
end

if ~isfield(opts, 'n')
    opts.n = n;
end
%% Convergence Parameters
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
    opts.qn= struct('m',[],'update',[],'resetMem',[],'rho',[],'S',[],'Y',[]);
end
if ~isfield(opts.qn, 'm') || isempty(opts.qn.m)
    opts.qn.m = n;
end
if ~isfield(opts.qn, 'update') || isempty(opts.qn.update)

    opts.qn.update = 'bfgs';
end
if ~isfield(opts.qn, 'S') || isempty(opts.qn.S)
    m = opts.qn.m;
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
