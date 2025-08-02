function [opts, info] = defaultOpts(opts,x0)
% initialize info
info = struct('kkt', [], ...
    'iter',[],...
    'flag', [],...
    'msg',[]);
n = numel(x0);
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

if ~exist('n', 'var') || isempty(n)
    n = 0;
end

if ~isfield(opts, 'max_iter')
    opts.max_iter = 100;
end

if ~isfield(opts, 'tol_rel')
    opts.tol_rel = 1e-4;
end

if ~isfield(opts, 'tol_abs')
    opts.tol_abs = 1e-6;
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

if ~isfield(opts, 'kappa')
    opts.kappa = struct('init', 'uniform',...
        'fwd', 'uniform', ...
        'bwd', 'uniform');
end

if ~isfield(opts, 'r')
    opts.r = 1;
else 
    assert(opts.r > 0, 'Memory/effective-rank or hessian must by positive');
    opts.r = min(opts.r, n);
end
opts.S = zeros(n, opts.r);
opts.Y = zeros(n, opts.r);

if ~isfield(opts, 'qnUpdate') || isempty(opts.qnUpdate)
    opts.qnUpdate = 'bfgs';
else 
    assert( strcmpi('bfgs', opts.qnUpdate) ...
        || strcmpi('sr1', opts.qnUpdate), 'Must choose SR1 or BFGS')
end

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
end % defaultOpts
