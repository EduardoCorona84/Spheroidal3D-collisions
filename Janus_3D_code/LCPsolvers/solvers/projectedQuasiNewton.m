function [ x, info] = P_L_BFGS( fg, x0, opts )
% Nic Rummel April 2025
if ~exist('opts','var') || isempty(opts)
    opts = defaultOpts();
end

% set up this solvers options
pqnOpt = pqn_solopt();
pqnOpt.algo = 'PLB';
pqnOpt.maxmem = opts.r;
pqnOpt.use_tolx = false; % |x_k - x_{k-1}| / |x_k| < tol                 
pqnOpt.use_tolo = false; % |f_k -f_{k-1}| < tol
pqnOpt.use_tolg = false; % norm(g_k, inf) < tol
pqnOpt.use_kkt = true;   % dot(g_k, x_k) < tol
pqnOpt.maxit = opts.max_iter;
pqnOpt.tolx = NaN;                   
pqnOpt.tolo = NaN;
pqnOpt.tol_relk = opts.tol_rel;
pqnOpt.tol_absk = opts.tol_abs;
pqnOpt.tolg = NaN;
pqnOpt.verbose = false; 
pqnOpt.errFcn = opts.errFcn;
% call to the solver wrapper
out = pqn_general(fg, x0, pqnOpt);
x = out.x;
info.iter = out.iter; 
% Make err match that of BBPGD
phi = min(out.x, out.grad) ;
info.kkt = 1/2*dot(phi, phi);

if ~isempty(opts.errFcn)
    info.errHist = out.errHist;
end

if info.iter == opts.max_iter
    info.flag = 8;
    info.msg = 'maxlimit';
    return 
end

info.flag = 6;
info.msg = 'local minima';


end
