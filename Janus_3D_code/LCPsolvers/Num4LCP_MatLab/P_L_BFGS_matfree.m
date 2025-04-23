function [ x, err,  iter, flag, convergence, msg] = P_L_BFGS_matfree( A, b, x0, max_iter, tol_rel, tol_abs, profile )
% Nic Rummel April 2025
if nargin<3
    x0 = zeros(N,1);
end
if nargin<4
    max_iter = floor(N/2);
end
if nargin<5
    tol_rel = 0.0001;
end
if nargin<6
    tol_abs = 10*eps; % Order of 10th of numerical precision seems okay
end
if nargin<7
    profile = true;
end
% set up this solvers options
opt = pqn_solopt();
opt.algo = 'PLB';
opt.use_tolx = true; % rel err : |x_k - x_{k-1}| / |x_k| < tol                 
opt.use_tolo = true; % abs err : |f_k -f_{k-1}| < tol
opt.use_tolg = true; % abs err : norm(g_k, inf) < tol
opt.use_kkt = true;  % abs err : dot(g_k, x_k) < tol
opt.tolx = tol_rel;                   
opt.tolo = tol_abs;
opt.tolk = tol_abs;
opt.tolg = tol_abs;
opt.verbose = false; 
% call to the solver wrapper
fgFcn = @(x) quadProg(x, A, b);
out = pqn_general(fgFcn, x0, opt);
x = out.x;
iter = out.iter; 
% Make err match that of BBPGD
phi = min(out.x, out.grad) ;
err = 1/2*dot(phi, phi);

if iter == max_iter
    flag = 8;
    msg = 'maxlimit';
    return 
end

flag = 6;
msg = 'local minima';

convergence = [];
if profile
    % TODO, technically, we would want this match the above err...
end

end

function [f,g] = quadProg(x, A, b)
Ax = A(x);
f = 1/2*dot(x, Ax) + dot(x,b);
if nargout < 2 
    return
end
g = Ax + b;
end