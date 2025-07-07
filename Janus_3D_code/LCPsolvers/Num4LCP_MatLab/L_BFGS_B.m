function [ x, err,  iter, flag, convergence, msg] = L_BFGS_B( A, b, x0, max_iter, tol_rel, tol_abs, profile )
% Nic Rummel April 2025
N = size(A,2);
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
% The error bounds are slightly different for this code than that of BBGD
% in particular
%  - 'fctr' parameter specifies to break when
%    |f_{k+1} - f_{k} / max(f_{k+1}, f_{k}, 1) < factr*eps()
%    so we choose to make this as loose as possible via tol_rel & tol_abs
%  - 'pgtol' parameter specifies to break when
%    max{|proj g_i | i = 1, ..., n} <= pgtol
%    this is effectively a bound on the inf-norm of the gradient, so we use
%    the abs_tol
fctr = max(tol_rel, tol_abs) / eps();
opts    = struct(...
    'x0', x0, ...
    'printEvery', Inf, ...
    'm', 20, ...
    'pgtol', tol_abs, ...
    'factr', fctr, ...
    'maxIts', max_iter, ...
    'maxTotalIts',max_iter ...
    );
f = @(x) 1/2 * dot(x, A*x) + dot(x,b);
g = @(x) A*x + b;
lb = zeros(N,1);
ub = Inf(N,1);
[x, ~, info] = lbfgsb( {f,g} , lb, ub, opts );

% info.err is a Kx2 where K is the number of iterations, and the first column
% is are the objective values (f_1, ..., f_K) and the second column is the 
% norm(g, inf). Thus the error esimate is:
try 
    fK = info.err(end,1);
    fKm1 = info.err(end-1,1);
    gNorm =  info.err(end,2);
    err = min(abs(fK  - fKm1) / max([fK, fKm1, 1]), gNorm);
catch 
    err = info.err(end,2);
end
iter = info.totalIterations;
convergence = [];
if iter == max_iter
    flag = 8;
    msg = 'maxlimit';
    return 
end

flag = 6;
msg = 'local minima';
if profile
    % TODO, technically, we would want this match the above err...
    convergence = min( ...
        [Inf, abs(diff(info.err(:,1)))'], ...
        info.err(:,2) ... 
   );
end

end