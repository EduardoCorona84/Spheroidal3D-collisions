% P_LBFGS -- This function solves the following optimization problem
%  min (1/2)x^TAX + x^Tb, subject to x \ge 0
% based on projected L-BFGS method in 'Tackling Box-Constrained
%   Optimization via a New Projected Quasi-Newton Approach' by Dongmin et
%   al
%
% Usage:
%       [ x err  iter flag convergence msg] = P_LBBFGS(A, b, x0, max_iter,
%           tol_rel, tol_abs, profile)
%
% A        -- matrix or mat-vec operator in the objective fun
% b        -- vector in the objective fun
% x0       -- Starting vector (useful for warm-starts) -- *can* be zero.
% max_iter -- maximum iterations
% tol_rel  -- relative error tolerance \|x_{k+1} - x_k\| / \|x_k\|
% abs_rel  -- absolute error tolerance \|x_{k+1} - x_k\|
% profile  -- (bool) whether to profile the optimizer
% errFcn   -- (fcn) optional error function
% max_mem  -- (bool) whether to profile the optimizer
% max_line_search  -- (int) used in Armijo Line Search to limit number of 
%                     function evaluations 
% sigma    -- (real) used in Armijo Line Search to determine significant
%             decreas in objective
% beta     -- (real) used in Armijo Line Search to scale step size in back
%             tracing
%
% Extended from code from Version 1.2 (c) 2009  Dongmin Kim  and Suvrit Sra
%
% April 2024 Nic Rummel
function [ x, err,  iter, flag, convergence, msg] = P_LBFGS(A, b, x0, ...
    max_iter, tol_rel, tol_abs, profile, errFcn, max_mem, ...
    max_line_search, sigma, beta)

if nargin < 2
    error('Too few arguments');
end

N    = length(b); % Number of variables
flag = 1;

%--- Make sure we got good working default values -------------------------
if ~exist('x0','var') || isempty(x0)
    x0 = zeros(N,1);
end
if ~exist('max_iter','var') || isempty(max_iter)
    max_iter = floor(N/2);
end
if ~exist('tol_rel','var') || isempty(tol_rel)
    tol_rel = 0.0001;
end
if ~exist('tol_abs','var') || isempty(tol_abs)
    tol_abs = 10*eps; % Order of 10th of numerical precision seems okay
end
if ~exist('profile','var') || isempty(profile)
    profile = true;
end
if ~exist('errFcn','var') || isempty(errFcn)
    errFcn = @(x) NaN; 
end
if ~exist('max_mem','var') || isempty(max_mem)
    max_mem = 20;
end
if ~exist('max_line_search','var') || isempty(max_line_search)
    max_line_search = 100; 
end
if ~exist('sigma','var') || isempty(sigma)
    sigma = 0.298;
end
if ~exist('beta','var') || isempty(beta)
    beta = 0.0498;
end

%--- to convert the flag return into a readable msg(flag)  ----------------
msg = {'preprocessing';  % flag = 1
    'iterating';      % flag = 2
    'relative';       % flag = 3
    'absolute';       % flag = 4
    'stagnation';     % flag = 5
    'local minima';   % flag = 6
    'nondescent';     % flag = 7
    'maxlimit'        % flag = 8
    };

%--- Make sure all values are valid ---------------------------------------
max_iter = max(max_iter,1);
tol_rel  = max(tol_rel,0);
tol_abs  = max(tol_abs,0);
x0       = max(0,x0);

%--- Setup values need while iterating ------------------------------------

convergence = [];
if profile
    convergence = zeros(max_iter,1); % Used when profiling to measure the convergence rate
end

% --- Start of algorithm --------------------------------------------------
fx_Ax = @(x,Ax) 1/2*dot(x, Ax) + dot(x,b) ;
gfx_Ax = @(Ax) Ax + b ;

%% ------------------------------------------------------
%  OTHER INITIALIZATION
%  ------------------------------------------------------
rho = ones(max_mem, 1);
alp = rho;
delx = zeros(length(x0), max_mem);
delg = delx;
last = 1;

x = x0;
Ax = A*x; 
oldx = x0;
obj = fx_Ax(x, Ax);
grad = gfx_Ax(Ax);
oldgrad = grad;
srch = -oldgrad;

[x, Ax, ~, ~] = pqn_line_search(A, x, obj, srch, oldgrad, fx_Ax, ...
    max_line_search, sigma, beta);
phi = min(grad, x);
err = 1/2 * dot(phi, phi);
grad = gfx_Ax(Ax);
obj = fx_Ax(x, Ax);
%% -----------------------------------------------------
%  The main iterative loop
%  -----------------------------------------------------
iter    = 1;           % Iteration counter
flag    = 2;           % 'iterating'
while iter < max_iter
    gp = find(x == 0 & grad > 0); % TODO: this is inefficient. No find is necessary leave as a mask

    % Compute L-BFGS.
    delx(:, last) = x - oldx;
    delg(:, last) = grad - oldgrad;
    delx(gp, :) = 0;
    delg(gp, :) = 0;

    oldx = x;
    oldgrad = grad;

    grad(gp) = 0;
    srch = grad;

    rho(last) = 1 / (delx(:, last)' * delg(:, last));
    pt = last;

    for i = 1 : min(iter, max_mem)
        alp(pt) = rho(pt) * delx(:, pt)' * srch;
        srch = srch - alp(pt) * delg(:, pt);
        pt = max_mem - mod(-pt + 1, max_mem);
    end

    srch = 1 / rho(last) / (delg(:, last)' * delg(:, last)) * srch;

    for i = 1 : min(iter, max_mem)
        pt = mod(pt, max_mem) + 1; % TODO : Notice that we are not actually throwing away memory here...
        b = rho(pt) * delg(:, pt)' * srch;
        srch = srch + (alp(pt) - b) * delx(:, pt);
    end

    last = mod(last, max_mem) + 1;

    srch = -srch;
    srch(gp) = 0;
    [x, Ax, ~, ~] = pqn_line_search(A, x, obj, srch, oldgrad, fx_Ax, ...
        max_line_search, sigma, beta);

    grad = gfx_Ax( Ax );
    obj  = fx_Ax(x, Ax);
    if profile
        convergence(iter) = errFcn(x) ; %#ok<AGROW>
    end
    % termination
    % updated to match the BBPGD in the other code.
    phi = min(grad, x);
    old_err = err; 
    err = 1/2*dot(phi,phi);
    if (abs(err - old_err) / abs(old_err)) < tol_rel  % Relative stopping criteria
        flag = 3;
        break;
    end
    if err < tol_abs   % Absolute stopping criteria
        flag = 4;
        break;
    end
    iter = iter + 1; 
end % of while

if iter >= max_iter
    flag = 8;
    iter = iter - 1;
end

msg = msg{flag};

end % P_LBFGS

%%  Armijo along projection arc.
%%
function [x, Ax, flag, step] = pqn_line_search(A, oldx, oldobj, srch, grad, ...
    fx_Ax, max_line_search, sigma, beta)

step = 1;
x = oldx + step * srch;
flag = -1;

for i = 1 : max_line_search
    x(x < 0) = 0;
    delx = oldx - x;
    Ax = A*x;
    fc = fx_Ax(x,Ax); 
    if oldobj - fc >= sigma * grad' * delx
        flag = 1;
        return;
    end

    step = step * beta;
    x = oldx + step * srch;
end

x = oldx;  
warning('Line search failed.') % we should never get to this line

end