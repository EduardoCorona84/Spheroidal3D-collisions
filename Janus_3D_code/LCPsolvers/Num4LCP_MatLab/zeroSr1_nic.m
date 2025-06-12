function [x, iter, errStruct] = zeroSr1_nic(fcnGrad, x0, opts)

opts = defaultOpts(opts, length(x0));

x_k = x0;
[f_k, grad_k] = fcnGrad(x0);
x_km1 = x_k;
grad_km1 = grad_k;
errStruct = struct( ...
        'err',zeros(opts.max_iter+1,1), ...
        'f',zeros(opts.max_iter+1,1), ...
        'xk',zeros(length(x_k), opts.max_iter+1), ...
        'reason','');
errStruct.err(1) = NaN;
errStruct.f(1) = f_k;
errStruct.xk(:,1) = x_k;
% prox-grad descent for initial step
p0 = - opts.kappa * opts.tau * grad_km1;
x_k = max(0, x_km1 + p0);
[f_k, grad_k] = fcnGrad(x_k);
% quasi-newton iterations
for k = 1:opts.max_iter
    [converged, errStruct] = checkConvergence(k, f_k, x_k, ...
        grad_k, errStruct, opts);
    if converged
        x = x_k;
        iter = k-1;
        return
    end

    s_k = x_k - x_km1;
    y_k = grad_k - grad_km1;
    [d_k, u_k, opts] = updateHk(k, s_k, y_k, opts);

    x_km1 = x_k;
    grad_km1 = grad_k;

    % quasi-newton step -k * H_k * g_k
    p = - opts.kappa * (d_k .* grad_km1 + u_k * dot(u_k, grad_k));

    xbar = x_km1 + p;
    x_k = prox_rank1(xbar, d_k, u_k);
    [f_k, grad_k] = fcnGrad(x_k);
end

warning('Maximum iterations reached')
x = x_k;
iter = opts.max_iter;

end % sr1_custom

function opts = defaultOpts(opts, N)

if ~isfield(opts, 'max_iter')
    opts.max_iter = 1000;
end

if ~isfield(opts, 'tol_rel')
    opts.tol_rel = 1e-4;
end

if ~isfield(opts, 'tol_abs')
    opts.tol_abs = 1e-11;
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
    % warning("Tau should be set as the 1/Lhat the the best estimate of the lipschitz constant of the the A mat");
    opts.kappa = 1;
end

end % defaultOpts

function [d, u, opts] = updateHk(k, s_k, y_k, opts)

assert(0 < opts.tau_min, "tau_min must be positive");
assert(opts.tau_min < opts.tau_max, "tau_max must be larger than tau_min");
assert(0 < opts.gamma && opts.gamma < 1, "gamma must be in (0,1)");

N = length(s_k);

tau_bb2 = dot(s_k,y_k) / norm(y_k,2)^2;
tau_bb2 = clip(tau_bb2, opts.tau_min, opts.tau_max);
if tau_bb2 == opts.tau_min
    warning('Convexity of cost function is stagnating')
end
d = opts.gamma* tau_bb2 * ones(N,1);

% TODO: this is a curvature check
if dot(s_k - d .* y_k, y_k) <= 1e-8 * norm(y_k,2)^2 * norm(s_k - d .* y_k,2)^2
    u = zeros(N,1);
else
    u = (s_k - d .* y_k) / sqrt(dot(s_k - d .* y_k, y_k));
end

end % updateHk

function x = prox_rank1(xbar, d, u)
if isempty(u)
     x = max(xbar, 0);
     return 
end
if all(xbar >= 0)
    x = xbar;
    return 
end
% Sherman-Morrison update: B = 1./d - vv^T
b = 1 ./ d;
denom = sqrt(1 + dot(u .* b, u));
v = b .* u / denom;
% Double checked this with the thm
alphas = - xbar ./ (d .* v);
alphas = sort(alphas);
N = length(xbar);
Ls = zeros(N,1);
for n = 1:N
    alpha = alphas(n);
    lambda = xbar + alpha * d .* v;
    mask = lambda > 0;
    Ls(n) = alpha + dot(v, xbar) - dot(v(mask), lambda(mask));
end
[ix, ~] = find(Ls < 0, 1, 'last');
if isempty(ix) % This means that all the constraints are active.
    alphastar = dot(v,xbar);
else
    alpha = alphas(ix);
    lambda = xbar + alpha * d .* v;
    mask = lambda > 0 | (lambda == 0 & d .* v > 0);
    alphastar = (dot(v,xbar) - dot(v(mask), xbar(mask))) / (dot(v(mask), d(mask).*v(mask)) - 1 );
end

% B = diag(1./d)-v*v';
% cvx_begin
%     variable z(N)
%     minimize( dot(xbar-z,B*(xbar-z)) )
%     subject to 
%         z >= 0
% cvx_end
% 
% 
% cvx_begin
%     variable a(1)
%     minimize( 1/2*sum_square(a) +  dot(a, v'*xbar))
%     subject to 
%         0 <= xbar ./ d + v *a
% cvx_end 
% 
% 
% % idiot checks
% lambdastar = xbar + alphastar * d .* v;
% assert(abs(alphastar + dot(v, xbar) - dot(v(mask), lambdastar(mask))) < 1e-12 )
% assert(alphas(ix) < alphastar && (ix == length(d) || alphastar < alphas(ix+1)) ) 

x = max(xbar + alphastar*d.*v, 0);
% x2 = max(xbar + a*d.*v, 0);
% 
% f = @(w) dot(xbar-w, (xbar-w)./d) - dot(v,xbar-w)^2;
% disp(f(x))
% disp(f(z))
% disp(f(x2))

end % prox_rank1


function [converged, errStruct] = checkConvergence(k, f_k, x_k, ...
    grad_k, errStruct, opts)
% kkt conditions / LCP being satisified is equivalen to
% the grad_k(i) = 0 perp  x_k(i) = 0
phi = min(grad_k,x_k);
old_err = errStruct.err(k);
err = 0.5*(phi'*phi);
errStruct.err(k+1) = err;
errStruct.f(k+1) = f_k;
errStruct.xk(:,k+1) = x_k;
converged = false;
if k == 1
    return 
end

if (abs(err - old_err) / abs(old_err)) < opts.tol_rel  % Relative stopping criteria
    errStruct.reason = "relative";
     converged = true;
elseif err < opts.tol_abs   % Absolute stopping criteria
    errStruct.reason = "absolute";
    converged = true;
end
if converged 
    errStruct.f = errStruct.f(1:k+1);
    errStruct.err = errStruct.err(1:k+1);
    errStruct.xk = errStruct.xk(:,1:k+1);
end

end % checkConvergence