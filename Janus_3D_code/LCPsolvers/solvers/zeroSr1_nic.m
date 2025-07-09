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
        iter = k;
        return
    end

    s = x_k - x_km1;
    y = grad_k - grad_km1;
    [h0, u, sigma, opts] = updateHk(s, y, opts);

    x_km1 = x_k;
    grad_km1 = grad_k;

    % quasi-newton step -k * H_k * g_k
    p = - h0 * grad_km1;
    if ~isempty(u)
        p = p - u * dot(u, grad_k);
    end
    % if isfield(opts, 'Q')
    %     % For QP, this is the optimal step length (see page 56 of N&W)
    %     kappa = -dot(p, grad_k)/ dot(p, opts.Q*p);
    %     % c1 = 1e-4;
    %     % c2 = 0.9;
    %     % [f_test, grad_test] = fcnGrad(x_k + kappa * p);
    %     % assert( f_test <= f_k + c1*kappa*dot(grad_k, p), 'Sufficient decrease condition not satisfied');
    %     % assert( dot(grad_test, p) >= c2 *dot(grad_k, p), 'Curvature condition not satisfied');
    % end
    y = x_km1 + opts.kappa*p;
    x_k = prox(y, h0, u, sigma, opts);
    % if isfield(opts, 'Q') && norm(x_k - y) / norm(y) > 100*eps()
    % q = x_k - x_km1;
    % eta = min(-dot(q, grad_k)/ dot(q, opts.Q*q),1);
    % x_k = x_km1 + eta*q;
    % end
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

function [h0, u, sigma, opts] = updateHk(s, y, opts)

assert(0 < opts.tau_min, "tau_min must be positive");
assert(opts.tau_min < opts.tau_max, "tau_max must be larger than tau_min");
assert(0 < opts.gamma && opts.gamma < 1, "gamma must be in (0,1)");

N = length(s);
denom = norm(y,2)^2;
if denom < 100*eps()
    h0 = opts.tau;
else
    tau_bb2 = dot(s,y) / denom;
    tau_bb2 = clip(tau_bb2, opts.tau_min, opts.tau_max);
    if tau_bb2 == opts.tau_min
        warning('Convexity of cost function is stagnating')
    end
    
    h0 = opts.gamma * tau_bb2;
end
% H_k = h0 * I + sigma * u * u'
delta = s - h0 .* y;
denom = dot(delta, y);
sigma = sign(denom);
if denom <= 1e-8 * norm(y,2)^2 * norm(s - h0 .* y,2)^2 
    u = [];
    return 
end

u = delta / sqrt(sigma*denom);

end % updateHk

function xstar = prox(y, h0, u, sigma, opts)
% Sherman-Morrison update: H = h0 * I  + sigma * uu^T 
% -> b0 = 1 / h0 ; B = 1/h0 * I - sigma * (u ./ h0)u ./h0)' / (1 + sigma * u' ./ h0 *u) 
denom = sqrt(1 + sigma * dot(u / h0, u));
w = u / h0 / denom; 
 
xstar = prox_rank1(y, h0, w, sigma, opts);
if ~isempty(u)
    n = length(y);
    H = eye(n)*h0 + sigma*u*u';
    B = eye(n)/h0 - sigma*w*w';
    assert(norm(H*B - eye(n)) < 1e-6)
end
end % prox

function [converged, errStruct] = checkConvergence(k, f_k, x_k, ...
    grad_k, errStruct, opts)
if k >= opts.max_iter 
    converged = true; 
    return
end
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