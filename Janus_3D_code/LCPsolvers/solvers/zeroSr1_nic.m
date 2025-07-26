function [x, info] = zeroSr1_nic(fcnGrad, x0, opts)
[opts, info] = defaultOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
kappa = 1;
x_k = x0;
x_km1 = NaN*ones(n,1);
grad_km1 = NaN*ones(n,1);
k = 0;
while true
    [f_k, grad_k] = fcnGrad(x_k);
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, kappa, info, opts);
    if converged
        x = x_k;
        break
    end
    
    s_k = x_k - x_km1;
    y_k = grad_k - grad_km1;
    x_km1 = x_k;
    grad_km1 = grad_k;
    [h0, u, sigma, opts] = updateHk(k, s_k, y_k, opts);
    % quasi-newton step direction
    p = - h0 * grad_km1;
    if ~isempty(u)
        p = p - u * dot(u, grad_k);
    end
    % step size direction
    kappa = stepSize(k, p, grad_k, opts); 
    x_k = prox(x_km1 + kappa * p, h0, u, sigma, opts);
    % Possibly a step length update after the projection
    q = x_k - x_km1;
    eta = min(1, stepSize(-1, q, grad_k, opts));
    x_k = x_km1 + eta*q;
    % Increment the number of iterations
    k = k + 1;
end

end % sr1_custom

function [h0, u, sigma, opts] = updateHk(k, s, y, opts)

if k == 0
    h0 = 1;
    u = [];
    sigma = NaN; 
    return 
end
assert(0 < opts.tau_min, "tau_min must be positive");
assert(opts.tau_min < opts.tau_max, "tau_max must be larger than tau_min");
assert(0 < opts.gamma && opts.gamma < 1, "gamma must be in (0,1)");

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
elseif dot(y,s) <= 1e-8 
    % curvature check
    u = [];
else
    u = delta / sqrt(sigma*denom);
end
end % updateHk

function xstar = prox(y, h0, u, sigma, opts, debug)
if ~exist('debug','var') || isempty(debug)
    debug = false;
end

if isempty(u)
    xstar = max(0, y);
    return 
end

% Sherman-Morrison update: H = h0 * I  + sigma * uu^T 
% -> b0 = 1 / h0 ; B = 1/h0 * I - sigma * (u ./ h0)u ./h0)' / (1 + sigma * u' ./ h0 *u) 
denom = sqrt(1 + sigma * dot(u / h0, u));
w = u / h0 / denom; 
xstar = prox_rank1(y, h0, w, sigma, opts);

if debug
    n = length(y);
    H = eye(n)*h0 + sigma*u*u';
    B = eye(n)/h0 - sigma*w*w';
    assert(norm(H*B - eye(n)) < 1e-6)
end
end % prox

function checkOpts(opts)
assert(strcmpi(opts.kappa.fwd, 'opt') ...
    || strcmpi(opts.kappa.fwd, 'uniform'),...
    ['Quasi Newtwon method should use a uniform step size of 1 because the '...
    'BB step size is baked into the hessian approximation (unconstrained) optimal size']);
end