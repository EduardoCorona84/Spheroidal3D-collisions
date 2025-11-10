function [x, info] = zeroSr1_nic(fg, x0, opts)
[opts, info] = defaultLCPOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
eta = 1;
x_k = x0;
if all(x0 == 0)
    Ax_k = zeros(n,1);
else 
    Ax_k = opts.A(x_k);
end
s = []; y = [];
k = 0;
[f_k, grad_k] = fg(x_k, Ax_k);
while true
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts);
    if converged
        x = x_k;
        break
    end
    k = k + 1;
    x_km1 = x_k; Ax_km1 = Ax_k; grad_km1 = grad_k;
    [h0, u, sigma, opts] = updateHk(s, y, opts);
    % quasi-newton step direction
    q = - h0 * grad_km1;
    if ~isempty(u)
        q = q - u * dot(u, grad_k);
    end
    % step size direction
    kappa = stepSize(k, q, x_km1, Ax_km1, opts); 
    x_k = prox(x_km1 + kappa * q, h0, u, sigma, opts);
    % Possibly a step length update after the projection
    p = x_k - x_km1;
    [eta, Ap] = stepSize(-1, p, x_km1, Ax_km1, opts);
    x_k = x_km1 + eta*p;
    Ax_k = Ax_km1 + eta*Ap;
    % Increment the number of iterations
    [f_k, grad_k] = fg(x_k, Ax_k);
    s = x_k - x_km1;
    y = grad_k - grad_km1;
end

end % sr1_custom

function [h0, u, sigma, opts] = updateHk(s, y, opts)

if isempty(s) || isempty(y)
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
assert(strcmpi(opts.stepSize.kappa, 'opt') ...
    || strcmpi(opts.stepSize.kappa, 'uniform'),...
    ['Quasi Newtwon method should use a uniform step size of 1 because the '...
    'BB step size is baked into the hessian approximation (unconstrained) optimal size']);
end