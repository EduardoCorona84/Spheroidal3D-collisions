function [x, info] = zeroSr1_nic(fg, x0, opts)
[opts, info] = defaultLCPOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
eta = 1;k = 0;
x_k = x0; Ax_k = zeros(n,1); s = []; y = [];
if any(x0 ~= 0)
    Ax_k = opts.A(x_k); s = x_k; y = Ax_k;
end
[f_k, grad_k] = fg(x_k, Ax_k);
while true
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts);
    if converged
        x = x_k;
        break
    end
    k = k + 1;
    x_km1 = x_k; Ax_km1 = Ax_k; f_km1 = f_k; grad_km1 = grad_k;
    [H, h0, u, sigma, opts] = updateHk(s, y, opts);
    % quasi-newton step direction
    q = - H(grad_km1);
    prox_k = @(xtilde) prox(xtilde, h0, u, sigma, opts);
    step_k = @(t, opts) fwdBwdstep(t, x_km1, Ax_km1, q, prox_k, fg, opts);
    [x_k, Ax_k, f_k, grad_k] = linesearch(x_km1, Ax_km1, f_km1, grad_km1, ...
        step_k, opts);
    s = x_k - x_km1;
    y = grad_k - grad_km1;
end

end % sr1_nic

function [H, h0, u, sigma, opts] = updateHk(s, y, opts)

if isempty(s) || isempty(y)
    H = @(g) g;
    h0 = 1;
    u = [];
    sigma = NaN; 
    return 
end

denom = norm(y,2)^2;
if denom < 100*eps()
    h0 = opts.tau;
else
    %% Constants from Stephens ProxQN paper
    gamma = 0.8;
    tau_min = 1e-14;
    tau_max = Inf;
    tau_bb2 = dot(s,y) / denom;
    tau_bb2 = clip(tau_bb2, tau_min, tau_max);
    if tau_bb2 == tau_min
        warning('Convexity of cost function is stagnating')
    end
    h0 = gamma * tau_bb2;
end
% H_k = h0 * I + sigma * u * u'
delta = s - h0 .* y;
denom = dot(delta, y);
sigma = sign(denom);
if (denom <= 1e-8 * norm(y,2)^2 * norm(s - h0 .* y,2)^2 || ...
        dot(y,s) <= 1e-8 )
    % The first check is for sufficient decrease the second is curvature 
    u = [];
    H = @(g) h0*g;
else
    u = delta / sqrt(sigma*denom);
    H = @(g) h0*g + sigma * u * dot(u, g);
end
end % updateHk

function xhat = prox(xtilde, h0, u, sigma, opts, debug)
if ~exist('debug','var') || isempty(debug)
    debug = false;
end

if isempty(u)
    xhat = max(0, xtilde);
    return 
end

% Sherman-Morrison update: H = h0 * I  + sigma * uu^T 
% -> b0 = 1 / h0 ; B = 1/h0 * I - sigma * (u ./ h0)u ./h0)' / (1 + sigma * u' ./ h0 *u) 
denom = sqrt(1 + sigma * dot(u / h0, u));
w = u / h0 / denom; 
xhat = prox_rank1(xtilde, h0, w, sigma, opts);

if debug
    n = length(xtilde);
    H = eye(n)*h0 + sigma*(u*u');
    B = eye(n)/h0 - sigma*(w*w');
    assert(norm(H*B - eye(n)) < 1e-6)
end
end % prox

function checkOpts(opts)
assert(strcmpi(opts.stepSize.fwd, 'opt') ...
    || strcmpi(opts.stepSize.fwd, 'uniform'),...
    ['Quasi Newtwon method should use a uniform step size of 1 because the '...
    'BB step size is baked into the hessian approximation (unconstrained) optimal size']);
end