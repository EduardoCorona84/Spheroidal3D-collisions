function [x, info] = bindingProxQuasiNewton(fg, x0, opts)
[opts, info] = defaultLCPOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
eta = 1;
x_k = x0;
x_km1 = NaN*ones(n,1);
grad_km1 = NaN*ones(n,1);
Ax_km1 = [];
Aq = [];
k = 0;
while true
    [f_k, grad_k, Ax_k] = fg(x_k, Ax_km1, Aq, eta);
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts);
    if converged
        x = x_k;
        break
    end
    s_k = x_k - x_km1;
    y_k = grad_k - grad_km1;
    x_km1 = x_k;
    grad_km1 = grad_k;
    Ax_km1 = Ax_k;
    bMask1 = (x_k == 0 & grad_k > 0);
    % first binding set
    [H, h0, U, V] = get_H_BFGS(k, s_k, y_k, opts, bMask1);
    % quasi-newton step direction
    p = -H(grad_k);
    % second binding set
    bMask2 = (x_k == 0 & p < 0);
    p(bMask1 | bMask2) = 0;
    % step size direction
    kappa = stepSize(k, p, x_km1, Ax_km1, opts); 
    x_k = prox(x_km1 + kappa * p, h0, U, V, opts, bMask1 | bMask2);    
    % Possibly a step length update after the projection
    q = x_k - x_km1;
    [eta, Aq] = stepSize(-1, q, x_km1, Ax_km1, opts);
    x_k = x_km1 + eta*q;
    % Increment the number of iterations
    k = k + 1;
end

end % proxQuasiNewton


function xstar = prox(y, h0, U, V, opts, bMask )

if isempty(U) && isempty(V)
    % Project with respect to the identity
    xstar = max(0, y);
    return 
elseif ~any(y < 0)
    xstar = y; 
    return 
end
% Make the proximal problem smaller by marginializing over the binding set
n = numel(y);
xstar = zeros(n,1);
y = y(~bMask);
U = U(~bMask,:);
V = V(~bMask,:);
if size(U,2) + size(V,2) == 1
    % Use special rank_1 case when it is possible
    % The sign on sigma is counter intuitive, but remember B = B0 + UU' - VV'
    if ~isempty(U)
        sigma = -1;
        w = U; 
    else 
        sigma = 1;
        w = V;
    end
    xstar(~bMask) = prox_rank1(y, h0, w, sigma, opts);
    return  
end
xstar(~bMask) = prox_rankr(y, h0, U, V, opts);

end % prox

function checkOpts(opts)
assert(strcmpi(opts.stepSize.kappa, 'opt') ...
    || strcmpi(opts.stepSize.kappa, 'uniform'),...
    ['Quasi Newtwon method should use a uniform step size of 1 because the '...
    'BB step size is baked into the hessian approximation (unconstrained) optimal size']);
assert(0 < opts.tau_min, "tau_min must be positive");
assert(opts.tau_min < opts.tau_max, "tau_max must be larger than tau_min");
assert(0 < opts.gamma && opts.gamma < 1, "gamma must be in (0,1)");
end