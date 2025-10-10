function [x, info] = proxQuasiNewton(fg, x0, opts)
[opts, info] = defaultLCPOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
x_k = x0;
if all(x0 == 0)
    Ax_k = zeros(n,1);
else 
    Ax_k = opts.A(x_k);
end
[f_k, grad_k] = fg(x_k, Ax_k);
eta = 1;
s = [];
y = [];
k = 0;
while true
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts);
    if converged
        x = x_k;
        break
    end
    k = k + 1;
    x_km1 = x_k;
    grad_km1 = grad_k;
    Ax_km1 = Ax_k;
    [H, h0, U, V] = updateHk(k, s, y, opts);
    % quasi-newton step direction
    p = -H(grad_k);
    % step size direction
    kappa = stepSize(k, p, x_km1, Ax_km1, opts); 
    x_k = x_km1 + kappa * p;
    x_k = prox(x_k, h0, U, V, opts);    
    % Possibly a step length update after the projection
    q = x_k - x_km1;
    [eta, Aq] = stepSize(-1, q, x_km1, Ax_km1, opts);
    x_k = x_km1 + eta*q;
    Ax_k = Ax_km1 + eta*Aq;
    % Increment the number of iterations
    
    [f_k, grad_k] = fg(x_k, Ax_k);
    s = x_k - x_km1;
    y = grad_k - grad_km1;
end

end % proxQuasiNewton

function [H, h0, U, V] = updateHk(k, s, y_k, opts)

switch lower(opts.qnUpdate)
    case 'bfgs'
        [H, h0, U, V] = get_H_BFGS(k, s, y_k, opts);
    case 'sr1'
        
    otherwise
        error([opts.qnUpdate ' update not implement'])
end
end % updateHk

function xstar = prox(y, h0, U, V, opts)
if isempty(U) && isempty(V)
    % Project with respect to the identity
    xstar = max(0, y);
elseif size(U,2) + size(V,2) == 1
    % The sign on sigma is counter intuitive, but remember B = B0 + UU' - VV'
    if ~isempty(U)
        sigma = -1;
        w = U; 
    else 
        sigma = 1;
        w = V;
    end
    xstar = prox_rank1(y, h0, w, sigma, opts);
    return  
else
    xstar = prox_rankr(y, h0, U, V, opts);
end
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