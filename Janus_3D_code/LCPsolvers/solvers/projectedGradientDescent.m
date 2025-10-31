function [ x, info] = projectedGradientDescent(fg, x0, opts )
% May 2018 Wen Yan
% Edited 2025 Nic Rummel  
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
[f_k, grad_k] = fg(x_k, Ax_k);
k = 0;
while true
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts);
    if converged
        x = x_k;
        return 
    end
    k = k + 1;
    x_km1 = x_k; Ax_km1 = Ax_k; grad_km1 = grad_k;
    % Gradient descent direction
    q = -grad_k;
    % Select step size
    kappa = stepSize(k, q, x_km1, Ax_km1,  opts, s, y); 
    % Gradient Descent Step
    x_k = x_km1 + kappa * q;
    % Projection to positive orthant
    x_k = max(0, x_k);
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

end

function checkOpts(opts)
assert(strcmpi(opts.stepSize.kappa, 'opt') ...
    || contains(lower(opts.stepSize.kappa), 'bb'),...
    'Projected Gradient Descent should use the the BB step size or the (unconstrained) optimal size');
end