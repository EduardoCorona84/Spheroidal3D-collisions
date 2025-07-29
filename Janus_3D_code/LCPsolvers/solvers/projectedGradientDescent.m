function [ x, info] = projectedGradientDescent(fg, x0, opts )
% May 2018 Wen Yan
% Edited 2025 Nic Rummel  
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
        return 
    end
    s_k = x_k - x_km1;
    y_k = grad_k - grad_km1;
    x_km1 = x_k; 
    grad_km1 = grad_k;
    Ax_km1 = Ax_k;
    % Gradient descent direction
    p = -grad_k;
    % Select step size
    kappa = stepSize(k, p, x_km1, Ax_km1,  opts, s_k, y_k); 
    % Gradient Descent Step
    x_k = x_km1 + kappa * p;
    % Projection to positive orthant
    x_k = max(0, x_k);
    % Possibly a step length update after the projection
    q = x_k - x_km1;
    [eta, Aq] = stepSize(-1, q, x_km1, Ax_km1, opts);
    x_k = x_km1 + eta*q;
    % Increment the number of iterations
    k = k + 1;
end

end

function checkOpts(opts)
assert(strcmpi(opts.stepSize.kappa, 'opt') ...
    || contains(lower(opts.stepSize.kappa), 'bb'),...
    'Projected Gradient Descent should use the the BB step size or the (unconstrained) optimal size');
end