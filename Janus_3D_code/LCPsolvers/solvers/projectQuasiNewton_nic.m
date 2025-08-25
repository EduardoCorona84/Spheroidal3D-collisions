function [x,info] = projectQuasiNewton_nic(fg, x0, opts)
% Nic Rummel 2025
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

%% -----------------------------------------------------
%  The main iterative loop
%  -----------------------------------------------------
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
    % first binding set
    bMask = x_k == 0 & grad_k > 0;
    H = get_H_BFGS(k, s_k, y_k, opts, bMask);
    p = -H(grad_k);
    % second binding set
    bMask = x_k == 0 & (grad_k > 0 | p < 0);
    p(bMask) = 0;
    
    % step size direction
    kappa = stepSize(k, p, x_km1, Ax_km1, opts);
    % projected quasi-newton step
    x_k = max(0, x_km1 + kappa * p);
    % Possibly a step length update after the projection
    q = x_k - x_km1;
    [eta,Aq] = stepSize(-1, q, x_km1, Ax_km1, opts);
    x_k = x_km1 + eta*q;
    k = k+1;
end % of while
end % projectedQuasiNewton_nic

function checkOpts(opts)
assert(strcmpi(opts.stepSize.kappa, 'opt') ...
    || strcmpi(opts.stepSize.kappa, 'uniform'),...
    ['Quasi Newtwon method should use a uniform step size of 1 because the '...
    'BB step size is baked into the hessian approximation (unconstrained) optimal size']);
end