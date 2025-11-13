function [ x, info] = projectedGradientDescent(fg, x0, opts )
% May 2018 Wen Yan
% Edited 2025 Nic Rummel  
[opts, info] = defaultLCPOpts(opts, x0);
checkOpts(opts)
n = numel(x0); eta = 1; k = 0; tau = 1;
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
        return 
    end
    k = k + 1;
    x_km1 = x_k; Ax_km1 = Ax_k; f_km1 = f_k; grad_km1 = grad_k;
    % Gradient descent direction
    tau = stepSize(k,-grad_km1,x_km1,Ax_km1,opts,s,y);
    q = -tau * grad_km1;
    % Select step size
    prox_k = @(xtilde) max(xtilde,0);
    step_k = @(t, opts) fwdBwdstep(t, x_km1, Ax_km1, q, prox_k, fg, opts);
    % Linesearch 
    [x_k, Ax_k, f_k, grad_k] = linesearch(x_km1, Ax_km1, f_km1, grad_km1, ...
        step_k, opts);
    s = x_k - x_km1;
    y = grad_k - grad_km1;
end

end

function checkOpts(opts)
assert(strcmpi(opts.stepSize.fwd, 'opt') ...
    || contains(lower(opts.stepSize.fwd), 'bb'),...
    'Projected Gradient Descent should use the the BB step size or the (unconstrained) optimal size');
end