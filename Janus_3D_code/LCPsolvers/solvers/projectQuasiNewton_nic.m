function [x,info] = projectQuasiNewton_nic(fg, x0, opts)
% Nic Rummel 2025
[opts, info] = defaultOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
kappa = 1;
x_k = x0;
x_km1 = NaN*ones(n,1);
grad_km1 = NaN*ones(n,1);
k = 0;
%% ------------------------------------------------------
%  OTHER INITIALIZATION
%  ------------------------------------------------------
rho = ones(opts.r, 1);
alp = rho;
delx = zeros(length(x0), opts.r);
delg = delx;
last = 1;
%% -----------------------------------------------------
%  The main iterative loop
%  -----------------------------------------------------
while true
    [f_k, grad_k] = fg(x_k);
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, kappa, info, opts);
    if converged
        x = x_k;
        break
    end
    if k <= 0
        p = -grad_k;
    else
        gp = find(x_k == 0 & grad_k > 0);
    
        % Compute L-BFGS.
        delx(:, last) = x_k - x_km1;
        delg(:, last) = grad_k - grad_km1;
        delx(gp, :) = 0;
        delg(gp, :) = 0;
    
        grad_k(gp) = 0;
        p = grad_k;
        rho(last) = 1 / (delx(:, last)' * delg(:, last));
        pt = last;
    
        for i = 1 : min(k, opts.r)
            alp(pt) = rho(pt) * delx(:, pt)' * p;
            p = p - alp(pt) * delg(:, pt);
            pt = opts.r - mod(-pt + 1, opts.r);
        end
    
        p = 1 / rho(last) / (delg(:, last)' * delg(:, last)) * p;
        
        for i = 1 : min(k, opts.r)
            pt = mod(pt, opts.r) + 1;
            b = rho(pt) * delg(:, pt)' * p;
            p = p + (alp(pt) - b) * delx(:, pt);
        end
    
        last = mod(last, opts.r) + 1;
    
        p = -p;
        p(gp) = 0;
    end
    x_km1 = x_k;
    grad_km1 = grad_k;
    % step size direction
    kappa = stepSize(k, p, grad_k, opts);
    % projected quasi-newton step
    x_k = max(0, x_km1 + kappa * p);
    % Possibly a step length update after the projection
    q = x_k - x_km1;
    eta = min(1, stepSize(-1, q, grad_k, opts));
    x_k = x_km1 + eta*q;
    k = k+1;
end % of while
end % projectedQuasiNewton_nic

function checkOpts(opts)
assert(strcmpi(opts.kappa.fwd, 'opt') ...
    || strcmpi(opts.kappa.fwd, 'uniform'),...
    ['Quasi Newtwon method should use a uniform step size of 1 because the '...
    'BB step size is baked into the hessian approximation (unconstrained) optimal size']);
end