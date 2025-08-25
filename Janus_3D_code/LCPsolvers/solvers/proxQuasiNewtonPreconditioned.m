function [x, info] = proxQuasiNewtonPreconditioned(fg, x0, opts)
[opts, info] = defaultLCPOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
eta = 1;
f_k = NaN; 
f_km1 = NaN; 
x_k = x0;
x_km1 = NaN*ones(n,1);
d_k = ones(n,1);
d_km1 = d_k;
r_k = NaN; 
r_km1 = NaN; 
grad_km1 = NaN*ones(n,1);
Ax_km1 = [];
Aq = [];
k = 0;

% alpha = 0.5; % Parameter for the convex combination of BFGS and the preconditioner
% TODO: update this to match eq 13 of Udell, I don't know why there is a
% diameter constraining the set of preconditioners
% beta = 1/ (4*opts.L^2); % step size for the precondtioner
% regret = @(f_k, f_km1) (f_k - opts.fstar) / (f_km1 - opts.fstar);
% grad_regret = @(f_km1, grad_k, grad_km1) 1 / (f_km1 - opts.fstar) * grad_k *grad_km1';
%% find the best diagonal preconditioner
A = opts.AA; 
b = opts.b;
cvx_begin SDP quiet
    variables tau(1) d(n)
    maximize tau
    subject to 
    diag(d) - tau*A == semidefinite(n)
    A - diag(d)   == semidefinite(n)
    0 <= d
    0<= tau <=1
cvx_end



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
    % f_km1 = f_k;
    x_km1 = x_k;
    grad_km1 = grad_k;
    Ax_km1 = Ax_k;
    % d_km1 = d_k;
    [H, h0, U, V] = updateHk(k, s_k, y_k, opts);
    
    %% Instead of finding the best step size, find the best diagonal preconditioner such that the constrain is not violated
    Hexp = H(eye(n));
    cvx_begin quiet
        variables alpha(1) tau(1) 
        maximize tau
        subject to 
        (alpha*Hexp + (1-alpha)*diag(d))- tau*A   == semidefinite(n)
        A - (alpha*Hexp + (1-alpha)*diag(d))   == semidefinite(n)
        0<= tau <=1    
        0<= alpha <=1  
    cvx_end
   
    
    
    % quasi-newton step direction
    % p = -H(grad_k);
    p = - (alpha*Hexp + (1-alpha)*diag(d))* grad_k;
    % step size direction
    kappa = stepSize(k, p, x_km1, Ax_km1, opts); 
    x_k = prox(x_km1 + kappa * p, h0, U, V, opts);  
    % Possibly a step length update after the projection
    q = x_k - x_km1;
    [eta, Aq] = stepSize(-1, q, x_km1, Ax_km1, opts);
    x_k = x_km1 + eta * q;
    
    % Increment the number of iterations
    k = k + 1;
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