function [x, info] = optDiagPrecond(fg, x0, opts, debug)
if ~exist('debug','var') || isempty(debug)
    debug = false;
end
[opts, info] = defaultLCPOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
eta = 1;
x_k = x0;
% z_k = x0;
% z_km1 = NaN*ones(n,1);
x_km1 = NaN*ones(n,1);
grad_km1 = NaN*ones(n,1);
Ax_km1 = [];
Aq = [];
k = 0;

%% find the best diagonal preconditioner
A = opts.AA; 
% b = opts.b; 
cvx_begin SDP quiet
    variables tau(1) d(n)
    maximize tau
    subject to 
    diag(d) - tau*A == semidefinite(n)
    A - diag(d) == semidefinite(n)
    0 <= d
cvx_end

if debug 
    disp(['\kappa(A): ' num2str(cond(A)) ])
    disp(['\sqrt(\kappa(A)): ' num2str(sqrt(cond(A))) ])
    B = diag(sqrt(1./ d))*A*diag(sqrt(1 ./d));
    % c = diag(sqrt(1./ d))*b;
    disp(['\kappa(d^{-1/2} A d^{-1/2}' num2str(cond(B))])
end

while true
    [f_k, grad_k, Ax_k] = fg(x_k, Ax_km1, Aq, eta);
    %% Scale the gradient so that it is for the preconditioned problem
    % grad_k = sqrt(1./d) .* grad_k;
    grad_k = 1./d .* grad_k;
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts);
    if converged
        x = x_k;
        break
    end
    
    % s_k = z_k - z_km1;
    s_k = x_k - x_km1;
    y_k = grad_k - grad_km1;
    x_km1 = x_k;
    % z_km1 = z_k;
    grad_km1 = grad_k;
    Ax_km1 = Ax_k;
    %% Take a step in the unconstrained problem for z
    p = - grad_km1;
    kappa = stepSize(k, p, [], [], opts, s_k, y_k); 
    x_k = x_km1 + kappa * p;
    % z_k = z_km1 + kappa * p;
    % x_k =  z_k ./ sqrt(d);
    x_k = max(x_k, 0);  
    %% After the projection take the optimal step in x
    q = x_k - x_km1;
    [eta, Aq] = stepSize(-1, q, x_km1, Ax_km1, opts);
    x_k = x_km1 + eta * q;
    % z_k = sqrt(d) .* x_k; 
    % Increment the number of iterations
    k = k + 1;
end

end % proxQuasiNewton

function checkOpts(opts)
% assert(strcmpi(opts.stepSize.kappa, 'opt') ...
%     || strcmpi(opts.stepSize.kappa, 'uniform'),...
%     ['Quasi Newtwon method should use a uniform step size of 1 because the '...
%     'BB step size is baked into the hessian approximation (unconstrained) optimal size']);
assert(0 < opts.tau_min, "tau_min must be positive");
assert(opts.tau_min < opts.tau_max, "tau_max must be larger than tau_min");
assert(0 < opts.gamma && opts.gamma < 1, "gamma must be in (0,1)");
end