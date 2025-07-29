function [x, info] = onlineScaledGradient(fg, x0, opts, debug)
if ~exist('debug','var') || isempty(debug)
    debug = false;
end
[opts, info] = defaultLCPOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
eta = 1;
x_k = x0;
d_k = ones(n,1);
% z_k = x0;
% z_km1 = NaN*ones(n,1);
x_km1 = NaN*ones(n,1);
grad_km1 = NaN*ones(n,1);
d_km1 = NaN*ones(n,1);
r_k = NaN*ones(n,1);
Ax_km1 = [];
Aq = [];
k = 0;

%% Extract some cheater info for now
L = opts.L;
b = opts.b;
xlb = opts.AA \ -b;
flb = dot(xlb, 1/2*opts.AA*xlb +b); 
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
    ixNnz = d_k > 0;
    grad_k(ixNnz) = d_k(ixNnz) .* grad_k(ixNnz);
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts);
    if converged
        x = x_k;
        break
    end
    %% Compute deltas for secant equations or BB stepsize
    s_k = x_k - x_km1;
    y_k = grad_k - grad_km1;
    %% Now update the memory terms
    f_km1 = f_k;
    x_km1 = x_k;
    grad_km1 = grad_k;
    Ax_km1 = Ax_k;
    %% Take a step in the unconstrained problem for z
    p = - grad_km1;
    kappa = stepSize(k, p, [], [], opts, s_k, y_k); 
    x_k = x_km1 + kappa * p;
    x_k = max(x_k, 0);  
    %% After the projection take the optimal step in x
    q = x_k - x_km1;
    [eta, Aq] = stepSize(-1, q, x_km1, Ax_km1, opts);
    x_k = x_km1 + eta * q;
    %% Update diagonal preconditioner 
    grad_k = Ax_km1 + eta*Aq + b;
    % if k > 1 
    %     ss_k = d_k - d_km1;
    %     yy_k = r_k - r_km1;
    %     beta = (ss_k'*ss_k)/(ss_k'*yy_k);
    % else 
        beta = 1/(4*L^2); 
    % end
    d_km1 = d_k;
    % r_km1 = r_k;
    % r_k = (opts.AA*(x_km1 - grad_km1)+b) .* grad_km1 / (f_km1 - flb);
    r_k = -grad_k .* grad_km1 / (f_km1 - flb);
    d_k = max(0,d_km1 - beta * r_k);
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