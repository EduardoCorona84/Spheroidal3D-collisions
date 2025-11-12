function [ x, info] = minmap_newton(fg, x0, opts, debug)
if ~exist('debug','var') || isempty(debug)
    debug = false;
end
[opts, info] = defaultLCPOpts(opts, x0);
% Check that options are set properly for this solver
checkOpts(opts);
% initialize
x_k = x0;
n = numel(x_k);
if all(x0 == 0)
    f_k = 0; 
    Ax_k = zeros(n,1);
    grad_k = opts.b;
else 
    Ax_k = opts.A(x_k);
    [f_k, grad_k] = fg(x_k,Ax_k);
end
x_km1=[]; grad_km1=[]; eta = 1;
k = 0;
while true
    if debug 
        fprintf('f_k = %.4g', f_k);
    end
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts, x_km1, grad_km1);
    if converged
        x = x_k;
        return
    end
    k = k + 1;
    x_km1 = x_k;  grad_km1 = grad_k; Ax_km1 = Ax_k;
    %--- Solve the Newton system --------------------------------------------
    S       = find(grad_k < x_km1);
    % Newton Step
    maxCGIter = 100;
    tolCG = 1e-4;
    [p,~] = cgs(@(x) Japp(opts.A,S,x),  min(grad_km1, x_km1), ...
        tolCG, maxCGIter, [],[], x_km1); % Guess the last iterate for initial guess to solution
    p = -p;
    % Optimal step size
    [eta, Ap] = stepSize(k, p, x_km1, Ax_km1, opts);
    % 
    x_k = x_km1 + eta*p;
    Ax_k = Ax_km1 + eta*Ap;
    [f_k, grad_k] = fg(x_k, Ax_k);
end

end

function y = Japp(A,S,x)

y = x; 
ytmp = A(x);
y(S,:) = ytmp(S,:); 

end

function checkOpts(opts)
    
end