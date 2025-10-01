function [x, info] = subspaceMin(fg, x0, opts)
[opts, info] = defaultLCPOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
AX = zeros(n,n);
X = zeros(n,n);
X(:,1) = x0;
x_km1 = NaN*ones(n,1);
grad_km1 = NaN*ones(n,1);
Ax_km1 = [];
Aq = [];
k = 1;
%% Do one step of Projected Gradient Descent to Get things going
[f_k, grad_k, Ax(:,k)] = fg(X(:,k), AX(:,k-1), Aq, eta);
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
k = k +1;
X(:,k) = X(:,k-1) + eta*q;
%% Start the subspace minimization
while true
    [f_k, grad_k, Ax(:,k)] = fg(X(:,k), AX(:,k-1), Aq, eta);
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts);
    if converged
        x = x_k;
        break
    end
    %% Form inner optimization problem 
    B = X(:,1:k)' * AX(:,1:k);
    c = X(:,1:k)' * b;
    cvx_begin 
        var z(k)
        min 0.5*dot(y,B*y + c)
        such that 
        0 <= X(:,1:k)'*z
    cvx_end
    x_k = X(:,1:k)*z;
    % Update storage
    k = k + 1;
    X(:,k) = x_k
    q = X(:,k) - X(:,k-1);
    eta = 1;
end

end % subspaceMin


function checkOpts(opts)
    return 
end