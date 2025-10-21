function [x, info] = multifidelity_quasi_newton(A, Ahat, b, x0, opts)
    % This is a quick and dirty implementation of a multifidelity prox quasi-Newton implementation. The idea is to use Ahat as the initial Hessian approximation and then solve the prox subproblem using CVX.
    % We do not have a fast way to solve the prox subproblem, so we will use CVX to solve it.
    %The main point is to just see how much a good initial Hessian approximation helps.
    n = numel(x0);
    x_k = x0;
    if all(x0 == 0)
        Ax_k = zeros(n,1);
    else 
        Ax_k = A*x_k;
    end
    f_k = (1/2)*x_k'*Ax_k + b'*x_k;
    grad_k = Ax_k + b;
    eta = 1;
    s = [];
    y = [];
    Ahat_old = Ahat;
    k = 0;
    while true
        [converged, info] = checkConvergence(k, f_k, x_k, ...
            grad_k, eta, info, opts);
        if converged
            x = x_k;
            break
        end
        k = k + 1;
        x_km1 = x_k;
        grad_km1 = grad_k;
        Ax_km1 = Ax_k;
        Ahat = update_Ahat(k, s, y, Ahat_old, opts);
        % quasi-newton step direction
        % step size direction
        %{
        x_k = prox(x_k, h0, U, V, opts);  
        we will replace prox with a CVX solve
        %}
        cvx begin
            variable x_cvx(n)
            minimize( (1/2)*quad_form(x_cvx - x_km1, Ahat) + dot((x_cvx - x_km1), grad_km1))
            subject to
                x_cvx >= 0
        cvx end
        x_k = x_cvx;
        % Possibly a step length update after the projection
        q = x_k - x_km1;
        [eta, Aq] = stepSize(-1, q, x_km1, Ax_km1, opts);
        x_k = x_km1 + eta*q;
        Ax_k = Ax_km1 + eta*Aq;
        % Increment the number of iterations
        f_k = (1/2)*x_k'*Ax_k + b'*x_k;
        grad_k = Ax_k + b;
        
        s = x_k - x_km1;
        y = grad_k - grad_km1;
        Ahat_old = Ahat;
    end


end

function A = update_Ahat(k, s, y_k, Ahat, opts)

switch lower(opts.qnUpdate)
    case 'bfgs'
        %store it densely for now as I am lazy
        switch lower(opts.memory_type)
            case 'dense'
                Ahat_s = Ahat*s;
                A = Ahat + (y_k*y_k') / (y_k'*s) - (Ahat_s*Ahat_s') / (s'*Ahat_s);
            otherwise
                error([opts.memory_type ' memory type not implemented'])
        end

    case 'sr1'
        
    otherwise
        error([opts.qnUpdate ' update not implement'])
end
end % updateHk