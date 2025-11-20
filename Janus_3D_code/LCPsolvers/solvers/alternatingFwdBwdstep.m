function [x_k,Ax_k,f_k,grad_k, eta, p, Ap] = alternatingFwdBwdstep(t, x_km1, Ax_km1, q, prox, fg, opts)
    % Same thing as fwdBwdstep, but for the alternating schemes. More work needs to be done to make this compatible with the wolfe/strong wolfe conditions (but this should evaluate, just not meaningful line search info).
    if t == 0
        x_k = x_km1;
        Ax_k = Ax_km1;
        [f_k, grad_k] = fg(x_k, Ax_k);
        eta = 0;
        p = q;
        Ap = zeros(numel(p),1);
        return 
    end
    % High Step
    xtilde = x_km1 + t * q;
    xhat = prox(xtilde);
    % Low step
    ztilde = xhat + opts.inner.damping*t * (Ax_km1 + opts.inner.A(xhat - x_km1) + opts.inner.b);
    zhat = prox(ztilde);
    % Optimal Step
    p = zhat - x_km1;
    [eta, Ap] = stepSize(-1, p, x_km1, Ax_km1, opts);
    x_k = x_km1 + eta*p;
    Ax_k = Ax_km1 + eta*Ap;
    [f_k, grad_k] = fg(x_k, Ax_k);
end