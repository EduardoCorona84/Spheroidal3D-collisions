function [x_k,Ax_k,f_k,grad_k, eta, p, Ap] = fwdBwdstep(t, x_km1, Ax_km1, q, prox, fg, opts)
    if t == 0
        x_k = x_km1;
        Ax_k = Ax_km1;
        [f_k, grad_k] = fg(x_k, Ax_k);
        eta = 0;
        p = q;
        Ap = zeros(numel(p),1);
        return 
    end
    xtilde = x_km1 + t * q;
    xhat = prox(xtilde);
    p = xhat - x_km1;
    [eta, Ap] = stepSize(-1, p, x_km1, Ax_km1, opts);
    x_k = x_km1 + eta*p;
    Ax_k = Ax_km1 + eta*Ap;
    [f_k, grad_k] = fg(x_k, Ax_k);
end