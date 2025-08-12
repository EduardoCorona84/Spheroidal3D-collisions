function [x, info] = outer_prox_qn(x0, grad, opts)
    %TO DO: Add support for tolerance based stopping criteria
    info.x_iters = zeros(length(x0), opts.inner.max_iter);
    info.matvecs = zeros(1, opts.inner.max_iter);
    for i=1:opts.inner.max_iter
        info.x_iters(:, i) = x0;
        grad_0 = grad(x0);
        info.matvecs(i) = info.matvecs(i) + 1;
        p = -opts.inner.prox.H(grad_0);
        x_1 = prox(x0 + p, opts.inner.prox);
        x0 = x_1;
    end
    x = x_1;

end