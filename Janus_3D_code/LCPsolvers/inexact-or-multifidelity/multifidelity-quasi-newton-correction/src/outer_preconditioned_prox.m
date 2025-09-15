function [x, info] = outer_preconditioned_prox(x0, grad, opts)

    %should make sure the indexing aligns with nic's convention
    info.iterHist = zeros(opts.max_iter, length(x0));
    info.matvecs = zeros(opts.max_iter, 1);
    for i = 1:opts.max_iter
        info.iterHist(i, :) = x0';
        grad_0 = grad(x0);
        info.matvecs(i) = info.matvecs(i) + 1;
        p = -opts.H(grad_0);
        x_1 = prox(x0 + p, opts);
        x0 = x_1;
    end
    x = x_1;

end