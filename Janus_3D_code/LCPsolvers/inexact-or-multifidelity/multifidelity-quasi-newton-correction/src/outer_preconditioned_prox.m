function [x, info] = outer_preconditioned_prox(x0, grad, Ax_km1, p_1, Q, opts)

    %should make sure the indexing aligns with nic's convention
    info.iterHist = zeros(opts.max_iter, length(x0));
    info.matvecs = zeros(opts.max_iter, 1);
    grad_0 = Ax_km1;
    for i = 1:opts.max_iter
        info.iterHist(i, :) = x0';
        switch(opts.gradient_mode)
            %gradient mode
            case 'full'
                grad_0 = grad(x0);
            case 'step'
                grad_0 = grad_0 + grad(p_1);
            case 'reuse'
                grad_0 = Ax_km1 + opts.b;
        end
        info.matvecs(i) = info.matvecs(i) + 1;
        p = -opts.H(grad_0);

        %now we see if we want to proceed with the low fidelity step

        descent = check_descent(grad_0, -p, Q, opts);

        if descent == true
            x_1 = prox(x0 + p, opts);
            x0 = x_1;
            grad_0 = grad_0 - opts.b;
            p_1 = x_1 - x0;
        else
            %reject step, we just return the current iterate
            x_1 = x0;
            break;
        end
    end
    x = x_1;

end