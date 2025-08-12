function info = update_info(type, info, outer_iters, x_outer_0, A_true_0, b, opts)
    switch lower(type)

        case 'outer'
            info.outer.x_iters(:, outer_iters) = x_outer_0;
            info.outer.matvecs(outer_iters) = info.outer.matvecs(outer_iters) + 1; 
            if opts.outer.store_evals
                info.outer.evals(:, outer_iters) = A_true_0;
            end
            if opts.outer.store_grad
                info.outer.grad_iters(:, outer_iters) = A_true_0 + b;
            end
            if opts.outer.store_f
                info.outer.f_iters(outer_iters) = (1/2)*x_outer_0'*A_true_0 + b'*x_outer_0;
            end

    end

end