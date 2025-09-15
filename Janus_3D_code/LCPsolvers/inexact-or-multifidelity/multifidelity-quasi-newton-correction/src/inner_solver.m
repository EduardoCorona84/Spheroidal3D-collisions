function [x_inner_1, info_inner, opts] = inner_solver(x_inner_0, opts, s_k, y_k)

    %This is a wrapper function for calling a desired inner solver.
    %create corrected gradient
    grad_inner = @(x) opts.inner.solver_opts.A(x) + opts.outer.correction_opts.update_matrix*x + opts.inner.solver_opts.b;

    switch lower(opts.inner.solver_opts.solver)
        case 'outer preconditioned prox'
            %In this case we are going to use the same prox evaluator as the outer steps.
            %So we set the inner opts to be the same as the outer, but we retain the low-fidelity A matrix and the max_iters
            A_low = opts.inner.solver_opts.A; 
            max_iter_inner = opts.inner.solver_opts.max_iter;
            solver_type = opts.inner.solver_opts.solver;
            opts.inner.solver_opts = opts.outer.solver_opts;
            opts.inner.solver_opts.A = A_low;
            opts.inner.solver_opts.max_iter = max_iter_inner;
            opts.inner.solver_opts.solver = solver_type;

            %pass to the solver
            [x_inner_1, info_inner] = outer_preconditioned_prox(x_inner_0, grad_inner, opts.inner.solver_opts);
            
        case 'bbpgd'
            x_inner_1 = x_inner_0 - ((s_k'*s_k)/(s_k'*y_k))*grad_inner(x_inner_0);

            x_inner_1 = max(x_inner_1, 0);

            %placeholder
            info_inner = [];

        case 'bb2'

    end

end