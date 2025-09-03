function [x_inner_1, info_inner, opts] = inner_solver(x_inner_0, opts)

    %This is a wrapper function for calling a desired inner solver.

    switch lower(opts.inner.solver)
        case 'outer preconditioned prox'
            %In this case we are going to use the same prox evaluator as the outer steps.
            %So we set the inner opts to be the same as the outer, but we retain the low-fidelity A matrix.
            temp = opts.inner.solver_opts.A; 
            opts.inner.solver_opts = opts.outer.solver_opts;
            opts.inner.solver_opts.A = temp;

            %create corrected gradient
            grad_inner = @(x) opts.inner.solver_opts.A(x) + opts.outer.correction_opts.update_matrix*x + opts.inner.solver_opts.b;

            %pass to the solver
            [x_inner_1, info_inner] = outer_preconditioned_prox(x_inner_0, grad_inner, opts.inner.solver_opts);
            
        case 'bb1'
            

        case 'bb2'

    end

end