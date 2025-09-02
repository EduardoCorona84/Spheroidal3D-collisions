function [x_inner_1, info_inner, opts] = inner_solver(x_inner_0, grad_inner, opts)

    %This is a wrapper function for calling a desired inner solver.

    switch lower(opts.inner.solver)
        case 'outer preconditioned prox'
            opts.inner.solver_opts = opts.outer.solver_opts;
            [x_inner_1, info_inner] = outer_preconditioned_prox(x_inner_0, grad_inner, opts.inner.solver_opts);
            
        case 'bb1'
            

        case 'bb2'

    end

end