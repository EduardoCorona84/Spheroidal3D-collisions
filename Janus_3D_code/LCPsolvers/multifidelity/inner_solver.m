function [x_inner_1, info_inner, opts] = inner_solver(x_inner_0, grad_inner, opts)

    %This is a wrapper function for calling a desired inner solver.

    switch lower(opts.inner.solver)
        case 'outer prox qn'
            opts.inner.prox = opts.outer.prox;
            [x_inner_1, info_inner] = outer_prox_qn(x_inner_0, grad_inner, opts);

    end

end



