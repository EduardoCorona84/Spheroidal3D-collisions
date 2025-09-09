function [x_k_half, opts] = outer_solver(k, x_km1, grad_km1, Ax_km1, s_k, y_k, opts)

    switch lower(opts.outer.solver_opts.solver)
        case 'bb1'
            x_k_half = x_km1 - ((s_k'*s_k)/(s_k'*y_k))*grad_km1;
        case 'bb2'
            x_k_half = x_km1 - ((s_k'*s_k)/(s_k'*y_k))*grad_km1;
        case 'prox'
            %update the inverse Hessian, we have a subsection of opts dedicated to the outer solver, that we pass to all these functions so its plug and play
            opts.outer.solver_opts = update_Hk(k, s_k, y_k, opts.outer.solver_opts);
            %take the prox step.
            p = -opts.outer.solver_opts.H(grad_km1);
            % step size direction
            kappa = step_size(1, p, x_km1, Ax_km1, opts.outer.solver_opts);

            x_k_half = x_km1 + kappa * p;
            %prox step, but once again storing the parameters in opts
            x_k_half = prox(x_k_half, opts.outer.solver_opts);

    end

end

%nic's code, but now we store these things in opts
function opts = update_Hk(k, s, y_k, opts)

    switch lower(opts.qnUpdate)
        case 'bfgs'
            opts = get_H_BFGS_multi(k, s, y_k, opts);
        case 'sr1'
            
        otherwise
            error([opts.qnUpdate ' update not implemented'])
    end
end % updateHk