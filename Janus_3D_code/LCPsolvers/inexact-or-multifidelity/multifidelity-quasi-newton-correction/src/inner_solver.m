function [x_inner_1, info_inner, opts] = inner_solver(x_inner_0, opts, s_k, y_k, info)

    %This is a wrapper function for calling a desired inner solver.
    %create corrected gradient
    switch(opts.outer.correction_opts.memory)
        case 'dense full'
            grad_inner = @(x) opts.inner.solver_opts.A(x) + opts.outer.correction_opts.update_matrix*x + opts.inner.solver_opts.b;
        
        case 'dense limited'

        case 'compact'
            switch lower(opts.outer.correction_opts.update)
                case 'sr1'
                    %use the compact representation
                    S = opts.outer.correction_opts.S(:, 1:opts.outer.correction_opts.curr_mem);
                    utility_matrix = opts.outer.correction_opts.utility_matrix(:, 1:opts.outer.correction_opts.curr_mem);

                    %compute the gradient using the compact representation
                    grad_inner = @(x) opts.inner.solver_opts.A(x) + utility_matrix*((utility_matrix'*S)\(utility_matrix'*x)) + opts.inner.solver_opts.b;
                case 'bfgs'
                    S = opts.outer.correction_opts.S(:, 1:opts.outer.correction_opts.curr_mem);
                    Y = opts.outer.correction_opts.Y(:, 1:opts.outer.correction_opts.curr_mem);
                    utility_matrix = opts.outer.correction_opts.utility_matrix(:, 1:opts.outer.correction_opts.curr_mem);

                    %using multisecant update
                    grad_inner = @(x) opts.inner.solver_opts.A(x) + Y*((Y'*S)\(Y'*x)) - utility_matrix*((S'*utility_matrix)\(utility_matrix'*x)) + opts.inner.solver_opts.b;

                case 'dfp'
                    S = opts.outer.correction_opts.S(:, 1:opts.outer.correction_opts.curr_mem);
                    Y = opts.outer.correction_opts.Y(:, 1:opts.outer.correction_opts.curr_mem);
                    utility_matrix = opts.outer.correction_opts.utility_matrix(:, 1:opts.outer.correction_opts.curr_mem);

                    %using multisecant update
                    grad_inner = @(x) opts.inner.solver_opts.A(x) + (Y - utility_matrix)*((Y'*S)\(Y'*x)) + Y*((Y'*S)\((Y - utility_matrix)'*x)) - Y*((Y'*S)\((Y - utility_matrix)'*S*((Y'*S)\(Y'*x)))) + opts.inner.solver_opts.b;

                otherwise
                    error('Unknown update option');
            end

    end

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