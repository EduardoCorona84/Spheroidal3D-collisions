function [x_inner_1, info_inner, opts] = inner_solver(x_inner_0, opts, s_k, y_k, Ax_km1, p_k, info)

    %This is a wrapper function for calling a desired inner solver.
    %We have a low fidelity input of x_k1/2 = x_k + p_k.
    %If the algorithm is convering we wexpect p_k to be small and the magnitude of x_k1/2 to be dominated by x_k.
    %Thus we want we to compute a low fidelity matrix multiply of A*x_k1/2, we can do this A*x_k + A_hat*p_k, as we already have A*x_k = Ax_km1.
    %We can then use this to compute a gradient step with the low fidelity matrix.
    %Thus we need to switch if we want to use a full inner gradient evaluation, or an inner gradient evaluation based on the step only. 
    %In intuition, the step only should basically always better.
    %For now, we will deal with this in the inner solver.
    %create corrected gradient
    if opts.outer.correction == true
        switch(opts.outer.correction_opts.memory)
            case 'dense full'
                grad_inner = @(x) opts.inner.solver_opts.A(x) + opts.outer.correction_opts.update_matrix*x + opts.inner.solver_opts.b;
                Q = [];
            
            case 'dense limited'

            case 'compact'
                switch lower(opts.outer.correction_opts.update)
                    case 'sr1'
                        %use the compact representation
                        S = opts.outer.correction_opts.S(:, 1:opts.outer.correction_opts.curr_mem);
                        [Q, R, idx] = qr(S, 0);
                        tol = 1e-2;
                        Q_indices = abs(diag(R)) >= tol * abs(R(1,1));
                        final_indices = idx(Q_indices);
                        Q = Q(:, Q_indices);

                        S = S(:, final_indices);

                        utility_matrix = opts.outer.correction_opts.utility_matrix(:, 1:opts.outer.correction_opts.curr_mem);

                        utility_matrix = utility_matrix(:, final_indices);

                        %compute the gradient using the compact representation
                        grad_inner = @(x) opts.inner.solver_opts.A(x) + utility_matrix*((utility_matrix'*S)\(utility_matrix'*x)) + opts.inner.solver_opts.b;
                    case 'bfgs'
                        S = opts.outer.correction_opts.S(:, 1:opts.outer.correction_opts.curr_mem);
                        [Q, R, idx] = qr(S, 0);
                        tol = 1e-2;
                        Q_indices = abs(diag(R)) >= tol * abs(R(1,1));
                        final_indices = idx(Q_indices);
                        Q = Q(:, Q_indices);
                        S = S(:, final_indices);
                        Y = opts.outer.correction_opts.Y(:, 1:opts.outer.correction_opts.curr_mem);
                        Y = Y(:, final_indices);
                        utility_matrix = opts.outer.correction_opts.utility_matrix(:, 1:opts.outer.correction_opts.curr_mem);
                        utility_matrix = utility_matrix(:, final_indices);

                        %using multisecant update
                        grad_inner = @(x) opts.inner.solver_opts.A(x) + Y*((Y'*S)\(Y'*x)) - utility_matrix*((S'*utility_matrix)\(utility_matrix'*x)) + opts.inner.solver_opts.b;

                    case 'dfp'
                        S = opts.outer.correction_opts.S(:, 1:opts.outer.correction_opts.curr_mem);
                        %Perform Column Selection to ensure numerical stability
                        [Q, R, idx] = qr(S, 0);
                        %determine columns to keep uisng a tolerance and the diagonal of R
                        tol = 1e-2;
                        Q_indices = abs(diag(R)) >= tol * abs(R(1,1));
                        final_indices = idx(Q_indices);
                        Q = Q(:, Q_indices);
                        S = S(:, final_indices);
                        Y = opts.outer.correction_opts.Y(:, 1:opts.outer.correction_opts.curr_mem);
                        Y = Y(:, final_indices);
                        utility_matrix = opts.outer.correction_opts.utility_matrix(:, 1:opts.outer.correction_opts.curr_mem);
                        utility_matrix = utility_matrix(:, final_indices);

                        %using multisecant update
                        grad_inner = @(x) opts.inner.solver_opts.A(x) + (Y - utility_matrix)*((Y'*S)\(Y'*x)) + Y*((Y'*S)\((Y - utility_matrix)'*x)) - Y*((Y'*S)\((Y - utility_matrix)'*S*((Y'*S)\(Y'*x)))) + opts.inner.solver_opts.b;

                    otherwise
                        error('Unknown update option');
                end
        end
    else
        %no correction, just use the low fidelity matrix
        Q = [];
        grad_inner = @(x) opts.inner.solver_opts.A(x) + opts.inner.solver_opts.b;
    end

    switch lower(opts.inner.solver_opts.solver)
        case 'outer preconditioned prox'
            %In this case we are going to use the same prox evaluator as the outer steps.
            %So we set the inner opts to be the same as the outer, but we retain the low-fidelity A matrix and the max_iters
            gradient_mode = opts.inner.solver_opts.gradient_mode;
            A_low = opts.inner.solver_opts.A; 
            max_iter_inner = opts.inner.solver_opts.max_iter;
            solver_type = opts.inner.solver_opts.solver;
            adaptive_opts = opts.inner.solver_opts.adaptive;
            opts.inner.solver_opts = opts.outer.solver_opts;
            opts.inner.solver_opts.A = A_low;
            opts.inner.solver_opts.max_iter = max_iter_inner;
            opts.inner.solver_opts.solver = solver_type;
            opts.inner.solver_opts.gradient_mode = gradient_mode;
            opts.inner.solver_opts.adaptive = adaptive_opts;

            %pass to the solver
            %gradient mode as step or reuse only works with 1 inner iteration right now
            [x_inner_1, info_inner] = outer_preconditioned_prox(x_inner_0, grad_inner, Ax_km1, p_k, Q, opts.inner.solver_opts);
            
        case 'bbpgd'
            x_inner_1 = x_inner_0 - ((s_k'*s_k)/(s_k'*y_k))*grad_inner(x_inner_0);

            x_inner_1 = max(x_inner_1, 0);

            %placeholder
            info_inner = [];

        case 'bb2'

    end

end