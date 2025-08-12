function [x, info] = multi_fidelity_solver(A, Ahat, b, x0, opts)

    %Check and Set Options
    [opts, info] = set_default_opts(x0, opts);

    %Use a low fidelity solver to warm start the high fidelity solver.
    if opts.outer.warm_start.enabled 
        [x_outer_0, info] = warm_start(Ahat, b, x0, info, opts);
    else
        x_outer_0 = x0;
    end

    %Compute high and low fidelity operators at the initial point 

    outer_iter = 1;
    if strcmp(opts.outer.step_init, 'zero')
        A_true_0 = x0;
        grad_outer_0 = b;
        Ab = A(b);
        step_size = ((b'*b)/(b'*Ab));
        x_outer_1 = -step_size * grad_outer_0;
        
        A_true_1 = step_size * Ab;
        grad_outer_1 = A_true_1 + b;
        info = update_info('outer', info, outer_iter, x_outer_1, A_true_1, b, opts);
    else
        A_true_0 = A(x_outer_0);
        grad_outer_0 = A_true_0 + b;
        info = update_info('outer', info, outer_iter, x_outer_0, A_true_0, b, opts);

        %Take a step with the high fidelity model 
        x_outer_1 = x_outer_0 - opts.outer.step_init*grad_outer_0;
        x_outer_1 = max(x_outer_1, 0);

        outer_iter = outer_iter + 1;
        A_true_1 = A(x_outer_1);
        grad_outer_1 = A_true_1 + b;
        info = update_info('outer', info, outer_iter, x_outer_1, A_true_1, b, opts);
    end

    while outer_iter < opts.outer.max_iter

        %compute the required information for the secant condition/qn update
        s_1 = x_outer_1 - x_outer_0;
        y_1 = A_true_1 - A_true_0;
        
        %Update the low fidelity operator to satisfy the secant equation of the high fidelity operator
        opts = update_low(Ahat, s_1, y_1, opts);
        grad_inner = @(x) Ahat(x) + opts.outer.update_matrix*x + b;

        %We now go to the inner_solver
        [x_inner_1, info_inner] = inner_solver(x_outer_1, grad_outer_1, grad_inner, s_1, y_1, outer_iter, opts);

        %update the low fidelity info
        info.inner{outer_iter} = info_inner;

        %update variables
        outer_iter = outer_iter + 1;
        x_outer_0 = x_outer_1;
        x_outer_1 = x_inner_1;
        A_true_0 = A_true_1;
        A_true_1 = A(x_outer_1);
        grad_outer_1 = A_true_1 + b;
        info = update_info('outer', info, outer_iter, x_outer_1, A_true_1, b, opts);

    end

    x = x_outer_1;
    info.outer_iter = outer_iter;

end

function opts = update_low(Ahat, s_1, y_1, opts)

    switch lower(opts.outer.low_update)
        case 'sr1'
            %compute quantity of the needed matrix
            quantity = y_1 - (Ahat(s_1) + opts.outer.update_matrix*s_1);
            opts.outer.update_matrix = opts.outer.update_matrix + ((quantity) * (quantity)' / (s_1' * quantity));

        case 'bfgs'
            %compute quantity of the needed matrix
            quantity = Ahat(s_1) + opts.outer.update_matrix*s_1;
            opts.outer.update_matrix = opts.outer.update_matrix - (quantity*quantity')/(s_1' * quantity) + (y_1*y_1')/(y_1' * s_1);

        case 'dfp'
            %From wikipedia, maybe should double check
           
        otherwise
            error('Unknown low fidelity update method: %s', opts.low_update);
    end

end

