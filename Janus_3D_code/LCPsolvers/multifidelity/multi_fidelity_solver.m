function [x, info] = multi_fidelity_solver(A, Ahat, b, x0, opts)

    %Check and Set Options
    [opts, info] = set_default_opts(x0, opts);

    %Use a low fidelity solver to warm start the high fidelity solver.
    if opts.outer.warm_start.enabled 
        [x_outer_0, info] = warm_start(Ahat, b, x0, info, opts);
    else
        x_outer_0 = x0;
    end

    %For only the warm start, we just return the warm started solution. 
    if opts.outer.warm_start.only
        x = x_outer_0;
        return;
    end

    %Compute high and low fidelity operators at the initial point 

    outer_iter = 1;

    A_true_0 = A(x_outer_0);
    grad_outer_0 = A_true_0 + b;
    f_0 = 1/2*x_outer_0'*A_true_0 + b'*x_outer_0;
    info = update_info('outer', info, outer_iter, x_outer_0, A_true_0, b, opts);

    %Take a gradient descent step with the high fidelity model 
    x_outer_1 = x_outer_0 - opts.outer.step_init*grad_outer_0;
    if opts.outer.projected == true
        x_outer_1 = max(x_outer_1, 0);
    end

    outer_iter = outer_iter + 1;
    %Now go to the loop

    while outer_iter <= opts.outer.max_iter

        %evaluate full fidelity operator
        A_true_1 = A(x_outer_1);    
        grad_outer_1 = A_true_1 + b;
        f_1 = 1/2*x_outer_1'*A_true_1 + b'*x_outer_1;
        info = update_info('outer', info, outer_iter, x_outer_1, A_true_1, b, opts);

        %check the sufficient decrease
        if ~strcmp(opts.outer.adaptive, 'none') && outer_iter > 2 %the first step can suck, so just ignore the sufficient decrease
            opts = update_adaptive(opts, f_1, f_0);
        end
        %update the function eval after we update the inner fidelity iterate count
        f_0 = f_1;

        %compute the required information for the secant condition/qn update
        s_1 = x_outer_1 - x_outer_0;
        y_1 = A_true_1 - A_true_0;

        %correct the low fidelity operator if requested
        if opts.outer.correction == true
            [opts, info] = update_low(Ahat, s_1, y_1, info, outer_iter, opts);
        end
        
        %Update the low fidelity operator to satisfy the secant equation of the high fidelity operator. If there is no correction, the matrix is just zeros.
        grad_inner = @(x) Ahat(x) + opts.outer.update_matrix*x + b;

        %Go to the outer solve.
        %In the prox case this updates the approximate inverse Hesssian and takes a step. 
        %This approximation is used with as the preconditioner and metric of the low fidelity steps. 

        [x_inner_0, opts] = outer_solver(x_outer_1, grad_outer_1, s_1, y_1, outer_iter, opts);
        %This x_inner_0 will be stored as the first iterate info_inner

        %We now go to the inner_solver
        if opts.inner.enabled
            [x_inner_1, info_inner, opts] = inner_solver(x_inner_0, grad_inner, opts);
        else
            x_inner_1 = x_inner_0;
            info_inner = [];
        end

        %update the low fidelity info
        info.inner{outer_iter} = info_inner;

        %update variables
        x_outer_0 = x_outer_1;
        x_outer_1 = x_inner_1;
        A_true_0 = A_true_1;
        outer_iter = outer_iter + 1;

    end

    x = x_outer_1;
    info.outer.iter = outer_iter - 1;

end

function [opts, info] = update_low(Ahat, s_1, y_1, info, outer_iter, opts)

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
            %not efficient, just for testing
            rho_1 = 1 / (y_1' * s_1);
            quantity = Ahat(s_1) + opts.outer.update_matrix*s_1;
            opts.outer.update_matrix = opts.outer.update_matrix - (quantity*y_1')*rho_1 - (y_1*quantity')*rho_1 + (y_1*s_1'*quantity*y_1')*rho_1^2 + (y_1*y_1')*rho_1;

        otherwise
            error('Unknown low fidelity update method: %s', opts.low_update);
    end

    if opts.outer.store_updates
        info.outer.updates{outer_iter} = opts.outer.update_matrix;
    end

end

function opts = update_adaptive(opts, f_1, f_0)
    %check if the iteration decrease
    if f_1 - f_0 < 0
        %do nothing if the objective decreased
        return
    else
        switch lower(opts.outer.adaptive)
            case 'high'
                opts.inner.enabled = false;
            case 'halving'
                if opts.inner.max_iter == 1
                    opts.inner.enabled = false;
                else
                    opts.inner.max_iter = floor(opts.inner.max_iter / 2);
                end
            otherwise
                error('Unknown adaptive method: %s', opts.outer.adaptive);
        end
    end
end

