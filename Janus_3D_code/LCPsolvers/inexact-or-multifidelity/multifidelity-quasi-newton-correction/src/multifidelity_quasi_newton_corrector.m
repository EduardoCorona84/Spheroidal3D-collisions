function [x, info] = multifidelity_quasi_newton_corrector(x0, fg, fg_low, opts, info)
    %x0 is the initial iterate
    %fg is a full fidelity function/gradient evaluation
    %fg_low is a low fidelity function/gradient evaluation
    %b is a dense vector

    n = numel(x0);
    eta = 1;
    x_k = x0;
    x_km1 = NaN*ones(n,1);
    grad_km1 = NaN*ones(n,1);
    Ax_km1 = [];
    Aq = [];
    f_km1 = [];
    k = 0;
    while true
        %compute gradient and function evaluation
        [f_k, grad_k, Ax_k] = fg(x_k, Ax_km1, Aq, eta);
        %convergence logic using high fidelity iterates (need to think about this more), we can also add adaptivity in here.
        [converged, opts, info] = check_convergence(k, f_k, f_km1, x_k, x_km1, info, opts);
        if converged
            x = x_k;
            break
        end

        %in my current implementation, the secant directions used for 
        %prox and the correction are the same, but we could change this in the future. (seems like not much gain)
        [s_k, y_k] = choose_direction(k, x_k, x_km1, Ax_k, grad_k, grad_km1, opts);

        x_km1 = x_k;
        grad_km1 = grad_k;
        Ax_km1 = Ax_k;
        f_km1 = f_k;

        %update the low fidelity operator using the high fidelity information 
        opts = update_low(k, fg_low, s_k, y_k, opts, info);

        %we have all the logic for the outer (high fidelity solver) in here. Return an intermediate step x_k_half and update opts associated with the outer solver
        [x_k_half, opts] = outer_solver(k, x_km1, grad_km1, Ax_km1, s_k, y_k, opts);

        %now take the next step with the inner solver
        if opts.inner.enabled == true && (k > 0 || opts.outer.warm_correction == true)
            [x_k, info_inner, opts] = inner_solver(x_k_half, opts, s_k, y_k);
            info.inner{k + 1} = info_inner;
        else
            %if no inner iteration, we leave the info struct empty
            x_k = x_k_half;
        end

        %take the optimal step size
        q = x_k - x_km1;
        [eta, Aq] = step_size(-1, q, x_km1, Ax_km1, opts.outer.solver_opts);
        x_k = x_km1 + eta*q;
        % Increment the number of iterations
        k = k + 1;
    end

end % proxQuasiNewton


function [s_k, y_k] = choose_direction(k, x_k, x_km1, Ax_k, grad_k, grad_km1, opts)
    %this function chooses the direction for the correction, either secant or the current matrix vector product
    if k == 0 
        if opts.outer.warm_correction == true
            s_k = x_k;
            y_k = Ax_k;
            return
        else 
            s_k = NaN;
            y_k = NaN;
            return
        end
    else
        switch lower(opts.outer.correction_opts.direction)
            case 'secant'
                s_k = x_k - x_km1;
                y_k = grad_k - grad_km1;
            case 'matvec'
                s_k = x_k;
                y_k = Ax_k;
            otherwise
                error('Unknown direction option');
        end
    end
end

function opts = update_low(k, fg_low, s_k, y_k, opts, info)

    if opts.outer.correction == false
        %do nothing
        return
    else
        if k == 0 && opts.outer.warm_correction == false
            %do nothing on the first iteration if not using warm correction
            return
        end
    end

    %in any other case, we update the low fidelity operator
    switch lower(opts.outer.correction_opts.update)
        case 'sr1'
            opts.outer.correction_opts = get_correction_SR1(k, fg_low, s_k, y_k, opts.outer.correction_opts, info);
        otherwise
            error('Unknown correction option');
    end
end
