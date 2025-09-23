function [x, info] = multifidelity_quasi_newton_corrector(x0, fg, fg_low, opts, info)
    %x0 is the initial iterate
    %fg is a full fidelity function/gradient evaluation
    %fg_low is a low fidelity function/gradient evaluation
    %b is a dense vector

    %add b correction if we warm started with a nonzero x0

    n = numel(x0);
    eta = 1;
    x_k = x0;
    x_km1 = NaN*ones(n,1);
    grad_km1 = NaN*ones(n,1);
    Ax_km1 = [];
    q = [];
    Aq = [];
    f_km1 = [];

    %for the pre iteration updates, we use k = -1 to indicate this.
    [s_k, y_k] = choose_direction(-1, x_k, x_km1, [], [], grad_km1, opts);
    opts = update_low(-1, fg_low, s_k, y_k, opts, info);

    k = 0;

    
    while true
        %compute gradient and function evaluation
        [f_k, grad_k, Ax_k] = fg(x_k, Ax_km1, Aq, eta);

        %add accept/reject logic here
        %convergence logic using high fidelity iterates (need to think about this more), we can also add adaptivity in here and accept/reject
        [converged, descent, opts, info] = check_convergence(k, f_k, f_km1, x_k, x_km1, grad_k, grad_km1, q, eta, info, opts);
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
        if opts.inner.enabled == true && (k > 0 || opts.outer.warm_correction == true) && descent == true
            [x_k, info_inner, opts] = inner_solver(x_k_half, opts, s_k, y_k, info);
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

    %Lets go through the cases:
    %If x0 = 0 (no warm start), then we cannot add any corrections the first term in any case as we have A*0 = 0; (even in the matvec updating case)
    %If x0 = a \neq 0 or warm start, then we can correct immediately using the matrix evaluation A*a = A*a, with s_k = a, y_k = A*a;
    %Similarly in the case that x0 = a \neq 0/warm, we may want to automatically add a correction of b as an independent matvec.


    %if we want to immediately add a b correction to our matrix
    if k == -1 
        if opts.outer.b_correction == true
            s_k = opts.outer.solver_opts.b;
            y_k = opts.outer.solver_opts.A(s_k);
            return
        else
            %do nothing
            s_k = [];
            y_k = [];
            return
        end
    %if we want to add a correction for the starting point assuming we have a nonzero x0
    elseif k == 0 
        if opts.outer.warm_correction == true
            s_k = x_k;
            y_k = Ax_k;
            return
        else 
            s_k = x_k - x_km1;
            y_k = grad_k - grad_km1;
            return
        end
    %all cases for k > 1 and up
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
        if k == -1 && opts.outer.b_correction == false
            %do nothing before the first iteration if we are not correcting along b.
            return
        elseif k == 0 && opts.outer.warm_correction == false
            %do nothing on the first itertion if we are not correcting for a warm start (if we have 0 x0, cannot correct)
            return
        end
    end

    %in any other case, we update the low fidelity operator
    switch lower(opts.outer.correction_opts.update)
        case 'sr1'
            opts.outer.correction_opts = get_correction_SR1(k, fg_low, s_k, y_k, opts.outer.correction_opts, info);
        case 'bfgs'
            opts.outer.correction_opts = get_correction_BFGS(k, fg_low, s_k, y_k, opts.outer.correction_opts, info);
        case 'dfp'
            opts.outer.correction_opts = get_correction_DFP(k, fg_low, s_k, y_k, opts.outer.correction_opts, info);
        otherwise
            error('Unknown correction option');
    end
end
