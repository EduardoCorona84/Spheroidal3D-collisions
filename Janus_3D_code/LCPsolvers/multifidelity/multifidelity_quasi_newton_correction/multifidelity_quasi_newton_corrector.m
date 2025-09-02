function [x, info] = multifidelity_quasi_newton_corrector(x0, fg, fg_low, opts)
    %x0 is the initial iterate
    %fg is a full fidelity function/gradient evaluation
    %fg_low is a low fidelity function/gradient evaluation
    %b is a dense vector

    %set default opts and initialize variables in opts and info
    [opts, info] = set_default_opts(opts, x0, fg);
    %check_opts(opts)
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

        %used in correction of low fidelity or adaptivity
        s_k = x_k - x_km1;
        y_k = grad_k - grad_km1;
        x_km1 = x_k;
        grad_km1 = grad_k;
        Ax_km1 = Ax_k;
        f_km1 = f_k;

        %update the low fidelity operator using the high fidelity information 
        if opts.outer.correction == true
            opts = update_low(k, fg_low, s_k, y_k, opts, info);            
        end

        %we have all the logic for the outer (high fidelity solver) in here. Return an intermediate step x_k_half and update opts associated with the outer solver
        [x_k_half, opts] = outer_solver(k, x_km1, grad_km1, Ax_km1, s_k, y_k, opts);

        %now take the next step with the inner solver
        if opts.inner.enabled == true
            [x_k, info.inner{k + 1}, opts] = inner_solver(x_k_half, fg_low, correction, opts);
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

function opts = update_low(k, fg_low, s_k, y_k, opts, info)
    %this function updates the low fidelity correction. We can follow much of nic's logic here, but we don't need the functions to be as robust as nic's
    %and we are targetting the operator itself, and not the inverse

    %also need separate qn memories for these methods, make sure to set this in opts and be careful with implementing these things.
    switch lower(opts.outer.correction_opts.update)
        case 'sr1'
            opts.outer.correction_opts = get_correction_SR1(k, fg_low, s_k, y_k, opts.outer.correction_opts, info);
        otherwise
            error('Unknown correction option');
    end

end