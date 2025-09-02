function [x, info] = multifidelity_wrapper(x0, fg, fg_low, opts)
    %this function takes in an initial iterate, high and low fidelity functions, and an opts struct and returns a final iterate and an info struct.

    if opts.warm_start.enabled == true
        [opts.warm_start.solver_opts, info] = default_LCP_opts(x0, opts.warm_start.solver_opts);
        [x, warm_info] = warm_start_solver(fg, fg_low, opts.warm_start.solver_opts);
        x0 = x;
    end

    [x, info] = multifidelity_quasi_newton_corrector(x0, fg, fg_low, opts);