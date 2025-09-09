function [x, info] = multifidelity_wrapper(fg, fg_low, x0, opts)
    %this function takes in an initial iterate, high and low fidelity functions, and an opts struct and returns a final iterate and an info struct.

    [opts, info] = set_default_opts(opts, x0, fg, fg_low);

    if opts.warm.enabled == true
        [x, warm_info] = warm_start_solver(fg_low, x0, opts.warm.solver_opts);
        x0 = x;
        info.warm = warm_info;
    end

    [x, info] = multifidelity_quasi_newton_corrector(x0, fg, fg_low, opts, info);

end