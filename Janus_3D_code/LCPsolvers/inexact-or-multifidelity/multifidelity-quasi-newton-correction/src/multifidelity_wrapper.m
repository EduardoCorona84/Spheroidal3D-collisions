function [x, info] = multifidelity_wrapper(fg, fg_low, x0, opts)
    %this function takes in an initial iterate, high and low fidelity functions, and an opts struct and returns a final iterate and an info struct.

    [opts, info] = set_default_opts(opts, x0, fg, fg_low);

    if opts.warm.enabled == true
        [x, warm_info] = warm_start_solver(fg_low, x0, opts.warm.solver_opts);
        x0 = x;
        info.warm = warm_info;
        if opts.outer.correction
            %compute Ab
            [~, ~, Ab] = fg(opts.outer.solver_opts.b, [], [], []);
            %update the low fidelity along the direction b
            [~, ~, A_low_b] = fg_low(opts.outer.solver_opts.b, [], [], []);
            quantity = Ab - A_low_b; %assuming the update matrix is initialized to zero
            opts.outer.correction_opts.update_matrix = ((quantity) * (quantity)' / (opts.outer.solver_opts.b' * quantity));
        end
    end



    [x, info] = multifidelity_quasi_newton_corrector(x0, fg, fg_low, opts, info);

end