function [x, info] = warm_start_solver(fg_low, opts.warm_start.solver_opts);
    
    switch lower(opts.warm_start.solver_opts.type)
        case 'bbpgd'
            [x, info] = %name of bbpgd solver

        case 'prox'
            [x, info] = %name of prox solver

            %etc, fill in for each solver.
    end
end