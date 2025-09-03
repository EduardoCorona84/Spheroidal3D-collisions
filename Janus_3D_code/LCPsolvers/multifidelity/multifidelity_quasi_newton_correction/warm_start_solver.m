function [x, info] = warm_start_solver(fg_low, x0, opts)
    
    addpath('../../solvers')
    switch lower(opts.type)
        case 'bbpgd'
            [x, info] = projectedGradientDescent(fg_low, x0, opts);

        case 'prox'
            [x, info] = proxQuasiNewton(fg_low, x0, opts);

            %etc, fill in for each solver.
    end
end