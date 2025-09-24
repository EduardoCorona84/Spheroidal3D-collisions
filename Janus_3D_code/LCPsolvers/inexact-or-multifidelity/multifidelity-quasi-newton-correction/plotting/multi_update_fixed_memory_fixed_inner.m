function data = multi_update_fixed_memory_fixed_inner(A, A_low, b, parameters, opts, x_ref)

    %This function takes in A, A_low, b, and a struct of parameters (info below) and then returns the data struct which contains the relative and absolute error for a variety of inner iterations compared to a reference solution x_ref if provided. If x_ref is not provided, then the function will compute a reference solution using a high fidelity solve.


    %parameters is a struct that contains the following
    %1. parameters.update, the update that is being used
    %2. parameters.memory.type, the memory type
    %2. parameters.memory.r, the memory size if applicable
    %3. parameters.inner_iters, the number of inner iterations to be used

    %also, can optionally pass an opts struct that will initialize the solver options
    %if opts is not provided, then we will use default options in set_default_opts

    if nargin < 5
        x_ref = find_reference_solution(A, b);
    end

    data = struct();
    for i = 1:length(parameters.inner_iters)
        update = parameters.update{i};
        opts.inner.enabled = true;
        opts.inner.solver_opts.max_iter = parameters.inner_iters(i);
        opts.outer.correction = true;
        opts.outer.storeIters = true;
        opts.outer.correction_opts.r = parameters.memory.r;
        opts.outer.correction_opts.update = update;
        opts.outer.correction_opts.memory = parameters.memory.type;
        opts.outer.solver_opts.tol_abs = 1e-16;
        opts.outer.solver_opts.tol_rel = 1e-16;
        opts.outer.solver_opts.A = @(x) A*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) A_low*x;
        opts.inner.solver_opts.b = b;
        opts.outer.max_iter = 100;
        x0 = zeros(length(b), 1);
        [~, info] = multifidelity_wrapper(create_fg(A, b), create_fg(A_low, b), x0, opts);

        %create error history
        abs_error = vecnorm(info.outer.iterHist' - x_ref);
        rel_error = abs_error/norm(x_ref);

        %populate the data struct
        data.(['update_' update]).abs_error = abs_error;
        data.(['update_' update]).rel_error = rel_error;
    end
end