function x_ref = find_reference_solution(A, b)
    %This function takes in A and b and returns a reference solution x_ref by solving the high fidelity problem to high accuracy.
    if nargin < 5
        %compute reference solution
        fg = create_fg(A, b);
        problem_size = length(b);
        opts.warm.enabled = false;
        opts.inner.enabled = false;
        opts.outer.correction = false;
        opts.outer.storeIters = true;
        opts.outer.storeFuncs = true;
        opts.outer.solver_opts.tol_abs = 1e-16;
        opts.outer.solver_opts.tol_rel = 1e-16;
        opts.outer.solver_opts.A = @(x) A*x;
        opts.outer.solver_opts.b = b;
        opts.outer.max_iter = 100;
        x0 = zeros(problem_size, 1);
        [~, info_high] = multifidelity_wrapper(fg, fg, x0, opts);

        %find reference solution
        [~ , x_ref_ind] = min(info_high.outer.funcHist);
        x_ref = info_high.outer.iterHist(x_ref_ind, :)';
    end
end