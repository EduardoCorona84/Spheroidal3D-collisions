function method_comparisons(A, A_low, b)
    addpath('../src/');
    %This function takes in a matrix A (or function handle) and a vector b and then compares the performance of 
    %1. Warm Start
    %2. Only high fidelity prox quasi-Newton
    %3. Alternating evaluations no correction
    %4. Alternating evaluations with correction (SR1)
    %5. Alternating evaluations with correction (BFGS)
    %6. Alternating evaluations with correction (DFP)

    %for now, we will just do high and alternating correction
    fg = create_fg(A, b);
    fg_low = create_fg(A_low, b);
    problem_size = length(b);

    %high
    opts.warm.enabled = false;
    opts.inner.enabled = false;
    opts.outer.correction = false;
    opts.outer.storeFuncs = true;
    opts.outer.solver_opts.tol_abs = 1e-16;
    opts.outer.solver_opts.tol_rel = 1e-16;
    opts.outer.solver_opts.A = @(x) A*x;
    opts.outer.solver_opts.b = b;
    opts.outer.max_iter = 100;
    x0 = zeros(problem_size, 1);
    [~, info_high] = multifidelity_wrapper(fg, fg_low, x0, opts);

    %find reference solution
    [~ , f_ref_ind] = min(info_high.outer.funcHist);
    x_ref = info_high.outer.iterHist(f_ref_ind, :)';

    %create error history for high, the construction has it so that every iteration is one matvec
    error_high = vecnorm(info_high.outer.iterHist' - x_ref)/norm(x_ref);

    %high
    opts.warm.enabled = true;
    opts.inner.enabled = false;
    opts.outer.correction = false;
    opts.outer.storeFuncs = true;
    opts.outer.solver_opts.tol_abs = 1e-16;
    opts.outer.solver_opts.tol_rel = 1e-16;
    opts.outer.solver_opts.A = @(x) A*x;
    opts.outer.solver_opts.b = b;
    opts.outer.max_iter = 100;
    x0 = zeros(problem_size, 1);
    [~, info_high_warm] = multifidelity_wrapper(fg, fg_low, x0, opts);

    error_high_warm = vecnorm(info_high_warm.outer.iterHist' - x_ref)/norm(x_ref);


    %alternating correction
    opts.warm.enabled = true;
    opts.inner.enabled = true;
    opts.outer.correction = true;
    opts.outer.storeFuncs = true;
    opts.outer.correction_opts.direction = 'matvec';
    opts.outer.correction_opts.skips = 'true';
    opts.outer.warm_correction = true;
    opts.outer.adaptive = 'retry';
    %opts.outer.solver_opts.solver = 'bbpgd';
    %opts.inner.solver_opts.solver = 'bbpgd';
    opts.outer.solver_opts.tol_abs = 1e-16;
    opts.outer.solver_opts.tol_rel = 1e-16;
    opts.outer.solver_opts.A = @(x) A*x;
    opts.outer.solver_opts.b = b;
    opts.inner.solver_opts.A = @(x) A_low*x;
    opts.inner.solver_opts.b = b;
    opts.outer.max_iter = 100;
    opts.inner.solver_opts.max_iter = 3;
    x0 = zeros(problem_size, 1);
    [~, info_corr] = multifidelity_wrapper(fg, fg_low, x0, opts);

    %create error history for corr
    error_corr = vecnorm(info_corr.outer.iterHist' - x_ref)/norm(x_ref);

    %plot the results
    semilogy(error_high, '-o', 'DisplayName', 'High Fidelity Only');
    hold on
    semilogy(error_high_warm, '-o', 'DisplayName', 'High Fidelity Warm Start');
    semilogy(error_corr, '-o', 'DisplayName', 'Multifidelity with Correction');
    xlabel('Matvecs');
    ylabel('Error');
    legend('Location', 'best');
    title('Method Comparison');
    hold off
end