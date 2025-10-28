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

    %high warm
    %{
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
    [~, info_high_warm] = multifidelity_wrapper(fg, fg_low, x0, opts);

    error_high_warm = vecnorm(info_high_warm.outer.iterHist' - x_ref)/norm(x_ref);
    %}


    %{
    opts.warm.enabled = false;
    opts.inner.enabled = false;
    opts.outer.correction = false;
    opts.outer.storeFuncs = true;
    opts.outer.correction_opts.direction = 'secant';
    opts.outer.correction_opts.r = 10;
    opts.outer.correction_opts.update = 'dfp';
    opts.outer.correction_opts.skips = 'true';
    opts.outer.correction_opts.memory = 'compact';
    opts.outer.warm_correction = true;
    opts.outer.adaptive = 'none';
    %opts.outer.solver_opts.solver = 'bbpgd';
    %opts.inner.solver_opts.solver = 'bbpgd';
    opts.outer.solver_opts.tol_abs = 1e-16;
    opts.outer.solver_opts.tol_rel = 1e-16;
    opts.outer.solver_opts.A = @(x) A*x;
    opts.outer.solver_opts.b = b;
    opts.inner.solver_opts.gradient_mode = 'reuse';
    opts.inner.solver_opts.A = @(x) A_low*x;
    opts.inner.solver_opts.b = b;
    opts.outer.max_iter = 100;
    opts.inner.solver_opts.max_iter = 1;
    x0 = zeros(problem_size, 1);
    [~, info_reuse] = multifidelity_wrapper(fg, fg_low, x0, opts);

    %create error history for corr
    error_reuse = vecnorm(info_reuse.outer.iterHist' - x_ref)/norm(x_ref);
    %}

    %alternating with gradient_step, no correction
    %{
    opts.warm.enabled = false;
    opts.inner.enabled = true;
    opts.outer.correction = false;
    opts.outer.storeFuncs = true;
    opts.outer.correction_opts.direction = 'secant';
    opts.outer.correction_opts.r = 10;
    opts.outer.correction_opts.update = 'dfp';
    opts.outer.correction_opts.skips = 'true';
    opts.outer.correction_opts.memory = 'compact';
    opts.outer.warm_correction = false;
    opts.outer.adaptive = 'none';
    %opts.outer.solver_opts.solver = 'bbpgd';
    %opts.inner.solver_opts.solver = 'bbpgd';
    opts.outer.solver_opts.tol_abs = 1e-16;
    opts.outer.solver_opts.tol_rel = 1e-16;
    opts.outer.solver_opts.A = @(x) A*x;
    opts.outer.solver_opts.b = b;
    opts.inner.solver_opts.gradient_mode = 'step';
    opts.inner.solver_opts.A = @(x) A_low*x;
    opts.inner.solver_opts.b = b;
    opts.outer.max_iter = 100;
    opts.inner.solver_opts.max_iter = 1;
    x0 = zeros(problem_size, 1);
    [~, info_step] = multifidelity_wrapper(fg, fg_low, x0, opts);
    

    %create error history for corr
    error_step = vecnorm(info_step.outer.iterHist' - x_ref)/norm(x_ref);
    %}

    %alternating correction dfp
    opts.warm.enabled = false;
    opts.inner.enabled = true;
    opts.outer.correction = true;
    opts.outer.storeFuncs = true;
    opts.inner.solver_opts.adaptive.enabled = false;
    opts.inner.solver_opts.adaptive.descent_parameter = 1/500;
    opts.outer.correction_opts.direction = 'secant';
    opts.outer.correction_opts.r = 20;
    opts.outer.correction_opts.update = 'dfp';
    opts.outer.correction_opts.skips = 'true';
    opts.outer.correction_opts.memory = 'compact';
    opts.outer.warm_correction = false;
    opts.outer.adaptive = 'none';
    %opts.outer.solver_opts.solver = 'bbpgd';
    %opts.inner.solver_opts.solver = 'bbpgd';
    opts.outer.solver_opts.tol_abs = 1e-16;
    opts.outer.solver_opts.tol_rel = 1e-16;
    opts.outer.solver_opts.A = @(x) A*x;
    opts.outer.solver_opts.b = b;
    opts.inner.solver_opts.gradient_mode = 'step';
    opts.inner.solver_opts.A = @(x) A_low*x;
    opts.inner.solver_opts.b = b;
    opts.outer.max_iter = 100;
    opts.inner.solver_opts.max_iter = 1;
    x0 = zeros(problem_size, 1);
    [~, info_dfp] = multifidelity_wrapper(fg, fg_low, x0, opts);

    %create error history for corr
    error_dfp = vecnorm(info_dfp.outer.iterHist' - x_ref)/norm(x_ref);


    %alternating correction dfp
    opts.warm.enabled = false;
    opts.inner.enabled = true;
    opts.outer.correction = true;
    opts.outer.storeFuncs = true;
    opts.inner.solver_opts.adaptive.enabled = true;
    opts.inner.solver_opts.adaptive.descent_parameter = 1;
    opts.outer.correction_opts.direction = 'secant';
    opts.outer.correction_opts.r = 20;
    opts.outer.correction_opts.update = 'dfp';
    opts.outer.correction_opts.skips = 'true';
    opts.outer.correction_opts.memory = 'compact';
    opts.outer.warm_correction = false;
    opts.outer.adaptive = 'none';
    %opts.outer.solver_opts.solver = 'bbpgd';
    %opts.inner.solver_opts.solver = 'bbpgd';
    opts.outer.solver_opts.tol_abs = 1e-16;
    opts.outer.solver_opts.tol_rel = 1e-16;
    opts.outer.solver_opts.A = @(x) A*x;
    opts.outer.solver_opts.b = b;
    opts.inner.solver_opts.gradient_mode = 'step';
    opts.inner.solver_opts.A = @(x) A_low*x;
    opts.inner.solver_opts.b = b;
    opts.outer.max_iter = 100;
    opts.inner.solver_opts.max_iter = 4;
    x0 = zeros(problem_size, 1);
    [~, info_dfp_adaptive] = multifidelity_wrapper(fg, fg_low, x0, opts);

    %create error history for corr
    error_dfp_adaptive = vecnorm(info_dfp_adaptive.outer.iterHist' - x_ref)/norm(x_ref);

    %alternating correction dfp
    %{
    opts.warm.enabled = false;
    opts.inner.enabled = true;
    opts.outer.correction = true;
    opts.outer.storeFuncs = true;
    opts.outer.correction_opts.direction = 'secant';
    opts.outer.correction_opts.r = 10;
    opts.outer.correction_opts.update = 'dfp';
    opts.outer.correction_opts.skips = 'true';
    opts.outer.correction_opts.memory = 'compact';
    opts.outer.warm_correction = false;
    opts.outer.adaptive = 'none';
    %opts.outer.solver_opts.solver = 'bbpgd';
    %opts.inner.solver_opts.solver = 'bbpgd';
    opts.outer.solver_opts.tol_abs = 1e-16;
    opts.outer.solver_opts.tol_rel = 1e-16;
    opts.outer.solver_opts.A = @(x) A*x;
    opts.outer.solver_opts.b = b;
    opts.inner.solver_opts.gradient_mode = 'step';
    opts.inner.solver_opts.A = @(x) A_low*x;
    opts.inner.solver_opts.b = b;
    opts.outer.max_iter = 100;
    opts.inner.solver_opts.max_iter = 2;
    x0 = zeros(problem_size, 1);
    [~, info_sr1] = multifidelity_wrapper(fg, fg_low, x0, opts);

    %create error history for corr
    error_sr1 = vecnorm(info_sr1.outer.iterHist' - x_ref)/norm(x_ref);

        %alternating correction dfp
    opts.warm.enabled = false;
    opts.inner.enabled = true;
    opts.outer.correction = true;
    opts.outer.storeFuncs = true;
    opts.outer.correction_opts.direction = 'secant';
    opts.outer.correction_opts.r = 10;
    opts.outer.correction_opts.update = 'dfp';
    opts.outer.correction_opts.skips = 'true';
    opts.outer.correction_opts.memory = 'compact';
    opts.outer.warm_correction = false;
    opts.outer.adaptive = 'none';
    %opts.outer.solver_opts.solver = 'bbpgd';
    %opts.inner.solver_opts.solver = 'bbpgd';
    opts.outer.solver_opts.tol_abs = 1e-16;
    opts.outer.solver_opts.tol_rel = 1e-16;
    opts.outer.solver_opts.A = @(x) A*x;
    opts.outer.solver_opts.b = b;
    opts.inner.solver_opts.gradient_mode = 'step';
    opts.inner.solver_opts.A = @(x) A_low*x;
    opts.inner.solver_opts.b = b;
    opts.outer.max_iter = 100;
    opts.inner.solver_opts.max_iter = 1;
    x0 = zeros(problem_size, 1);
    [~, info_bfgs] = multifidelity_wrapper(fg, fg_low, x0, opts);

    %create error history for corr
    error_bfgs = vecnorm(info_bfgs.outer.iterHist' - x_ref)/norm(x_ref);
    %}

    %plot the results
    semilogy(error_high, '-o', 'DisplayName', 'High Fidelity Only');
    hold on
    %semilogy(error_reuse, '-o', 'DisplayName', 'High Fidelity with Reuse');
    %semilogy(error_step, '-o', 'DisplayName', 'Multifidelity with Step');
    semilogy(error_dfp, '-o', 'DisplayName', 'Multifidelity with DFP');
    semilogy(error_dfp_adaptive, '-o', 'DisplayName', 'Multifidelity with DFP Adaptive');
    %semilogy(error_sr1, '-o', 'DisplayName', 'Multifidelity with SR1');
    %semilogy(error_bfgs, '-o', 'DisplayName', 'Multifidelity with BFGS');
    xlabel('Matvecs');
    ylabel('Error');
    legend('Location', 'best');
    title('Method Comparison');
    hold off
end