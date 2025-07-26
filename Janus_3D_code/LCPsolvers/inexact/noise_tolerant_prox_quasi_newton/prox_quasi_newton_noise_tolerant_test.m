function prox_quasi_newton_noise_tolerant_test(problem_size)
    
    addpath(genpath('../matrix_utilities'));
    matrices = construct_test_matrices(problem_size);
    noise_levels = [1e0, 1e-2, 1e-4, 1e-6]; 
    opts.noise_control = false;
    opts.noise_control_line_search = 'line search';
    opts.noise_control_lengthening = 'proximal';
    opts.noise_control_line_search_parameter = 1.25;
    opts.noise_control_parameter = 1;
    opts.noise = 1e-1;
    opts.r = 20;
    opts.max_iter = 1000;
    opts.tol_rel = 1e-12;
    opts.tol_abs = 1e-12;


    A = matrices.linear50.matrix;
    disp(class(A));
    b = (rand(problem_size, 1) -1/2)*problem_size;
    A_noisy = noise_maker(@(x) A*x);
    A_specific = @(x) A_noisy(x, 1e-1);
    gradient = construct_noisy_gradient(A_specific, b);
    opts.Q = A_specific;
    [x, iterProx, errStructProx] = prox_quasi_newton_noise_tolerant(gradient, zeros(problem_size, 1), opts);
    cvx_begin
        variable x_true(problem_size)
        minimize(0.5 * quad_form(x_true, A) + b' * x_true)
        subject to
            x_true >= 0
    cvx_end


    errors_prox = vecnorm(errStructProx.xk - x_true)/norm(x_true);
    plot(1:iterProx, errors_prox(1:iterProx), 'o-');
    hold on;

    opts.noise_control_lengthening = 'simple';
    [x, iterSimple, errStructSimple] = prox_quasi_newton_noise_tolerant(gradient, zeros(problem_size, 1), opts);
    errors_simple = vecnorm(errStructSimple.xk - x_true)/norm(x_true);
    plot(1:iterSimple, errors_simple(1:iterSimple), 'x-');
    %common parameters

end