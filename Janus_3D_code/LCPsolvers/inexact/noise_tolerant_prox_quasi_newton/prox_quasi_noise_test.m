addpath('../matrix_utilities/');

problem_size = 1000;
matrices = construct_test_matrices(problem_size);
noise_level = 1e-2;
b = (rand(problem_size, 1) - 0.5)*problem_size;
objective = @(x) 0.5 * x' * matrices.linear50.matrix * x + b' * x;
%find reference solution with cvx
cvx_begin quiet
    variable x_ref(problem_size)
    minimize 0.5 * dot(x_ref, matrices.linear50.matrix * x_ref) + dot(b, x_ref)
    subject to
        x_ref >= 0;
cvx_end

%optimal objective value
optimal_value = objective(x_ref);

fcnGrad = quadratic_gradient_maker(matrices.linear50.matrix, b, noise_level);
A_noisy = noise_maker(matrices.linear50.matrix);
opts.A = @(x) A_noisy(x, noise_level);
opts.kappa.fwd = 'opt'; 
opts.r = 20;
opts.tol_rel = 1e-20;
opts.tol_abs = 1e-20;

opts.noise_control = false;
[x_opt, info_opt] = prox_quasi_newton_noise_tolerant(fcnGrad, zeros(problem_size, 1), opts);

opt_evals = zeros(info_opt.iter,1);
for i = 1:info_opt.iter
    opt_evals(i) = objective(info_opt.iterHist(i, :)');
end
opt_errors = abs(opt_evals - optimal_value);

opts.noise_control = true;
opts.noise = noise_level;
opts.noise_control_type = 'simple';
opts.noise_control_parameter = 0.1;
opts.noise_control_max_iter = 100;
opts.noise_control_line_search_param = 2;

[x_simple, info_simple] = prox_quasi_newton_noise_tolerant(fcnGrad, zeros(problem_size, 1), opts);
simple_evals = zeros(info_simple.iter,1);
for i = 1:info_simple.iter
    simple_evals(i) = objective(info_simple.iterHist(i, :)');
end
simple_errors = abs(simple_evals - optimal_value);

opts.noise_control_type = 'proximal';
[x_prox, info_prox] = prox_quasi_newton_noise_tolerant(fcnGrad, zeros(problem_size, 1), opts);

prox_evals = zeros(info_prox.iter,1);
for i = 1:info_prox.iter
    prox_evals(i) = objective(info_prox.iterHist(i, :)');
end
prox_errors = abs(prox_evals - optimal_value);

semilogy(1:info_opt.iter, opt_errors, 'DisplayName', 'No Noise Control');
hold on;
semilogy(1:info_simple.iter, simple_errors, 'DisplayName', 'Simple Noise Control');
semilogy(1:info_prox.iter, prox_errors, 'DisplayName', 'Proximal Noise Control');
xlabel('Iteration');
ylabel('Objective Value Error');
title('Noise Tolerant Proximal Quasi-Newton');
legend('Location', 'best');
grid on;
hold off;



