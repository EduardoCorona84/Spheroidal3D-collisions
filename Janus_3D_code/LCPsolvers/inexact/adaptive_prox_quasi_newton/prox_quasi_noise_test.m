addpath('../matrix_utilities/');
addpath('../gmres_utilities/');

problem_size = 1000;
matrices = construct_test_matrices(problem_size);
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
x_optimal_value = objective(x_ref);
restart = 2;
gmres_max_iter = 200;

fcnGrad = construct_gmres_gradient(matrices.linear50.inverse, b, matrices.linear50.matrixNorm, restart, gmres_max_iter);
opts.A = matrix_construct_gmres(matrices.linear50.inverse, matrices.linear50.matrixNorm, restart, gmres_max_iter);
opts.kappa.fwd = 'uniform'; 
opts.kappa.bwd = 'uniform'; 
opts.kappa.init = 'uniform';
opts.noise_adaptive = true;
opts.noise = 1e-2;
opts.noise_control_parameter = 1;
opts.noise_control_reduction_parameter = 0.01;
opts.r = 20;
opts.tol_rel = 1e-20;
opts.tol_abs = 1e-20;

[x, info] = prox_quasi_newton_adaptive_noise_control(fcnGrad, zeros(problem_size, 1), opts);
gmres_iters = cumsum(info.gmres_iters(1:info.iter));

%calculate error 


objective_values = zeros(info.iter, 1);
for i = 1:info.iter
    objective_values(i) = objective(info.iterHist(i,:)');
end
objective_errors_cvx = abs(objective_values - x_optimal_value)/abs(x_optimal_value);
objective_errors_us = abs(objective_values - objective(x))/abs(objective(x));
semilogy(gmres_iters, objective_errors_cvx, 'o-');
hold on;
semilogy(gmres_iters, objective_errors_us, 'x-');
xlabel('GMRES Iterations');
ylabel('Relative Error');
title('Adaptive Prox Quasi-Newton with Noise Control');
grid on;
disp(objective(x))
disp(x_optimal_value);
disp(objective(x) - x_optimal_value);
if (objective(x) > x_optimal_value)
    disp('cvx better')
else
    disp('us better')
end
