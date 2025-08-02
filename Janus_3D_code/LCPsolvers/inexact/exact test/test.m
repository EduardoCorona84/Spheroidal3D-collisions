addpath("../matrix_utilities/")
problem_size = 1000;
matrices = construct_test_matrices(problem_size);
A = matrices.linear50.matrix;
A = (A + A') / 2; % make sure it's symmetric
b = (rand(problem_size, 1) - 1/2)*problem_size;

cvx_begin quiet
    cvx_precision best
    variable x_ref(problem_size)
    minimize 0.5 * dot(x_ref, A * x_ref) + dot(b, x_ref)
    subject to
        x_ref >= 0;
cvx_end

x_true = A\(-b);
objective = @(x) 0.5 * x' * A * x + b' * x;


fcnGrad = create_grad(A, b);
x0 = zeros(problem_size, 1);
opts.tol_rel = 1e-16;
opts.tol_abs = 1e-16;
opts.max_iter = 2000;
%[x, info] = proxQuasiNewton(fcnGrad, x0, opts);

opts.A = @(x) A*x;
opts.kappa.init = 'uniform';
opts.kappa.fwd = 'opt';
opts.kappa.bwd = 'opt';
[x_new, info_new] = projectedGradientDescent(fcnGrad, x0, opts);

%errors = vecnorm(info.iterHist' - x_ref)/norm(x_ref);


errors_new = vecnorm(info_new.iterHist' - x_true)/norm(x_true);
semilogy(errors_new)
hold on
errors_ref = vecnorm(info_new.iterHist' - x_ref)/norm(x_ref);
disp(objective(x_new));
disp(objective(x_ref));
disp(objective(x_new) < objective(x_ref));