addpath("../matrix_utilities");
addpath("../gmres_utilities");
problem_size = 100;
matrices = construct_test_matrices(problem_size);
rng(42);
b = (rand(problem_size, 1) - 1/2)*100;
matrix = diag(matrices.flat.eigenvalues);
opts.tol_rel = 1e-30;
opts.tol_abs = 1e-30;
opts.storeIts = true;
opts.max_iter = 400;
opts.noise_control.enabled = false;
opts.noise_adaptive.quad = false;
opts.noise_control.type = 'none';
opts.r = 10;
opts.kappa.init = 'uniform';
opts.kappa.fwd = 'uniform';
opts.kappa.bwd = 'uniform';

objective = @(x) (1/2)*x'*matrix*x + b'*x;
fcnGrad = create_gradient(matrix, b);
[~, info_ref] = proxQuasiNewton(fcnGrad, zeros(size(b)), opts);
objective_ref = zeros(info_ref.iter, 1);
for i = 1:info_ref.iter
    objective_ref(i) = objective(info_ref.iterHist(i, :)');
end

restart = 5;
gmres_max_iter = 200;
ref_val = min(objective_ref);

fcnGrad = construct_gmres_gradient(diag(1./diag(matrix)), b, max(diag(matrix)), 1/min(diag(matrix)), restart, gmres_max_iter);

opts.noise_control.noise = 1e0;
[~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
error_cur = zeros(info.iter, 1);
for j = 1:info.iter
    error_cur(j) = abs(objective(info.iterHist(j, :)') - ref_val)/abs(ref_val);
end
iter_data = cumsum(info.gmres_iters(1:info.iter));
semilogy(iter_data, error_cur, 'DisplayName', 'Fixed Noise Level: 1e-2');
semilogy(1:info.iter, error_cur, 'o', 'HandleVisibility','off');

disp(2);

