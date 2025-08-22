%We want to compare the iterates of nic's solver and mine, for the same problem.
addpath('../test_data/')
fname = '4x4x4_amphi_500';
load([fname '.mat'], ...
    'A_list', 'A_diag_list', 'b_list');

A = A_list{200};
A = (A + A')/2;
%Ahat = A_diag_list{200};
%Ahat = (Ahat + Ahat')/2;

b = b_list{150};
addpath('../solvers/')
addpath('../inexact/matrix_utilities/')
% Reference solution with no noise
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

fcnGrad = create_nic_gradient(A, b);
[~, nic_info] = proxQuasiNewton(fcnGrad, zeros(size(b)), opts);
nic_iters = nic_info.iterHist';
disp(size(nic_iters));
clear opts;

opts.outer.projected = true;
opts.outer.warm_start.enabled = false;
opts.outer.warm_start.only = false;
opts.inner.enabled = false;
opts.outer.correction = false;
Afunc = @(x) A*x;
Ahat_func = @(x) Ahat*x;
[~, my_info] = multi_fidelity_solver(Afunc, Ahat_func, b, zeros(size(b)), opts);
my_iters = my_info.outer.x_iters(:, 1:size(nic_iters, 2));
disp(size(my_iters));

error_iters = vecnorm(nic_iters -my_iters)./vecnorm(nic_iters);
disp(size(vecnorm(nic_iters)))
disp(size(error_iters));
figure;
plot(error_iters, 'o-', 'DisplayName', 'Nic Solver');
hold on
xlabel('Iteration');
ylabel('Objective Value');
legend();
title('Comparison of Solver Iterations');