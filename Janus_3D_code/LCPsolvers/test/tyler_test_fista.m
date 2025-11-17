addpath('../inexact-or-multifidelity/utilities');
addpath('../solvers');
% Just testing FISTA
problem_size = 1000;
matrices = construct_test_matrices(problem_size);
A = matrices.exp.matrix;
b = (randn(problem_size,1) - 1)*1000;
opts.A = @(x) A*x;
opts.b = b;
opts.errFcn = @(x) (1/2)*dot(min(x, opts.A(x) + opts.b), min(x, opts.A(x) + opts.b));
opts.max_iter = 1000;
opts.solver = 'fista';
opts.linesearch.budget = 1;
opts.linesearch.c1 = 0.1;
opts.linesearch.c2 = 0.1; 
opts.linesearch.tol = 1e-1;
opts.acceleration.method = 'adaptive';
opts.acceleration.restart = true; % Without restart the adaptive momentum methods seem much worse.
opts.acceleration.curvature = 'accelerated';
opts.stepSize.L = 3;
fg = @(x, Ax) quadraticLoss(x, A, b, Ax);

opts.stepSize.fwd = 'lipschitz';
opts.stepSize.bwd = 'uniform';
[x, info_adaptive_L_uniform] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'lipschitz';
opts.stepSize.bwd = 'opt';
[x, info_adaptive_L_opt] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'bb1';
opts.stepSize.bwd = 'uniform';
[x, info_adaptive_bb_uniform] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'bb1';
opts.stepSize.bwd = 'opt';
[x, info_adaptive_bb_opt] = fista(fg, zeros(problem_size,1), opts);

% now for known
opts.acceleration.method = 'known';
opts.acceleration.cond = 3;

opts.stepSize.fwd = 'lipschitz';
opts.stepSize.bwd = 'uniform';
[x, info_known_L_uniform] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'lipschitz';
opts.stepSize.bwd = 'opt';
[x, info_known_L_opt] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'bb1';
opts.stepSize.bwd = 'uniform';
[x, info_known_bb_uniform] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'bb1';
opts.stepSize.bwd = 'opt';
[x, info_known_bb_opt] = fista(fg, zeros(problem_size,1), opts);

[x, info_known] = fista(fg, zeros(problem_size,1), opts);

semilogy(info_adaptive_L_uniform.errHist, '-o', 'DisplayName', 'Adaptive FISTA Lipschitz and Uniform'); 
hold on;
semilogy(info_adaptive_L_opt.errHist, '-o', 'DisplayName', 'Adaptive FISTA Lipschitz and Optimal'); 
semilogy(info_adaptive_bb_uniform.errHist, '-o', 'DisplayName', 'Adaptive FISTA BB1 and Uniform'); 
semilogy(info_adaptive_bb_opt.errHist, '-o', 'DisplayName', 'Adaptive FISTA BB1 and Optimal'); 
hold on;
semilogy(info_known_L_uniform.errHist, '-x', 'DisplayName', 'Known FISTA Lipschitz and Uniform');
semilogy(info_known_L_opt.errHist, '-x', 'DisplayName', 'Known FISTA Lipschitz and Optimal');
semilogy(info_known_bb_uniform.errHist, '-x', 'DisplayName', 'Known FISTA BB1 and Uniform');
semilogy(info_known_bb_opt.errHist, '-x', 'DisplayName', 'Known FISTA BB1 and Optimal');
hold off;
title('FISTA on LCP Problem');
xlabel('Iteration');
ylabel('Objective Error');
legend('Location', 'Best');