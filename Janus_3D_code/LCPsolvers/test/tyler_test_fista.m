
addpath('../solvers');
% Just testing FISTA
problem_size = 1000;
U = randn(problem_size, problem_size);
[Q, ~] = qr(U);
D = diag(linspace(1, 5, problem_size));
matrix = Q*D*Q';
A = matrix;
disp(eig(A));
b = (randn(problem_size,1) - 1)*10;
opts.A = @(x) A*x;
opts.b = b;
opts.errFcn = {@(x) (1/2)*dot(min(x, opts.A(x) + opts.b), min(x, opts.A(x) + opts.b)), @(x) (1/2)*opts.A(x)'*x + opts.b'*x};
opts.max_iter = 1000;
opts.solver = 'fista';
opts.storeIts = true;
opts.linesearch.budget = 1;
opts.linesearch.c1 = 0.1;
opts.linesearch.c2 = 0.1; 
opts.linesearch.tol = 1e-1;
opts.acceleration.method = 'adaptive';
opts.acceleration.cond = 5;
opts.acceleration.restart = true; % Without restart the adaptive momentum methods seem much worse.
opts.acceleration.curvature = 'feasible';
opts.acceleration.direction = 'feasible';
opts.stepSize.L = 5;
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


opts.acceleration.curvature = 'accelerated';

opts.stepSize.fwd = 'lipschitz';
opts.stepSize.bwd = 'uniform';
[x, info_adaptive_L_uniform_a] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'lipschitz';
opts.stepSize.bwd = 'opt';
[x, info_adaptive_L_opt_a] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'bb1';
opts.stepSize.bwd = 'uniform';
[x, info_adaptive_bb_uniform_a] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'bb1';
opts.stepSize.bwd = 'opt';
[x, info_adaptive_bb_opt_a] = fista(fg, zeros(problem_size,1), opts);

%now for known
opts.acceleration.method = 'known';
opts.acceleration.restart = false;
opts.acceleration.curvature = 'feasible';

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


opts.acceleration.curvature = 'accelerated';
opts.stepSize.fwd = 'lipschitz';
opts.stepSize.bwd = 'uniform';
[x, info_known_L_uniform_a] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'lipschitz';
opts.stepSize.bwd = 'opt';
[x, info_known_L_opt_a] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'bb1';
opts.stepSize.bwd = 'uniform';
[x, info_known_bb_uniform_a] = fista(fg, zeros(problem_size,1), opts);
opts.stepSize.fwd = 'bb1';
opts.stepSize.bwd = 'opt';
[x, info_known_bb_opt_a] = fista(fg, zeros(problem_size,1), opts);

% adaptive with feasible
semilogy(info_adaptive_L_uniform.errHist(:, 1), '-o', 'DisplayName', 'Adaptive FISTA Lipschitz and Uniform'); 
hold on;
semilogy(info_adaptive_L_opt.errHist(:, 1), '-o', 'DisplayName', 'Adaptive FISTA Lipschitz and Optimal'); 
semilogy(info_adaptive_bb_uniform.errHist(:, 1), '-o', 'DisplayName', 'Adaptive FISTA BB1 and Uniform and feasible'); 
semilogy(info_adaptive_bb_opt.errHist(:, 1), '-o', 'DisplayName', 'Adaptive FISTA BB1 and Optimal and feasible');  
% adaptive with accelerated
semilogy(info_adaptive_bb_uniform_a.errHist(:, 1), '-*', 'DisplayName', 'Adaptive FISTA BB1 and Uniform and accelerated curvature'); 
semilogy(info_adaptive_bb_opt_a.errHist(:, 1), '-*', 'DisplayName', 'Adaptive FISTA BB1 and Optimal and accelerated curvature');  

%known with feasible

semilogy(info_known_L_uniform.errHist(:, 1), '-x', 'DisplayName', 'Known FISTA Lipschitz and Uniform');
semilogy(info_known_L_opt.errHist(:, 1), '-x', 'DisplayName', 'Known FISTA Lipschitz and Optimal');
semilogy(info_known_bb_uniform.errHist(:, 1), '-x', 'DisplayName', 'Known FISTA BB1 and Uniform');
semilogy(info_known_bb_opt.errHist(:, 1), '-x', 'DisplayName', 'Known FISTA BB1 and Optimal');
%known with accelerated
semilogy(info_known_bb_uniform_a.errHist(:, 1), '-.', 'DisplayName', 'Known FISTA BB1 and Uniform and accelerated');
semilogy(info_known_bb_opt_a.errHist(:, 1), '-.', 'DisplayName', 'Known FISTA BB1 and Optimal and accelerated');

hold off;
title('FISTA on LCP Problem');
xlabel('Iteration');
ylabel('Objective Error');
legend('Location', 'Best');
