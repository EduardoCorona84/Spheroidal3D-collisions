addpath('../inexact-or-multifidelity/utilities/');
addpath('../solvers/');
n = 100;
matrices = construct_test_matrices(n);
b = 2*(randn(n, 1) - 1/2);

fg = create_fg(matrices.linear50.matrix, b);
opts.solver = 'apgd';
opts.b = b;
opts.A = @(x) matrices.linear50.matrix * x;
[x, info] = accelerated_prox_grad(fg, zeros(n, 1), opts);