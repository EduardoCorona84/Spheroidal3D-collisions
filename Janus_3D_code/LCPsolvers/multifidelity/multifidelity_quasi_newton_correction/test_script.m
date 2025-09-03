%in this script we want to test the current code versus nic's implementation of prox quasi-Newton (all high) and then we want to test versus my previous implementation of alternating low and high fidelity with correction (and just does it run)

%add code to generate test matrices
addpath('../../inexact/matrix_utilities/')

%add nic's solvers
addpath('../../solvers/')

problem_size = 100;

matrices = construct_test_matrices(problem_size);
%we will just test everything on this matrix
the_matrix = matrices.linear50.matrix;

b = (rand(problem_size, 1) - 1/2)*10;

%generate function handle
fg = create_fg(the_matrix, b);

%dummy low fidelity
fg_low = create_fg(eye(problem_size), b);
%1. First Case: Compare High against High (only high fidelity evaluations)

%create opts for me
opts.inner.enabled = false;
x0 = zeros(problem_size, 1);
[x_new, info_new] = multifidelity_wrapper(x0, fg, fg_low, opts);

%create opts for nic

[x_nic, info_nic] = proxQuasiNewton(x0, fg, opts_nic);



%2. Second Case: Compare Alternating with Correction (Old) versus Alternating with Correction (New)

%add old code
addpath('../LCPsolvers/multifidelity/')

[x_old, info_old] = multifidelity_solver(x0, fg, opts_nic);

[x_]