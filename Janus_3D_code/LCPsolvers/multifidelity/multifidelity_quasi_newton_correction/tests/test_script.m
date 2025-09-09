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

%create opts for nic
%same solver type for both
opts_nic.max_iter = 10;
opts_nic.solver = 'prox';
opts_nic.A = @(x) the_matrix*x;
opts_nic.b = b;
opts_nic.storeIts = true;
[x_nic, info_nic] = proxQuasiNewton(fg, x0, opts_nic);

opts.warm.enabled = false;
opts.inner.enabled = true;
opts.outer.solver_opts.A = @(x) the_matrix*x;
opts.outer.solver_opts.b = b;
opts.outer.max_iter = 10;
x0 = zeros(problem_size, 1);
[x_new, info_new] = multifidelity_wrapper(fg, fg_low, x0, opts);





%2. Second Case: Compare Alternating with Correction (Old) versus Alternating 

%add old code
