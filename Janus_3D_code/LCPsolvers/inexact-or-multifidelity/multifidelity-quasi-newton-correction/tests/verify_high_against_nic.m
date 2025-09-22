%In this script, we want to verify that our high fidelity only prox solver and nic's prox quasi solver have the same outputs.

%add nic's code
addpath('../../../solvers/');

%add my code
addpath('../src/')

%add required utilities
addpath('../../utilities/')


problem_size = 100;

%generate test matrices
matrices = construct_test_matrices(problem_size);
b = (rand(problem_size, 1) - 1/2)*10;

%dummy low fidelity function as we are only testing high
fg_low = create_fg(eye(problem_size), b);

%now go through each test matrix and compare output on mine and nic's code
spectrums = fieldnames(matrices);

for i = 1:length(matrices)
    A = matrices.(spectrums{i}).matrix;
    fg = create_fg(A, b);
    
    %set nic's opts and run code
    opts_nic.solver = 'prox';
    opts_nic.max_iter = 100;
    opts_nic.A = @(x) A*x;
    opts_nic.b = b;
    opts_nic.storeIts = true;
    x0 = zeros(problem_size, 1);
    [~, info_nic] = proxQuasiNewton(fg, x0, opts_nic);

    %set my opts and run code
    opts.warm.enabled = false;
    opts.inner.enabled = false;
    opts.outer.correction = false;
    opts.outer.max_iter = 100;
    opts.outer.solver_opts.A = @(x) A*x;
    opts.outer.solver_opts.b = b;
    x0 = zeros(problem_size, 1);
    [~, info_me] = multifidelity_wrapper(fg, fg_low, x0, opts);

    %find min iters from the two methods
    min_iters = min(size(info_nic.iterHist, 1), size(info_me.outer.iterHist, 1));
    %now double check that they have the same outputs7
    assert(isequal(info_nic.iterHist(1:min_iters, :), info_me.outer.iterHist(1:min_iters, :)), 'iterHist do not match');

end

%if successful, print a message
disp('All tests passed!');




