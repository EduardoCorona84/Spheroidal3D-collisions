function testing_prox(problem_size)

   
    %This function tests the multifidelity prox solver with fixed low fidelity steps
    matrices = construct_test_matrices(problem_size);
    A = matrices.linear50.matrix;
    A_func = @(x) A*x;
    Ahat = matrices.linear50.matrix + diag(rand(problem_size, 1));
    Ahat_func = @(x) Ahat*x;
    b = (rand(problem_size, 1) -1/2)*problem_size;

    opts.outer.store_f = true;

    %take 1 low between high
    [x, info_multi] = multi_fidelity_solver(A_func, Ahat_func, b, zeros(problem_size, 1), opts);

    %take only high
    opts.inner.enabled = false;
    [x, info_single] = multi_fidelity_solver(A_func, Ahat_func, b, zeros(problem_size, 1), opts);

    ref = min(info_single.outer.f_iters);
    semilogy(abs((info_multi.outer.f_iters - ref)/(ref)));
    hold on;
    semilogy(abs((info_single.outer.f_iters - ref)/(ref)), 'r--');
    legend('Multi Fidelity', 'Single Fidelity');
    title('Multi-Fidelity Solver Performance');
    xlabel('Iteration');
    ylabel('Function Evaluations');
    hold off;

end