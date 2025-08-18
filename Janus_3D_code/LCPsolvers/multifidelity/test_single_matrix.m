function test_single_matrix(A, Ahat, b, Amat, Ahat_mat)
    %This function tests a single matrix and its approximation via various methods of multifidelity solve:
    %1. High Fidelity Only
    %2. Only low fidelity (warm start, 20 iters)
    %3. Warm Start into high fidelity
    %4. Alternating low and high fidelity with no correction
    %5. Alternating low and high fidelity with correction (SR1)
    %6. Alternating low and high fidelity with correction and warm start (20 iters).
    %7. Alternating low and high fidelity with high fidelity once criterion is reached (not sufficient decrease, similar) and correction and warm start (20 iters)
    %8. low fidelity warm start (20 iters), into alternating high and low fidelity, but now initial is 4 low to 1 high, halved when criteria is not met.

    %first find the reference value of the true argmin (we just do this by running a ton of iterations with only high fidelity, just reduces to prox quasi-Newton.) Checked against cvx's solution for a sanity check. (Note cvx's solution often has a worse objective value when the matrix is truly symmetric, so we use the prox quasi-Newton solver instead.)

    problem_size = size(b, 1);
    error_cell = cell(8, 1);

    %sanity check cvx
    cvx_begin quiet
        variable x_cvx(problem_size)
        minimize(1/2 * quad_form(x_cvx, Amat) + b'*x_cvx)
        subject to
            x_cvx >= 0
    cvx_end
    
 

    objective = @(x) (1/2)*x'*A(x) + b'*x;
    
    %1. High Fidelity Only
    opts.outer.warm_start.enabled = false;
    opts.outer.warm_start.only = false;
    opts.inner.enabled = false;
    opts.outer.correction = false;
    [~, info] = multi_fidelity_solver(A, Ahat, b, zeros(problem_size, 1), opts);
    %find x_ref
    objective_values = zeros(1, info.outer.iter);
    for i = 1:info.outer.iter
        objective_values(i) = objective(info.outer.x_iters(:, i));
    end
    [min_val, min_ind] = min(objective_values);
    if (min_val > objective(x_cvx))
        error('CVX solution is the minimum of the objective function, check your problem setup.');
    end
    x_ref = info.outer.x_iters(:, min_ind);

    error_cell{1} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);
    
    %{
    sanity check against cvx
    cvx_begin quiet
        variable x_cvx(problem_size)
        minimize(1/2 * quad_form(x_cvx, A(x_cvx)) + b'*x_cvx)
        subject to
            x_cvx >= 0
    cvx_end
    %}

    %2. Only low fidelity (warm start, 100 iters)
    %opts.outer.warm_start.enabled = true;
    %opts.outer.warm_start.max_iter = 100;
    %opts.outer.warm_start.only = true;
    %[x_warm, ~] = multi_fidelity_solver(A, Ahat, b, zeros(problem_size, 1), opts);

    cvx_begin 
        variable x_warm(problem_size)
        minimize(1/2 * quad_form(x_warm, Ahat_mat) + b'*x_warm)
        subject to
            x_warm >= 0
    cvx_end

    %compute error for the final low fidelity warm start iterate
    error_cell{2} = norm(x_warm - x_ref)/norm(x_ref);


    %3. Warm Start into High Fidelity Only, can just do this by passing x_warm
    opts.outer.warm_start.enabled = false;
    opts.outer.warm_start.only = false;
    opts.outer.correction = false;
    opts.inner.enabled = false;
    [~, info] = multi_fidelity_solver(A, Ahat, b, x_warm, opts);

    %compute error for all the high fidelity iterates
    error_cell{3} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);

    %4. Alternating low and high fidelity with no correction
    %disable adaptive high fidelity switching
    opts.outer.adaptive = 'none';
    opts.inner.enabled = true;
    opts.inner.max_iter = 1; %take 1 low fidelity step between high fidelities

    [~, info] = multi_fidelity_solver(A, Ahat, b, zeros(problem_size, 1), opts);

    %compute error for all the high fidelity iterates.
    error_cell{4} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);

    %5. Alternating low and high fidelity with correction (SR1)
    opts.outer.correction = true;
    opts.outer.store_updates = true;
    opts.outer.low_update = 'sr1';

    [~, info] = multi_fidelity_solver(A, Ahat, b, zeros(problem_size, 1), opts);

    %{
    error_matrices = zeros(1, info.outer.iter);
    for i = 2:info.outer.iter
        error_matrices(i) = norm(Amat - (Ahat_mat + info.outer.updates{i}))/norm(Amat);
    end

    semilogy(error_matrices);
    %}


    %compute error for all the high fidelity iterates.
    error_cell{5} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);

    %6. Alternating low and high fidelity with correction and warm start (100 iters).
    %We just do this by passing x_warm as the starting iterate

    [~, info] = multi_fidelity_solver(A, Ahat, b, x_warm, opts);

    %compute error for all the high fidelity iterates.
    error_cell{6} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);

    %7. Alternating low and high fidelity with correction, warm start, and adaptive high fidelity switching
    opts.outer.adaptive = 'high';

    [~, info] = multi_fidelity_solver(A, Ahat, b, x_warm, opts);

    %compute error for all the high fidelity iterates.
    error_cell{7} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);

    %8. Four low fidelities for every 1 high fidelity, with adaptive halving of low fidelity evals once criterion is reached. With warm start and correction
    opts.outer.adaptive = 'halving';
    opts.inner.max_iter = 4; %take 4 low fidelity steps between high fidelity

    [~, info] = multi_fidelity_solver(A, Ahat, b, x_warm, opts);

    %compute error for all the high fidelity iterates.
    error_cell{8} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);

    %Now plot the error results

    methods = {'Only High Fidelity', 'Only Low Fidelity', 'Warm Start High Fidelity', 'Alternating Low/High', 'Alternating Low/High (SR1)', 'Warm Start Alternating Low/High (SR1)', 'Warm Start Adaptive High Fidelity (SR1)', 'Warm Start Halving High Fidelity (SR1)'};

    
    for i = 1:8
        if i == 2
            semilogy(error_cell{i}*ones(1, 200), 'DisplayName', methods{i});
            hold on;
        end
        [val, ind] = min(error_cell{i});
        if val < 1e-10
            semilogy(error_cell{i}(1:90), 'DisplayName', methods{i});
            hold on;
        else
            semilogy(error_cell{i}, 'DisplayName', methods{i});
            hold on;
        end
    end
    hold off;
    legend show;
    xlabel('High Fidelity Matvecs')
    ylabel('Relative Error')
    relative_matrix_error = norm(Amat - Ahat_mat)/norm(Amat);
    title_str = sprintf('Multifidelity Solver Error Comparison\nRelative Matrix Error: %.2e', relative_matrix_error);
    title(title_str);
    saveas(gcf, 'multifidelity_sr1_comparison.svg')



end