function update_comparisons(A, Ahat, b, Amat, Ahat_mat)
    % This function updates the comparisons between the original and approximated matrices
    % by computing the relative error between them.
    problem_size = size(Amat, 2);

    relative_matrix_error = norm(Amat - Ahat_mat)/norm(Amat);
    fprintf('Relative Matrix Error: %.2e\n', relative_matrix_error);
    %show the error in matrix norm of the three update methods
    opts.outer.correction = true;
    opts.outer.store_updates = true;
    opts.outer.low_update = 'sr1';
    opts.outer.adaptive = 'none';
    opts.inner.enabled = true;
    opts.inner.max_iter = 1; %take 1 low fidelity step between high fidelities

    %1 SR1
    %2 BFGS
    %3 DFP
    methods = {'sr1', 'bfgs', 'dfp'};
    error_matrices = zeros(3, 200);
    for i = 1:3
        opts.outer.low_update = methods{i};
        [~, info] = multi_fidelity_solver(A, Ahat, b, zeros(problem_size, 1), opts);
        for j = 2:info.outer.iter
            error_matrices(i, j) = norm(Amat - (Ahat_mat + info.outer.updates{j}))/norm(Amat);
        end
    end

    %now plot
    for i = 1:3
        semilogy(error_matrices(i, :), 'DisplayName', methods{i});
        hold on;
    end
    legend('show');
    xlabel('High Fidelity Matvecs');
    ylabel('Relative Error of Matrix Approximation');
    relative_matrix_error = norm(Amat - Ahat_mat)/norm(Amat);
    title_str = sprintf('Comparison of Updating Methods with Alternating Steps\nInitial Relative Matrix Error: %.2e', relative_matrix_error);
    title(title_str);
    saveas(gcf, 'update_comparions_alternating.svg')

    hold off;

    %Now do the same but only high fidelity
    opts.inner.enabled = false;
    error_matrices = zeros(3, 200);
    for i = 1:3
        opts.outer.low_update = methods{i};
        [~, info] = multi_fidelity_solver(A, Ahat, b, zeros(problem_size, 1), opts);
        for j = 2:info.outer.iter
            error_matrices(i, j) = norm(Amat - (Ahat_mat + info.outer.updates{j}))/norm(Amat);
        end
    end

    %now plot
    for i = 1:3
        semilogy(error_matrices(i, :), 'DisplayName', methods{i});
        hold on;
    end
    legend('show');
    xlabel('High Fidelity Matvecs');
    ylabel('Relative Error of Matrix Approximation');
    relative_matrix_error = norm(Amat - Ahat_mat)/norm(Amat);
    title_str = sprintf('Comparison of Updating Methods High Fidelity Only\nInitial Relative Matrix Error: %.2e', relative_matrix_error);
    title(title_str);
    saveas(gcf, 'update_comparions_high_fidelity.svg');

end