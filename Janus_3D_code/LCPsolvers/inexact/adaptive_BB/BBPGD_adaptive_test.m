function BBPGD_adaptive_test(problem_size) 

    addpath(genpath('../matrix_utilities'));
    addpath(genpath('../gmres_utilities'));

    matrices = construct_test_matrices(problem_size);
    noise_levels = [1e0, 1e-2, 1e-4, 1e-6]; 
    
    % common parameters
    opt_params.opt_tol = 1e-12;
    opt_params.noise_control_param = 0.5;
    opt_params.noise_shrink_param = 0.5;
    opt_params.projection = true;
    opt_max_iter = 500;


    %gmres parameters
    restart = 5;
    gmres_max_iter = 200;

    b = 10*2*(rand(problem_size, 1) - 0.5);
    matrix_names = fieldnames(matrices);
    
    figure;
    for i = 1:length(matrix_names)
        %find subplot
        subplot(2, 2, i);
        
        %make matrix
        matrix_data = matrices.(matrix_names{i});
        %Construct GMRES, we want to compute matvecs of A, which means we pass the inverse to gmres.
        A_gmres = matrix_construct_gmres(matrix_data.inverse, matrix_data.matrixNorm, restart, gmres_max_iter);

        x_errors = cell(5, 1);
        noise_kkt_errors = cell(5, 1);
        cum_gmres = cell(5, 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %Find cvx true solution
        cvx_begin quiet
                variable x_true(problem_size)
                minimize 0.5*dot(x_true, matrix_data.matrix * x_true) + dot(x_true, b)
                subject to
                0 <= x_true
        cvx_end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for j = 1:4
            opt_params.adaptive = false;
            opt_params.noise = noise_levels(j);
            [~, x_iters, gmres_iters, kkt_errors, ~, ~, iters] = BBPGD_adaptive(A_gmres, b, zeros(problem_size, 1), opt_max_iter, opt_params);

            x_errors{j} = vecnorm(x_iters(1:iters - 1) - x_true)/norm(x_true);
            noise_kkt_errors{j} = kkt_errors(1:iters - 1);
            cum_gmres{j} = cumsum(gmres_iters(1:iters - 1));
        end

        opt_params.adaptive = true;
        opt_params.noise = 1e-1;
        [~, ~, gmres_iters, kkt_errors, ~, ~, iters] = BBPGD_adaptive(A_gmres, b, zeros(problem_size, 1), opt_max_iter, opt_params);

        x_errors{5} = vecnorm(x_iters(1:iters - 1) - x_true)/norm(x_true);
        noise_kkt_errors{5} = kkt_errors(1:iters - 1);
        cum_gmres{5} = cumsum(gmres_iters(1:iters - 1));
        % Plot
        for k = 1:4
            semilogy(cum_gmres{k}, noise_kkt_errors{k}, 'DisplayName', sprintf('Noise level: %e', noise_levels(k)), 'LineWidth', 1.5);
            hold on;
        end
        semilogy(cum_gmres{5}, noise_kkt_errors{5}, 'DisplayName', 'Adaptive', 'LineWidth', 1.5, 'LineStyle', '--');
        title(sprintf('Matrix: %s', matrix_names{i}));
        xlabel('Cumulative GMRES Iterations');
        ylabel('KKT Error');
        legend('Location', 'best');
        grid on;
        hold off;
    end
    sgtitle('BBPGD Adaptive Test Results');
    legend('show');
    filename = sprintf('BBPGD_adaptive_test_%d.pdf', problem_size);
    exportgraphics(gcf, filename, 'ContentType', 'vector');
    pause(4);

    figure;
    for i = 1:length(matrix_names)
        %find subplot
        subplot(2, 2, i);
        for k = 1:4
            semilogy(1:length(noise_kkt_errors{k}), noise_kkt_errors{k}, 'DisplayName', sprintf('Noise level: %e', noise_levels(k)), 'LineWidth', 1.5);
            hold on;
        end
        semilogy(1:length(noise_kkt_errors{5}), noise_kkt_errors{5}, 'DisplayName', 'Adaptive', 'LineWidth', 1.5, 'LineStyle', '--');
        title(sprintf('Matrix: %s', matrix_names{i}));
        xlabel('Iteration');
        ylabel('KKT Error');
        legend('Location', 'best');
        grid on;
        hold off;
    end
    sgtitle('BBPGD Adaptive Test Results - KKT Error');
    legend('show'); 
    filename = sprintf('BBPGD_adaptive_test_kkt_%d.pdf', problem_size);
    exportgraphics(gcf, filename, 'ContentType', 'vector');

    figure;
    for i = 1:length(matrix_names)
        %find subplot
        subplot(2, 2, i);
        for k = 1:4
            semilogy(1:length(x_errors{k}), x_errors{k}, 'DisplayName', sprintf('Noise level: %e', noise_levels(k)), 'LineWidth', 1.5);
            hold on;
        end
        semilogy(1:length(x_errors{5}), x_errors{5}, 'DisplayName', 'Adaptive', 'LineWidth', 1.5, 'LineStyle', '--');
        title(sprintf('Matrix: %s', matrix_names{i}));
        xlabel('Iteration');
        ylabel('Relative Error');
        legend('Location', 'best');
        grid on;
        hold off;
    end
    sgtitle('BBPGD Adaptive Test Results - Relative Error');
    legend('show'); 
    filename = sprintf('BBPGD_adaptive_test_relerror_%d.pdf', problem_size);
    exportgraphics(gcf, filename, 'ContentType', 'vector');
end