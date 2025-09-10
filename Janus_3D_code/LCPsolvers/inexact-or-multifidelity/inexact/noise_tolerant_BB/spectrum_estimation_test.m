function spectrum_estimation_test(problem_size)
    % This function runs a test for the spectrum estimation methods on the constructed test matrices.
    %Specifically, we see how well BB1 estimates the min and max spectrum in the presence of noise.
    % Add the folder and all its subdirectories
    addpath(genpath('matrix_utilities'));

    matrices = construct_test_matrices(problem_size); % Construct test matrices of size 100
    noise_levels = [1e0, 1e-2, 1e-4, 1e-6]; 
    opt_params.opt_tol = 1e-16;
    opt_params.method = 'none';
    opt_params.noise_control_param = 0.5; %somewhat high, no idea what good values are
    opt_params.line_search_param = 2; %doubles step length each time. I think this is reasonable
    opt_params.projection = false; %no projection for this test

    %For this test,
    b = 10*2*(rand(problem_size, 1) - 0.5); %no idea how to choose b
    max_iter = 200;
    matrix_names = fieldnames(matrices);
    for i = 1:4 %iterate over the matrices
        figure;
        matrix = matrices.(matrix_names{i}).matrix;
        A_noisy = noise_maker(@(x) matrix * x);
        
        % Store results for all noise levels
        mu_results = cell(4, 1);
        L_results = cell(4, 1);

        for j = 1:4 %iterate over all the noise levels
            noise = noise_levels(j);
            opt_params.noise = noise;
            fprintf('Testing matrix: %s with noise level: %e\n', matrix_names{i}, noise);
            [~, ~, ~, ~, mu_estimates, L_estimates, iters] = BBPGD_noise_tolerant(A_noisy, b, zeros(problem_size, 1), max_iter, opt_params);
            iters = iters - 1;
            mu_estimates = mu_estimates(1:iters);
            L_estimates = L_estimates(1:iters);
            
            % Calculate relative errors
            true_mu = 1/matrices.(matrix_names{i}).inverseNorm;
            true_L = matrices.(matrix_names{i}).matrixNorm;
            
            mu_relative_error = abs(mu_estimates - true_mu) / abs(true_mu);
            L_relative_error = abs(L_estimates - true_L) / abs(true_L);
            
            % Store results
            mu_results{j} = mu_relative_error;
            L_results{j} = L_relative_error;
        end
        
        % Plot mu estimates
        subplot(2, 1, 1);
        for j = 1:4
            semilogy(mu_results{j}, 'DisplayName', sprintf('Noise Level %e', noise_levels(j)));
            hold on;
        end
        title(sprintf('Min Spectrum (mu) Relative Error for %s', matrix_names{i}));
        xlabel('Iteration');
        ylabel('Relative Error (log scale)');
        legend('show');
        grid on;
        hold off;
        
        % Plot L estimates
        subplot(2, 1, 2);
        for j = 1:4
            semilogy(L_results{j}, 'DisplayName', sprintf('Noise Level %e', noise_levels(j)));
            hold on;
        end
        title(sprintf('Max Spectrum (L) Relative Error for %s', matrix_names{i}));
        xlabel('Iteration');
        ylabel('Relative Error (log scale)');
        legend('show');
        grid on;
        hold off;
        
        % Save the figure as PDF
        filename = sprintf('spectrum_estimation_%s.pdf', matrix_names{i});
        exportgraphics(gcf, filename, 'ContentType', 'vector');
        
        pause(10); % Pause to visualize the plot before moving to the next matrix
    end
end