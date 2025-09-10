function method_comparison_test(problem_size)
    % This function compares different noise-tolerant methods across matrices and noise levels
    % Add the folder and all its subdirectories
    addpath(genpath('matrix_utilities'));
    
    matrices = construct_test_matrices(problem_size);
    noise_levels = [1e0, 1e-2, 1e-4, 1e-6]; 
    methods = {'none', 'affine true', 'affine estimate', 'lengthening'};
    
    % Common parameters
    opt_params.opt_tol = 1e-16;
    opt_params.noise_control_param = 0.5;
    opt_params.line_search_param = 1.1;
    opt_params.projection = false;
    
    b = 10*2*(rand(problem_size, 1) - 0.5);
    max_iter = 1000;
    matrix_names = fieldnames(matrices);
    
    % Loop over each matrix
    for i = 1:length(matrix_names)
        figure;
        matrix_data = matrices.(matrix_names{i});
        matrix = matrix_data.matrix;
        A_noisy = noise_maker(@(x) matrix * x);
        
        % True solution for error calculation
        x_true = matrix \ (-b);
        
        % Loop over each noise level (create subplots)
        for j = 1:4
            % Create subplots with space at top for title and bottom for legend
            h = subplot(2, 2, j);
            % Adjust position to leave room at top and bottom
            pos = get(h, 'Position');
            pos(2) = pos(2) + 0.08;  % Move up by 8%
            pos(4) = pos(4) - 0.15; % Reduce height by 15% (for both title and legend space)
            set(h, 'Position', pos);
            
            noise = noise_levels(j);
            opt_params.noise = noise;
            
            % Store results for all methods
            method_errors = cell(4, 1);
            
            % Loop over each method
            for k = 1:4
                opt_params.method = methods{k};
                
                % For 'affine true' method, we need to pass true spectrum values
                if strcmp(methods{k}, 'affine true')
                    opt_params.L = matrix_data.matrixNorm;
                    opt_params.mu = 1/matrix_data.inverseNorm;
                end
                
                fprintf('Testing matrix: %s, noise: %e, method: %s\n', ...
                    matrix_names{i}, noise, methods{k});

                [~, x_iters, ~, ~, ~, ~, iters] = BBPGD_noise_tolerant(A_noisy, b, zeros(problem_size, 1), max_iter, opt_params);

                % Calculate error at each iteration
                errors = zeros(iters-1, 1);
                for iter = 1:iters-1
                    errors(iter) = norm(x_iters(:, iter) - x_true) / norm(x_true);
                end
                
                method_errors{k} = errors;
            end
            
            % Plot all methods for this noise level
            for k = 1:4
                loglog(method_errors{k}, 'LineWidth', 1.5);
                hold on;
            end
            
            title(sprintf('Noise Level: %e', noise_levels(j)));
            xlabel('Iteration (log scale)');
            ylabel('Relative Error (log scale)');
            grid on;
            hold off;
        end
        
        % Add overall title for the figure with better positioning
        sgtitle(sprintf('Method Comparison for Matrix: %s', matrix_names{i}), 'FontSize', 14);
        
        % Create a legend positioned in the bottom space we created
        lgd = legend(methods, 'Position', [0.2, 0.02, 0.6, 0.08]);
        lgd.FontSize = 10;
        lgd.Orientation = 'horizontal';
        
        % Save the figure
        filename = sprintf('method_comparison_%s.pdf', matrix_names{i});
        exportgraphics(gcf, filename, 'ContentType', 'vector');
        
        pause(2); % Brief pause to see the plot
    end
end