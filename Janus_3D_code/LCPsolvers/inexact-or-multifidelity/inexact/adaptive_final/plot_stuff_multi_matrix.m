function plot_stuff_multi_matrix(matrices, b)
    %Takes in multiple matrices, creates a plot of adaptive vs not adaptive for various noise levels
    %matrices should be a struct with fields like matrices.exp, matrices.well_conditioned, etc.
    
    noise_levels = [1e0, 1e-2, 1e-4, 1e-6, 1e-8];
    addpath("../gmres_utilities");
    
    % Get matrix names from the struct
    matrix_names = fieldnames(matrices);
    if length(matrix_names) ~= 4
        error('Expected exactly 4 matrices in the input struct');
    end
    
    % Create subplots for 4 matrices (2x2 grid)
    figure('Position', [100, 100, 1200, 800]);
    
    for mat_idx = 1:4
        subplot(2, 2, mat_idx);
        
        current_matrix = matrices.(matrix_names{mat_idx});
        
        % Reference solution with no noise
        opts.tol_rel = 1e-30;
        opts.tol_abs = 1e-30;
        opts.storeIts = true;
        opts.max_iter = 400;
        opts.noise_control.enabled = false;
        opts.noise_adaptive.quad = false;
        opts.noise_control.type = 'none';
        opts.r = 10;
        opts.kappa.init = 'uniform';
        opts.kappa.fwd = 'uniform';
        opts.kappa.bwd = 'uniform';
        
        objective = @(x) (1/2)*x'*current_matrix.matrix*x + b'*x;
        fcnGrad = create_gradient(current_matrix.matrix, b);
        [~, info_ref] = proxQuasiNewton(fcnGrad, zeros(size(b)), opts);
        objective_ref = zeros(info_ref.iter, 1);
        for i = 1:info_ref.iter
            objective_ref(i) = objective(info_ref.iterHist(i, :)');
        end
        [ref_val, ref_ind] = min(objective_ref);
        x_ref = info_ref.iterHist(ref_ind, :)';

        restart = 5;
        gmres_max_iter = 200;

        % Create Noisy Function Eval
        fcnGrad = construct_gmres_quadratic_gradient(current_matrix.inverse, b, current_matrix.matrixNorm, current_matrix.inverseNorm, restart, gmres_max_iter);
        error_data = cell(length(noise_levels), 1);
        iter_data = cell(length(noise_levels), 1);

        % Test fixed noise levels
        for i = 1:length(noise_levels)
            opts.noise_control.noise = noise_levels(i);
            [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
            error_cur = vecnorm(info.iterHist(1:info.iter, :)' - x_ref)/norm(x_ref);
            error_data{i} = error_cur;
            iter_data{i} = cumsum(info.gmres_iters(1:info.iter));
        end

        % Test adaptive approach
        opts.noise_control.enabled = true;
        opts.noise_control.type = 'perturbed_adaptive_no_skip';
        opts.noise_control.parameter = 1;
        opts.adjust = 1/2;
        opts.noise_control.noise = 1/2;

        [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
        error_skip_perturbed = vecnorm(info.iterHist(1:info.iter, :)' - x_ref)/norm(x_ref);
        iter_skip_perturbed = cumsum(info.gmres_iters(1:info.iter));

        % Plot for this matrix
        semilogy(iter_skip_perturbed, error_skip_perturbed, 'LineWidth', 2, 'DisplayName', 'Adaptive with No Skips')
        hold on
        for i = 1:length(noise_levels)
            semilogy(iter_data{i}, error_data{i}, 'DisplayName', ['Fixed: ' num2str(noise_levels(i))])
        end
        xlabel('GMRES Iterations');
        ylabel('Relative Objective Value Error');
        title(['Matrix: ' strrep(matrix_names{mat_idx}, '_', ' ')]);
        legend('Location', 'best');
        hold off;
    end

    sgtitle('Adaptive vs Fixed Noise Control, 1/2 - Multiple Matrices');
    saveas(gcf, 'adaptive_vs_fixed_multi_matrix_1_2.svg');
    hold off;
    
    %do the same but with 1/10 as adjust param
    for mat_idx = 1:4
        subplot(2, 2, mat_idx);
        
        current_matrix = matrices.(matrix_names{mat_idx});
        
        % Reference solution with no noise
        opts.tol_rel = 1e-30;
        opts.tol_abs = 1e-30;
        opts.storeIts = true;
        opts.max_iter = 400;
        opts.noise_control.enabled = false;
        opts.noise_adaptive.quad = false;
        opts.noise_control.type = 'none';
        opts.r = 10;
        opts.kappa.init = 'uniform';
        opts.kappa.fwd = 'uniform';
        opts.kappa.bwd = 'uniform';
        
        objective = @(x) (1/2)*x'*current_matrix.matrix*x + b'*x;
        fcnGrad = create_gradient(current_matrix.matrix, b);
        [~, info_ref] = proxQuasiNewton(fcnGrad, zeros(size(b)), opts);
        objective_ref = zeros(info_ref.iter, 1);
        for i = 1:info_ref.iter
            objective_ref(i) = objective(info_ref.iterHist(i, :)');
        end
        [ref_val, ref_ind] = min(objective_ref);
        x_ref = info_ref.iterHist(ref_ind, :)';

        restart = 5;
        gmres_max_iter = 200;

        % Create Noisy Function Eval
        fcnGrad = construct_gmres_quadratic_gradient(current_matrix.inverse, b, current_matrix.matrixNorm, current_matrix.inverseNorm, restart, gmres_max_iter);
        error_data = cell(length(noise_levels), 1);
        iter_data = cell(length(noise_levels), 1);

        % Test fixed noise levels
        for i = 1:length(noise_levels)
            opts.noise_control.noise = noise_levels(i);
            [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
            error_cur = vecnorm(info.iterHist(1:info.iter, :)' - x_ref)/norm(x_ref);
            error_data{i} = error_cur;
            iter_data{i} = cumsum(info.gmres_iters(1:info.iter));
        end

        % Test adaptive approach
        opts.noise_control.enabled = true;
        opts.noise_control.type = 'perturbed_adaptive_no_skip';
        opts.noise_control.parameter = 1;
        opts.adjust = 1/10;
        opts.noise_control.noise = 1/2;

        [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
        error_skip_perturbed = vecnorm(info.iterHist(1:info.iter, :)' - x_ref)/norm(x_ref);
        iter_skip_perturbed = cumsum(info.gmres_iters(1:info.iter));

        % Plot for this matrix
        semilogy(iter_skip_perturbed, error_skip_perturbed, 'LineWidth', 2, 'DisplayName', 'Adaptive with No Skips')
        hold on
        for i = 1:length(noise_levels)
            semilogy(iter_data{i}, error_data{i}, 'DisplayName', ['Fixed: ' num2str(noise_levels(i))])
        end
        xlabel('GMRES Iterations');
        ylabel('Relative Objective Value Error');
        title(['Matrix: ' strrep(matrix_names{mat_idx}, '_', ' ')]);
        legend('Location', 'best');
        hold off;
    end
    
    sgtitle('Adaptive vs Fixed Noise Control, 1/10 - Multiple Matrices');
    saveas(gcf, 'adaptive_vs_fixed_multi_matrix_1_10.svg');
end