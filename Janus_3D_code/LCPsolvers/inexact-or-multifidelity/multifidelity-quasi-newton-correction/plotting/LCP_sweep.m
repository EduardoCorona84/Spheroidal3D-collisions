function LCP_sweep()
    %this function sweeps through different test/example matrices of the LCP, creates a low fidelity version, and runs various high/multi fidelity methods to compare performance
    addpath('../../utilities/');
    addpath('../src/');
    addpath('../data/');
    addpath('../plotting/');
    try %#ok<TRYNC>
        rng('default')
    end
    rng(2);
    fname = '../data/all_data';
    load([fname '.mat'], ...
        'A_list',  'b_list');


    load('../data/noisy_matrices_rel_error_5.00e-02.mat', 'A_noisy_list', 'errors_list', 'noise_level');
    total_runs = length(A_list);
    methods = {'High Fidelity Only', 'Low Fidelity Only', 'Warm Start High Fidelity', 'Alternating Low/High No Correction', 'Alternating Low/High with Correction (SR1, Full)', 'Alternating Low/High with Correction (SR1, 2)', 'Warm Alternating Low/High with Correction (SR1, Full)', 'Warm Alternating Low/High with Correction (SR1, 2)'};
    tolerance_levels = [1e0, 1e-2, 1e-4, 1e-6, 1e-8];
    errors = cell(8, 1);

    %create a min_iters array, with dimensions (total_runs, methods, tolerances)
    min_iters = zeros(total_runs, length(methods), length(tolerance_levels));


    good_runs = [];
    for runs = 1:total_runs
        %get the matrix and vector
        Amat = A_list{runs};
        Amat = (Amat + Amat')/2; %symmetrize the matrix
        %now get noisy A
        Ahat_mat = A_noisy_list{runs};
        Ahat_mat = (Ahat_mat + Ahat_mat')/2; %symmetrize the matrix

        b = b_list{runs};

        fg = create_fg(Amat, b);
        fg_low = create_fg(Ahat_mat, b);
        problem_size = length(b);

        %1. High Fidelity Only
        opts.warm.enabled = false;
        opts.inner.enabled = false;
        opts.outer.correction = false;
        opts.outer.storeIts = true;
        opts.outer.storeFuncs = true;
        opts.outer.solver_opts.tol_abs = 1e-16;
        opts.outer.solver_opts.tol_rel = 1e-16;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        x0 = zeros(problem_size,1);

        [~, info_high] = multifidelity_wrapper(fg, fg, x0, opts);
        %find x_ref
        [~ , f_ref_ind] = min(info_high.outer.funcHist);
        x_ref = info_high.outer.iterHist(f_ref_ind, :)';

        %create error history for high, the construction has it so that every iteration is one matvec
        % Find exact match, ignoring case
        index = find(strcmpi(methods, 'High Fidelity Only'));
        errors{index} = vecnorm(info_high.outer.iterHist' - x_ref)/norm(x_ref);
        disp(cond(Amat));
        disp(eig(Amat));
        disp(eig(Ahat_mat));

        clear opts;
        %2 Low Fidelity Only/Sovling with CVX
        cvx_begin 
            variable x_warm(problem_size)
            minimize(1/2 * quad_form(x_warm, Ahat_mat) + b'*x_warm)
            subject to
                x_warm >= 0
        cvx_end
        %compute error for the final low fidelity warm start iterate
        index = find(strcmpi(methods, 'Low Fidelity Only'));
        errors{index} = norm(x_warm - x_ref)/norm(x_ref);

        clear opts;
        %3. Warm Start into High Fidelity Only, can just do this by passing x_warm or using the wrapper to call a warm start solver. We will just pass x_warm to save cost.
        opts.warm.enabled = false;
        opts.inner.enabled = false;
        opts.outer.correction = false;
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-16;
        opts.outer.solver_opts.tol_rel = 1e-16;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        [~, info] = multifidelity_wrapper(fg, fg, x_warm, opts);
        %compute error for all the high fidelity iterates
        index = find(strcmpi(methods, 'Warm Start High Fidelity'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);

        clear opts;
        %4. Alternating low and high fidelity with no correction
        opts.warm.enabled = false;
        opts.inner.enabled = true;
        opts.inner.max_iter = 1; %take 1 low fidelity step between high fidelities
        %no correction
        opts.outer.correction = false;
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-16;
        opts.outer.solver_opts.tol_rel = 1e-16;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) A_hat_mat*x;
        opts.inner.solver_opts.b = b;
        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);

        index = find(strcmpi(methods, 'Alternating Low/High No Correction'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);

        clear opts;
        %5 Alternating low and high fidelity with correction (SR1, full)
        opts.warm.enabled = false;
        opts.inner.enabled = true;
        opts.inner.max_iter = 1; %take 1 low fidelity step between high fidel
        opts.outer.correction = true;
        opts.outer.correction_opts.update = 'SR1';
        opts.outer.correction_opts.memory = 'dense full';
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-16;
        opts.outer.solver_opts.tol_rel = 1e-16;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) A_hat_mat*x;
        opts.inner.solver_opts.b = b;

        index = find(strcmpi(methods, 'Alternating Low/High with Correction (SR1, Full)'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);

        clear opts;
        %6. Alternating low and high fidelity with (SR1, 2)
        opts.warm.enabled = false;
        opts.inner.enabled = true;
        opts.inner.max_iter = 1; %take 1 low fidelity step between high fidel
        opts.outer.correction = true;
        opts.outer.correction_opts.update = 'SR1';
        opts.outer.correction_opts.memory = 'compact';
        opts.outer.correction_opts.direction = 'matvec';
        opts.outer.correction_opts.persistent_first_update = true;
        opts.outer.correction_opts.r = 2;
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-16;
        opts.outer.solver_opts.tol_rel = 1e-16;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) A_hat_mat*x;
        opts.inner.solver_opts.b = b;

        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);

        index = find(strcmpi(methods, 'Alternating Low/High with Correction (SR1, 2)'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);

        clear opts;
        %7. Warm Alternating Low/High with Correction (SR1, Full)
        opts.warm.enabled = true;
        opts.outer.warm_correction = true;
        opts.outer.b_correction = true;
        opts.inner.enabled = true;
        opts.inner.max_iter = 1; %take 1 low fidelity step between high fidel
        opts.outer.correction = true;
        opts.outer.correction_opts.update = 'SR1';
        opts.outer.correction_opts.memory = 'dense full';
        opts.outer.correction_opts.direction = 'secant';
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-16;
        opts.outer.solver_opts.tol_rel = 1e-16;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) A_hat_mat*x;
        opts.inner.solver_opts.b = b;

        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);

        index = find(strcmpi(methods, 'Warm Alternating Low/High with Correction (SR1, Full)'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);

        clear opts
        %8 Warm Alternating Low/High with Correction (SR1, 2)
        opts.warm.enabled = true;
        opts.outer.warm_correction = true;
        opts.outer.b_correction = true;
        opts.inner.enabled = true;
        opts.inner.max_iter = 1; %take 1 low fidelity step between high fidel
        opts.outer.correction = true;
        opts.outer.correction_opts.update = 'SR1';
        opts.outer.correction_opts.memory = 'compact';
        opts.outer.correction_opts.direction = 'matvec';
        opts.outer.correction_opts.persistent_first_update = true;
        opts.outer.correction_opts.r = 2;
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-16;
        opts.outer.solver_opts.tol_rel = 1e-16;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) A_hat_mat*x;
        opts.inner.solver_opts.b = b;

        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);

        index = find(strcmpi(methods, 'Warm Alternating Low/High with Correction (SR1, 2)'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);

        %Now from each run find the min iteration corresponding to a specific error level
        %i corresponds to the method, j the tolerance to reach
        %for each method, check over each noise level if it met that level
        for i = 1:length(methods)
            for j = 1:length(tolerance_levels)
                %find the iters of method i that reached less than tolerance level j
                valid_iters = find(errors{i} < tolerance_levels(j));
                if isempty(valid_iters)
                    %if no iters, we assign inf to that entry of min iters (it did not happen)
                    min_iters(run, i, j) = Inf;
                else
                    %if not empty, we take the min of those iters
                    min_iters(run, i, j) = min(valid_iters);
                end
            end
        end

        %we check if the min iters are greater than 0. This should always happen, unless there was no collision (for checking, I should outot the size of this on the graph. We should expect 300 (or very close to it)).
        if all(min_iters(run, :, :) > 0)
            good_runs = [good_runs run];
        end
        %now repeat for the next mc

    end


    %go through each tolerance we want to meet
    for j = 1:length(tolerance_levels)
        figure('Name', sprintf('Tolerance %.0e', tolerance_levels(j)), 'NumberTitle', 'off');

        % Extract data for this tolerance level across all methods. This is a matrix with size mcGood x 8 (all the methods). Each row is a run, and column is the method for that run.
        data_for_plot = squeeze(min_iters(good_runs, :, j));

        % Calculate success percentages for each method.
        %Create a vector with the columns corresponding to each method.
        success_percentages = zeros(1, size(data_for_plot, 2));
        for method = 1:size(data_for_plot, 2)
            total_runs = size(data_for_plot, 1);
            %data_for_plot(:, method) is a vector, corresponding to the min_iter of each run for that method. A run is successful if the min_iter is not infinity, so  this sums all the sucessful runs.
            successful_runs = sum(isfinite(data_for_plot(:, method)));
            success_percentages(method) = (successful_runs / total_runs) * 100;
        end
        
        % copilot did this, I don't entirely get it. We should just remove these data points?
        %This just makes the infs NaNs, which I don't think does anything
        data_for_plot(~isfinite(data_for_plot)) = NaN;
        
        % Create boxchart, which is more modern than boxplot
        b = boxchart(data_for_plot);
        
        % Set labels and colors
        ax = gca;
        ax.XTickLabel = methods;
        

        title(['Performance at Relative Tolerance Level: ' tolerance_levels{j}]);
        subtitle(sprintf('Low Fidelity Matrix Relative Error: %.2e', errors_list{run}.spectral.rel));
        ylabel('Minimum High Fidelity Matvecs');
        xlabel('Method');
        grid on;
        
        % Add success rate labels above each box
        ylims = ylim;
        for method = 1:length(methods)
            text(method, ylims(2) * 0.95, sprintf('%.0f%%', success_percentages(method)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
        
        % Rotate x-axis labels if they overlap
        xtickangle(45);

        % Save the figure as an SVG and fig
        base_filename = sprintf('/plots/tolerance_%s_performance', tolerance_levels{j});
        saveas(gcf, [base_filename '.svg']);
        saveas(gcf, [base_filename '.fig']);

    end

end