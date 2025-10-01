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


    load('../data/noisy_matrices_rel_error_1.00e-01.mat', 'A_noisy_list', 'errors_list', 'noise_target_list');
    total_runs = length(A_list);
    methods = {'High Fidelity Only', 'Alternating Low/High No Correction', 'Alternating Low/High with Low Fidelity Step', 'Alternating Low/High with High Reuse', 'Alternating Low/High with Correction (SR1, 10, secant)', 'Alternating Low/High with Correction (BFGS, 10, secant)', 'Alternating Low/High with Correction (DFP, 10, secant)'};
    tolerance_levels = [1e0, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12];
    errors = cell(length(methods), 1);

    %create a min_iters array, with dimensions (total_runs, methods, tolerances)
    min_iters = zeros(total_runs, length(methods), length(tolerance_levels));
    sizes = zeros(total_runs, 1);


    good_runs = [];
    for runs = 1:5
        fprintf('Starting Run %d of %d\n', runs, total_runs);
        %get the matrix and vector
        Amat = A_list{runs};
        Amat = (Amat + Amat')/2; %symmetrize the matrix
        %now get noisy A
        Ahat_mat = A_noisy_list{runs};
        Ahat_mat = (Ahat_mat + Ahat_mat')/2; %symmetrize the matrix

        b = b_list{runs};
        sizes(runs) = length(b);

        fg = create_fg(Amat, b);
        fg_low = create_fg(Ahat_mat, b);
        problem_size = length(b);

        %1. High Fidelity Only
        opts.warm.enabled = false;
        opts.outer.max_iter = 50;
        opts.inner.enabled = false;
        opts.outer.correction = false;
        opts.outer.storeIts = true;
        opts.outer.storeFuncs = true;
        opts.outer.adaptive = 'none';
        opts.outer.solver_opts.tol_abs = 1e-30; 
        opts.outer.solver_opts.tol_rel = 1e-30;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        x0 = zeros(problem_size,1);

        [~, info_high] = multifidelity_wrapper(fg, fg, x0, opts);
        %find x_ref
        [f_ref , f_ref_ind] = min(info_high.outer.funcHist);
        x_ref = info_high.outer.iterHist(f_ref_ind, :)';

        %create error history for high, the construction has it so that every iteration is one matvec
        % Find exact match, ignoring case
        index = find(strcmpi(methods, 'High Fidelity Only'));
        errors{index} = vecnorm(info_high.outer.iterHist' - x_ref)/norm(x_ref);
        clear opts;

        %2. Alternating low and high fidelity with no correction
        opts.warm.enabled = false;
        opts.inner.enabled = true;
        opts.outer.max_iter = 50;
        opts.inner.solver_opts.max_iter = 1;
        opts.inner.gradient_mode = 'full';
        %no correction
        opts.outer.correction = false;
        opts.outer.storeIts = true;
        opts.outer.adaptive = 'none';
        opts.outer.solver_opts.tol_abs = 1e-30;
        opts.outer.solver_opts.tol_rel = 1e-30;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) Ahat_mat*x;
        opts.inner.solver_opts.b = b;
        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);

        index = find(strcmpi(methods, 'Alternating Low/High No Correction'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);
        clear opts;



        %3. Alternating low and high fidelity with low fidelity step
        opts.warm.enabled = false;
        opts.inner.enabled = true;
        opts.outer.max_iter = 50;
        opts.inner.solver_opts.max_iter = 1;
        opts.inner.gradient_mode = 'step';
        %no correction
        opts.outer.correction = false;
        opts.outer.storeIts = true;
        opts.outer.adaptive = 'none';
        opts.outer.solver_opts.tol_abs = 1e-30;
        opts.outer.solver_opts.tol_rel = 1e-30;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) Ahat_mat*x;
        opts.inner.solver_opts.b = b;
        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);

        index = find(strcmpi(methods, 'Alternating Low/High with Low Fidelity Step'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);
        clear opts;

        %4. Alternating Low/High with High Reuse
        opts.warm.enabled = false;
        opts.inner.enabled = true;
        opts.outer.max_iter = 50;
        opts.inner.solver_opts.max_iter = 1;
        opts.inner.gradient_mode = 'reuse';
        %no correction
        opts.outer.correction = false;
        opts.outer.adaptive = 'none';
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-30;
        opts.outer.solver_opts.tol_rel = 1e-30;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) Ahat_mat*x;
        opts.inner.solver_opts.b = b;
        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);

        index = find(strcmpi(methods, 'Alternating Low/High with High Reuse'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);
        clear opts;

        %5. Alternating Low/High with Correction (SR1, 10)
        opts.warm.enabled = false;
        opts.inner.enabled = true;
        opts.outer.max_iter = 50;
        opts.inner.gradient_mode = 'step';
        opts.outer.correction = true;
        opts.outer.adaptive = 'none';
        opts.outer.correction_opts.update = 'SR1';
        opts.outer.correction_opts.memory = 'compact';
        opts.outer.correction_opts.r = 10;
        opts.outer.correction_opts.direction = 'secant';
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-30;
        opts.outer.solver_opts.tol_rel = 1e-30;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) Ahat_mat*x;
        opts.inner.solver_opts.b = b;
        opts.inner.solver_opts.max_iter = 1;
        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);

        index = find(strcmpi(methods, 'Alternating Low/High with Correction (SR1, 10, secant)'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);
        clear opts;

        
        %6. Alternating Low/High with Correction (BFGS, 10)
        opts.warm.enabled = false;
        opts.inner.enabled = true;
        opts.outer.max_iter = 50;
        opts.outer.warm_correction = false;
        opts.inner.gradient_mode = 'step';
        opts.outer.correction = true;
        opts.outer.adaptive = 'none';
        opts.outer.correction_opts.update = 'bfgs';
        opts.outer.correction_opts.memory = 'compact';
        opts.outer.correction_opts.r = 10;
        opts.outer.correction_opts.direction = 'secant';
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-30;
        opts.outer.solver_opts.tol_rel = 1e-30;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) Ahat_mat*x;
        opts.inner.solver_opts.b = b;
        opts.inner.solver_opts.max_iter = 1;
        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);


        index = find(strcmpi(methods, 'Alternating Low/High with Correction (BFGS, 10, secant)'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);
        clear opts;

        %7. Alternating Low/High with Correction (DFP, 10)
        opts.warm.enabled = false;
        opts.inner.enabled = true;
        opts.outer.max_iter = 50;
        opts.inner.gradient_mode = 'step';
        opts.outer.correction = true;
        opts.outer.adaptive = 'none';
        opts.outer.correction_opts.update = 'dfp';
        opts.outer.correction_opts.memory = 'compact';
        opts.outer.correction_opts.r = 10;
        opts.outer.correction_opts.direction = 'secant';
        opts.outer.storeFuncs = true;
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-30;
        opts.outer.solver_opts.tol_rel = 1e-30;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        opts.inner.solver_opts.A = @(x) Ahat_mat*x;
        opts.inner.solver_opts.b = b;
        opts.inner.solver_opts.max_iter = 1;
        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);

        index = find(strcmpi(methods, 'Alternating Low/High with Correction (DFP, 10, secant)'));
        errors{index} = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);
        clear opts;


        %Now from each run find the min iteration corresponding to a specific error level
        %i corresponds to the method, j the tolerance to reach
        %for each method, check over each noise level if it met that level
        for i = 1:length(methods)
            for j = 1:length(tolerance_levels)
                %find the iters of method i that reached less than tolerance level j
                valid_iters = find(errors{i} < tolerance_levels(j));
                if isempty(valid_iters)
                    %if no iters, we assign inf to that entry of min iters (it did not happen)
                    min_iters(runs, i, j) = Inf;
                else
                    %if not empty, we take the min of those iters
                    min_iters(runs, i, j) = min(valid_iters);
                end
            end
        end

        %we check if the min iters are greater than 0. This should always happen, unless there was no collision (for checking, I should outot the size of this on the graph. We should expect 300 (or very close to it)).
        if all(min_iters(runs, :, :) > 0)
            good_runs = [good_runs runs];
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
            total_good_runs = size(data_for_plot, 1);
            %data_for_plot(:, method) is a vector, corresponding to the min_iter of each run for that method. A run is successful if the min_iter is not infinity, so  this sums all the sucessful runs.
            successful_runs = sum(isfinite(data_for_plot(:, method)));
            success_percentages(method) = (successful_runs / total_good_runs) * 100;
        end
        
        % copilot did this, I don't entirely get it. We should just remove these data points?
        %This just makes the infs NaNs, which I don't think does anything
        data_for_plot(~isfinite(data_for_plot)) = NaN;
        
        % Create boxchart, which is more modern than boxplot
        b = boxchart(data_for_plot);
        
        % Set labels and colors
        ax = gca;
        ax.XTickLabel = methods;
        

        title(['Performance at Relative Tolerance Level: ' num2str(tolerance_levels(j))]);
        subtitle(sprintf('Low Fidelity Matrix Relative Error: %.2e', noise_target_list{1}));
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
        base_filename = sprintf('./plots/tolerance_%s_performance1.00e-01', num2str(tolerance_levels(j)));
        saveas(gcf, [base_filename '.svg']);
        saveas(gcf, [base_filename '.fig']);

    end

end