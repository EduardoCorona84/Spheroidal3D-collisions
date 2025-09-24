function LCP_sweep(noise_level)
    %this function sweeps through different test/example matrices of the LCP, creates a low fidelity version, and runs various high/multi fidelity methods to compare performance
    addpath('../../utilities/');
    addpath('../src/');
    addpath('../data/');
    addpath('../plotting/');
    try %#ok<TRYNC>
        rng('default')
    end
    rng(2);
    fname = '../data/all_data.mat';
    load([fname '.mat'], ...
        'A_list',  'b_list');


    load('../data/noisy_matrices_rel_error_1.00e-02.mat', 'A_noisy_list', 'errors_list', 'noise_level');
    total_runs = length(A_list);
    methods = {'High Fidelity Only', 'Low Fidelity Only', 'Warm Start High Fidelity', 'Alternating Low/High No Correction', 'Alternating Low/High with Correction (SR1, Full)', 'Alternating Low/High with Correction (SR1, 2)', 'Warm Alternating Low/High with Correction (SR1, Full)', 'Warm Alternating Low/High with Correction (SR1, 2)'};
    tolerance_levels = [1e0, 1e-2, 1e-4, 1e-6, 1e-8];
    errors = struct(methods);

    %create a min_iters array, with dimensions (total_runs, methods, tolerances)
    min_iters = zeros(total_runs, length(methods), length(noise_levels));


    good_runs = [];
    for runs = 101:105
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
        errors.HighFidelityOnly = vecnorm(info_high.outer.iterHist' - x_ref)/norm(x_ref);

   
        clear opts;
        %2 Low Fidelity Only/Sovling with CVX
        cvx_begin 
            variable x_warm(problem_size)
            minimize(1/2 * quad_form(x_warm, A_hats{noise_level}) + b'*x_warm)
            subject to
                x_warm >= 0
        cvx_end
        %compute error for the final low fidelity warm start iterate
        errors.LowFidelityOnly = norm(x_warm - x_ref)/norm(x_ref);

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
        errors.WarmStartHighFidelity = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);

        clear opts;
        %4. Alternating low and high fidelity with no correction
        opts.warm.enabled = false;
        opts.inner.enabled = true;
        opts.inner.max_iter = 1; %take 1 low fidelity step between high fidelities
        opts.outer.correction = false;
        opts.outer.storeIts = true;
        opts.outer.solver_opts.tol_abs = 1e-16;
        opts.outer.solver_opts.tol_rel = 1e-16;
        opts.outer.solver_opts.A = @(x) Amat*x;
        opts.outer.solver_opts.b = b;
        [~, info] = multifidelity_wrapper(fg, fg_low, zeros(problem_size, 1), opts);

        errors. = vecnorm(info.outer.iterHist' - x_ref)/norm(x_ref);

        clear opts;
        [~, info] = multi_fidelity_solver(A, Ahat, b, zeros(problem_size, 1), opts);



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


        %Now from each run find the min iteration corresponding to a specific error level
        %i corresponds to the method, j the tolerance to reach
        %for each method, check over each noise level if it met that level
        for i = 1:8
            for j = 1:length(noise_levels)
                %find the iters of method i that reached less than noise level j
                valid_iters = find(error_cell{i} < noise_levels(j));
                if isempty(valid_iters)
                    %if no iters, we assign inf to that entry of min iters (it did not happen)
                    min_iters(true_mc, i, j) = Inf;
                else
                    %if not empty, we take the min of those iters
                    min_iters(true_mc, i, j) = min(valid_iters);
                end
            end
        end

        %we check if the min iters are greater than 0. This should always happen, unless there was no collision (for checking, I should outot the size of this on the graph. We should expect 300 (or very close to it)).
        if all(min_iters(true_mc, :, :) > 0)
            mcGood = [mcGood; true_mc];
        end
        %now repeat for the next mc

    end

    method_names = {'Only High Fidelity', 'Only Low Fidelity', 'Warm Start High Fidelity', 'Alternating Low/High', 'Alternating Low/High (SR1)', 'Warm Start Alternating Low/High (SR1)', 'Warm Start Adaptive High Fidelity (SR1)', 'Warm Start Halving High Fidelity (SR1)'};
    tolerance_names = {'1e0', '1e-2', '1e-4', '1e-6', '1e-8'};

    %go through each tolerance we want to meet
    for j = 1:length(noise_levels)
        figure('Name', ['Tolerance ' tolerance_names{j}], 'NumberTitle', 'off');

        % Extract data for this tolerance level across all methods. This is a matrix with size mcGood x 8 (all the methods). Each row is a run, and column is the method for that run.
        data_for_plot = squeeze(min_iters(mcGood, :, j));
        
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
        data_for_plot(~isfinite(data_for_plot)) = NaN;
        
        % Create boxchart, which is more modern than boxplot
        b = boxchart(data_for_plot);
        
        % Set labels and colors
        ax = gca;
        ax.XTickLabel = method_names;
        

        title(['Performance at Relative Tolerance Level: ' tolerance_names{j}]);
        ylabel('Minimum GMRES Iterations');
        xlabel('Method');
        grid on;
        
        % Add success rate labels above each box
        ylims = ylim;
        for method = 1:length(method_names)
            text(method, ylims(2) * 0.95, sprintf('%.0f%%', success_percentages(method)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
        
        % Rotate x-axis labels if they overlap
        xtickangle(45);

        % Save the figure as an SVG file
        filename = ['performance_qn_tolerance_' tolerance_names{j} '.svg'];
        saveas(gcf, filename);
end