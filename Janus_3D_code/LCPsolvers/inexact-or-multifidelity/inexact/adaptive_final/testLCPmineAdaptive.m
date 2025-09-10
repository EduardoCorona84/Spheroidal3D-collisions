function results = testLCPmine()
addpath('../gmres_utilities');
addpath('../matrix_utilities/');
addpath('../../test_data/')
try %#ok<TRYNC>
    rng('default')
end
rng(2);
fname = '4x4x4_amphi_500';
load([fname '.mat'], ...
    'A_list', 'b_list');
opts = struct( ...
    'max_iter',100, ...
    'tol_rel',1e-30, ...
    'tol_abs',1e-30, ... 
    'kappa', struct('init','uniform',...
        'fwd','uniform',...
        'bwd','uniform'),... 
    'r', 20, ...
    'qnUpdate', 'bfgs', ...
    'storeIts', true...
);
MC = length(A_list);
min_iters = zeros(300, 6, 5);
noise_levels = [1e0, 1e-2, 1e-4, 1e-6, 1e-8];

mcGood = [];
for mc = 101:400
    true_mc = mc - 100;
    % disp(['mc = ' num2str(mc)])
    A = A_list{mc};
    b = b_list{mc};
    problem_size = size(A,2);
    x0 = zeros(problem_size,1);
    %find reference for this problem
    opts.tol_rel = 1e-30;
    opts.tol_abs = 1e-30;
    opts.max_iter = 500;
    opts.noise_control.enabled = false;
    opts.noise_adaptive.quad = false;
    opts.noise_control.type = 'none';
    opts.r = 20;
    %symmetrize the matrix
    A = (A + A')/2;
    objective = @(x) (1/2)*x'*A*x + b'*x;
    restart = 5;
    gmres_max_iter = 200;
    %invert the matrix to make the gmres solver. I don't love this.
    invA = inv(A);
    invA = (invA + invA')/2;
    fcnGrad = construct_gmres_quadratic_gradient(invA , b, norm(A, 2), norm(invA, 2), restart, gmres_max_iter);
    opts.noise_control.noise = 1e-30; % No noise for reference
    [~, info_ref] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
    objective_ref = zeros(info_ref.iter, 1);
    for i = 1:info_ref.iter
        objective_ref(i) = objective(info_ref.iterHist(i, :)');
    end
    [ref_f_val, ref_ind] = min(objective_ref);
    ref_x_val = info_ref.iterHist(ref_ind, :)';

    %Now do fixed noise levels
    fcnGrad = construct_gmres_quadratic_gradient(invA, b, norm(A, 2), norm(invA, 2), restart, gmres_max_iter);
    error_data = cell(6, 1);
    iter_data = cell(6, 1);
    %1 = 1e0 1/10 start skips
    %2 = 1e0 1/5 start skips
    %3 = 1e0 1/2 start skips
    %4 = 1e-1 start 1/10 skips
    %5 = 1e-1 start 1/5 skips
    %6 = 1e-1 start 1/2 skips

    %Now do adaptive
    opts.noise_control.enabled = true;
    opts.noise_control.type = 'perturbed_adaptive_skip';
    opts.noise_control.parameter = 1e0;
    opts.noise_control.noise = 1e0;
    opts.adjust = 1/10;

    [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
    for i = 1:info.iter
        error_data{1}(i) = norm(info.iterHist(i, :)' - ref_x_val)/norm(ref_x_val);
    end
    iter_data{1} = cumsum(info.gmres_iters(1:info.iter));

    opts.noise_control.noise = 1e0;
    opts.adjust = 1/5;

    [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
    for i = 1:info.iter
        error_data{2}(i) = norm(info.iterHist(i, :)' - ref_x_val)/norm(ref_x_val);
    end
    iter_data{2} = cumsum(info.gmres_iters(1:info.iter));

    opts.noise_control.noise = 1e0;
    opts.adjust = 1/2;

    [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
    for i = 1:info.iter
        error_data{3}(i) = norm(info.iterHist(i, :)' - ref_x_val)/norm(ref_x_val);
    end
    iter_data{3} = cumsum(info.gmres_iters(1:info.iter));

    opts.noise_control.noise = 1e-1;
    opts.adjust = 1/10;

    [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
    for i = 1:info.iter
        error_data{4}(i) = norm(info.iterHist(i, :)' - ref_x_val)/norm(ref_x_val);
    end
    iter_data{4} = cumsum(info.gmres_iters(1:info.iter));

    opts.noise_control.noise = 1e-1;
    opts.adjust = 1/5;

    [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
    for i = 1:info.iter
        error_data{5}(i) = norm(info.iterHist(i, :)' - ref_x_val)/norm(ref_x_val);
    end
    iter_data{5} = cumsum(info.gmres_iters(1:info.iter));

    opts.noise_control.noise = 1e-1;
    opts.adjust = 1/2;

    [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
    for i = 1:info.iter
        error_data{6}(i) = norm(info.iterHist(i, :)' - ref_x_val)/norm(ref_x_val);
    end
    iter_data{6} = cumsum(info.gmres_iters(1:info.iter));


    %Now from each run find the min gmres iteration corresponding to a specific error level
    %i corresponds to the noise level of the method, j the tolerance to reach
    %for each method, check over each noise level if it met that level
    for i = 1:6
        for j = 1:length(noise_levels)
            %find the iters of method i that reached less than noise level j
            valid_iters = iter_data{i}(error_data{i} < noise_levels(j));
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

method_names = {'1e0, 1/10', '1e0, 1/5', '1e0, 1/2 ', '1e-1, 1/10 ', '1e-1, 1/5 ', '1e-1, 1/2 '};
tolerance_names = {'1e0', '1e-2', '1e-4', '1e-6', '1e-8'};

%go through each tolerance we want to meet
for j = 1:length(noise_levels)
    figure('Name', ['Tolerance ' tolerance_names{j}], 'NumberTitle', 'off');
    
    % Extract data for this tolerance level across all methods. This is a matrix with size mcGood x 6 (all the methods). Each row is a run, and column is the method for that run.
    data_for_plot = squeeze(min_iters(mcGood, :, j));
    disp(size(data_for_plot));
    
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
    filename = ['performance_adaptives_' tolerance_names{j} '.svg'];
    saveas(gcf, filename);
end





