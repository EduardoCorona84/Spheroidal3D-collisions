function results = testLCPmine()
addpath('../gmres_utilities');
addpath('../matrix_utilities/');
try %#ok<TRYNC>
    rng('default')
end
rng(2);
fname = 'all_data';
load([fname '.mat'], ...
    'A_list', 'b_list');
opts = struct( ...
    'max_iter',100, ...
    'tol_rel',1e-12, ...
    'tol_abs',1e-12, ... 
    'kappa', struct('init','uniform',...
        'fwd','uniform',...
        'bwd','uniform'),... 
    'r', 20, ...
    'qnUpdate', 'bfgs', ...
    'storeIts', true...
);
MC = length(A_list);
min_iters = zeros(400, 6, 5);
noise_levels = [1e0, 1e-2, 1e-4, 1e-6, 1e-8];

mcGood = [];
for mc = 51:450
    true_mc = mc - 50;
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
    objective = @(x) (1/2)*x'*A*x + b'*x;
    restart = 5;
    gmres_max_iter = 200;
    invA = inv(A);
    invA = (invA + invA')/2;
    fcnGrad = construct_gmres_gradient(invA , b, norm(A, 2), norm(invA, 2), restart, gmres_max_iter);
    opts.noise_control.noise = 1e-30; % No noise for reference
    [~, info_ref] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
    objective_ref = zeros(info_ref.iter, 1);
    for i = 1:info_ref.iter
        objective_ref(i) = objective(info_ref.iterHist(i, :)');
    end
    ref_val = min(objective_ref);

    %Now do fixed noise levels
    fcnGrad = construct_gmres_gradient(invA, b, norm(A, 2), norm(invA, 2), restart, gmres_max_iter);
    error_data = cell(length(noise_levels), 1);
    iter_data = cell(length(noise_levels), 1);

    for i = 1:length(noise_levels)
        opts.noise_control.noise = noise_levels(i);
        [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
        error_cur = zeros(info.iter, 1);
        for j = 1:info.iter
            error_cur(j) = abs(objective(info.iterHist(j, :)') - ref_val)/abs(ref_val);
        end
        error_data{i} = error_cur;
        iter_data{i} = cumsum(info.gmres_iters(1:info.iter));
    end

    %Now do adaptive
    opts.noise_control.enabled = true;
    opts.noise_control.type = 'perturbed_adaptive';
    opts.noise_control.parameter = 1;
    opts.noise_control.noise = 1e-1;

    [~, info] = proxQuasiNewtonAdaptive(fcnGrad, zeros(size(b)), opts);
    error_skip_perturbed = zeros(info.iter, 1);
    for i = 1:info.iter
        error_skip_perturbed(i) = abs(objective(info.iterHist(i, :)') - ref_val)/abs(ref_val);
    end
    iter_skip_perturbed = cumsum(info.gmres_iters(1:info.iter));

     %Now from each run find the min gmres iteration corresponding to a specific error level
    %i corresponds to the tolerance, j corresponds to the run
    for i = 1:length(noise_levels)
        for j = 1:length(noise_levels)
            valid_iters = iter_data{i}(error_data{i} < noise_levels(j));
            if isempty(valid_iters)
                min_iters(true_mc, i, j) = Inf; % or NaN, or some default value
            else
                min_iters(true_mc, i, j) = min(valid_iters);
            end
        end
    end
    for j = 1:length(noise_levels)
        valid_iters = iter_skip_perturbed(error_skip_perturbed < noise_levels(j));
        if isempty(valid_iters)
            min_iters(true_mc, end, j) = Inf; % or NaN, or some default value
        else
            min_iters(true_mc, end, j) = min(valid_iters);
        end
    end

    if all(min_iters(true_mc, :, :) > 0)
        mcGood = [mcGood; true_mc];
    end
    %now repeat for the next mc

end

method_names = {'1e0', '1e-2', '1e-4', '1e-6', '1e-8', 'Adaptive'};
tolerance_names = {'1e0', '1e-2', '1e-4', '1e-6', '1e-8'};

for j = 1:length(noise_levels)
    figure('Name', ['Tolerance ' tolerance_names{j}], 'NumberTitle', 'off');
    
    % Extract data for this tolerance level across all methods
    data_for_plot = squeeze(min_iters(mcGood, :, j));
    
    % Calculate success percentages for each method
    success_percentages = zeros(1, size(data_for_plot, 2));
    for method = 1:size(data_for_plot, 2)
        total_runs = size(data_for_plot, 1);
        successful_runs = sum(isfinite(data_for_plot(:, method)));
        success_percentages(method) = (successful_runs / total_runs) * 100;
    end
    
    % Replace Inf values with NaN for better visualization
    data_for_plot(~isfinite(data_for_plot)) = NaN;
    
    % Create boxplot with filled boxes
    h = boxplot(data_for_plot, 'labels', method_names);
    
    % Fill the boxes with colors
    box_handles = findobj(gca, 'Tag', 'Box');
    colors = [0.8 0 0;     % Dark red
              0 0 0.8;     % Dark blue  
              0 0.6 0;     % Dark green
              0.6 0 0.6;   % Dark magenta
              0 0.6 0.8;   % Dark cyan
              0.4 0.2 0.6]; % Dark purple
    for i = 1:length(box_handles)
        color_idx = length(box_handles) - i + 1; % Reverse order
        patch(get(box_handles(i), 'XData'), get(box_handles(i), 'YData'), ...
              colors(color_idx, :), 'FaceAlpha', 0.5);
    end

    title(['Performance at Tolerance Level: ' tolerance_names{j}]);
    ylabel('Minimum GMRES Iterations');
    xlabel('Method');
    grid on;
    
    % Add success rate labels above each box
    ylims = ylim;
    for method = 1:length(method_names)
        text(method, ylims(2) * 0.95, sprintf('%.0f%%', success_percentages(method)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', ...
            'Color', colors(method, :));
    end
    
    % Rotate x-axis labels if they overlap
    xtickangle(45);
end



%Now make a bar plot of the results.


