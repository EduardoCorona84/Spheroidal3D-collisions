function test_correction(A, b)
%In this code, we want to test how helpful the quasi-Newton corrections are to the performance of the solver. 
%We will have a single test matrix, and various levels of noisy versions of this test matrix (all made to be symmetric positive definite).
%For each of these, we will compare how useful the quasi-Newton correction is in actually helping convergence. 
%This will show if there is some benefit to using this idea of corrections, as the problem with the current test case is that the approximate matrix we have is so innaccurate, the corrections are largely useless. 
%1. High Fidelity Only, reference solution
A_func = @(x) A * x;
objective = @(x) (1/2)*x'*A_func(x) + b'*x;
opts.outer.projected = true;
opts.outer.warm_start.enabled = false;
opts.outer.warm_start.only = false;
opts.inner.enabled = false;
opts.outer.correction = false;
[~, info] = multi_fidelity_solver(A_func, A_func, b, zeros(size(b)), opts);
%find x_ref
objective_values = zeros(1, info.outer.iter);
for i = 1:info.outer.iter
    objective_values(i) = objective(info.outer.x_iters(:, i));
end
[~, min_ind] = min(objective_values);
x_ref = info.outer.x_iters(:, min_ind);
x_ref_rel_errors =vecnorm(info.outer.x_iters(:, 1:min_ind) - x_ref)/norm(x_ref);

noise_levels = [4, 3, 2, 1, 1/2, 1e-1, (1/2*1e-1), 1e-2, (1/2*1e-2), 1e-3, (1/2*1e-3)];
matrix_rel_errors = zeros(1, length(noise_levels));

%1 is alternating no correction
%2 is alternating correction sr1
%3 is alternating correction bfgs
%4 is alternating correction dfp
iterate_rel_errors = cell(length(noise_levels), 4);

%for every noise level, we want to test the effectiveness of the only high fidelity, alternating low and high fidelity, and corrected alternating low and high fidelity
for i = 1:length(noise_levels)
    matrix_noise_level = noise_levels(i);
    noise = randn(size(A));
    noise = (noise + noise')/2;
    noise = (noise/norm(noise, "fro"))*matrix_noise_level;
    A_noisy = A + noise;
    %find the spectral norm relative error, for reference
    matrix_rel_errors(i) = norm(A_noisy - A) / norm(A);
    A_noisy_func = @(x) A_noisy * x;
    %now test alternating high and low fidelity
    opts.outer.correction = false;
    opts.outer.adaptive = 'none';
    opts.inner.enabled = true;
    opts.inner.max_iter = 1; %take 1 low fidelity step between high fidelities
    [~, info] = multi_fidelity_solver(A_func, A_noisy_func, b, zeros(size(b)), opts);
    iterate_rel_errors{i, 1} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);
    %now test alternating high and low fidelity with correction
    opts.outer.correction = true;
    opts.outer.low_update = 'sr1';
    [~, info] = multi_fidelity_solver(A_func, A_noisy_func, b, zeros(size(b)), opts);
    iterate_rel_errors{i, 2} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);
    opts.outer.correction = true;
    opts.outer.low_update = 'bfgs';
    [~, info] = multi_fidelity_solver(A_func, A_noisy_func, b, zeros(size(b)), opts);
    iterate_rel_errors{i, 3} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);
    opts.outer.correction = true;
    opts.outer.low_update = 'dfp';
    [~, info] = multi_fidelity_solver(A_func, A_noisy_func, b, zeros(size(b)), opts);
    iterate_rel_errors{i, 4} = vecnorm(info.outer.x_iters - x_ref)/norm(x_ref);
end

%now plot results, different figure for each noise level
for i = 1:length(noise_levels)
    figure;
    semilogy(x_ref_rel_errors, 'k', 'DisplayName', 'Reference');
    hold on;
    semilogy(iterate_rel_errors{i, 1}, 'r', 'DisplayName', 'No Correction');
    semilogy(iterate_rel_errors{i, 2}, 'b', 'DisplayName', 'With Correction SR1');
    semilogy(iterate_rel_errors{i, 3}, 'g', 'DisplayName', 'With Correction BFGS');
    semilogy(iterate_rel_errors{i, 4}, 'm', 'DisplayName', 'With Correction DFP');
    title(['Noise Level: ', num2str(noise_levels(i)), ', Relative Error', 'Matrix Rel Error ', num2str(matrix_rel_errors(i))]);
    xlabel('High Fidelity Matvecs');
    ylabel('Relative Error');
    legend show;
    grid on;
    saveas(gcf, ['noise_level_', num2str(noise_levels(i)), 'inner_SR1_1.svg']);
    hold off;
end