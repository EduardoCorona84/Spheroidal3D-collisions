%add necessary path to source code
addpath('../src/');
%add necessary path to the utility file
addpath('../../utilities/');

%create test matrices
problem_size = 100;
matrices = construct_test_matrices(problem_size);

%create noisy versions of a test matrix, we'll just use the linear50 for this.
noise_levels = linspace(1, 20, 20); %no idea what relative error this will actually produce
noisy_matrices = cell(length(noise_levels), 2);
for i = 1:length(noise_levels)
    %create the noisy matrices for each noise level, with the associated error struct.
    [noisy_matrices{i, 1}, noisy_matrices{i, 2}] = construct_noisy_matrix(matrices.linear50.matrix, noise_levels(i));
end

%identify which matrices were successfully created
successful = cellfun(@(x) ~isempty(x), noisy_matrices(:, 1));
noisy_matrices = noisy_matrices(successful, :);

b = (2*rand(problem_size, 1) - 1);

%pass the test matrices to the method comparison function
%saved .fig files in the plots folder
method_comparisons(matrices.linear50.matrix, noisy_matrices{20,1}, b)

disp(noisy_matrices{20,2}.spectral.rel*100);