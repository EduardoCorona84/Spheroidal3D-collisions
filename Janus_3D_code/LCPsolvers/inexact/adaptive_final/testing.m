addpath("../matrix_utilities");
problem_size = 100;
matrices = construct_test_matrices(problem_size);
b = (rand(problem_size, 1) - 1/2)*100;

% Select 4 matrices for comparison
matrix_names = fieldnames(matrices);
selected_matrices = struct();
for i = 1:min(4, length(matrix_names))
    selected_matrices.(matrix_names{i}) = matrices.(matrix_names{i});
end

plot_stuff_multi_matrix(selected_matrices, b);


