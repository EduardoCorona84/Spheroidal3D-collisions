addpath('../inexact/matrix_utilities/')
matrices = construct_test_matrices(400);
A = matrices.exp.matrix;
b = (rand(size(A, 1), 1) - 1/2)*10;

test_correction(A, b);