addpath('../inexact/matrix_utilities/')
matrices = construct_test_matrices(200);
A = matrices.linearSize.matrix;
b = (rand(size(A, 1), 1) - 1/2)*10;

test_correction(A, b);