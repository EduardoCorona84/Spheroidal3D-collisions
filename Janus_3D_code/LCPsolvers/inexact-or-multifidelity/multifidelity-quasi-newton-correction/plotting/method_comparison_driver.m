%add necessary path to source code
addpath('../src/');
%add necessary path to the utility file
addpath('../../utilities/');

%create test matrices
problem_size = 100;
matrices = construct_test_matrices(problem_size);

%create noisy versions of a test matrix, we'll just use the linear50 for this.
noise_levels = linspace(1, 40, 40); %no idea what relative error this will actually produce
fname = '../data/all_data';
load([fname '.mat'], 'A_list',  'b_list');
A = A_list{1};
A = (A + A')/2; %make sure symmetric
b = b_list{1};

[noisy_matrix, error_data] = construct_noisy_matrix(A, 0.01);

%pass the test matrices to the method comparison function
%saved .fig files in the plots folder
method_comparisons(A, noisy_matrix, b)

disp(error_data.spectral.rel*100);