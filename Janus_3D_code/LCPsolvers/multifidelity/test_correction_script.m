addpath('../test_data/')
fname = '4x4x4_amphi_500';
load([fname '.mat'], ...
    'A_list', 'A_diag_list', 'b_list');

A = A_list{150};
A = (A + A')/2;

b = b_list{150};

test_correction(A, b);