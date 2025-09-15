addpath('../test_data/')
fname = '4x4x4_amphi_500';
load([fname '.mat'], ...
    'A_list', 'A_diag_list', 'b_list');

A = A_list{200};
A = (A + A')/2;
Ahat = A_diag_list{200};
Ahat = (Ahat + Ahat')/2;

b = b_list{200};

A_func = @(x) A * x;
Ahat_func = @(x) Ahat * x;

test_single_matrix(A_func, Ahat_func, b, A, Ahat)

update_comparisons(A_func, Ahat_func, b, A, Ahat)

