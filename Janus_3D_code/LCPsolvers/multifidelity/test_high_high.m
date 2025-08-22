%in this file, we want to test how the multifidelity does with the low fidelity operator being the high fidelity operator (the correction matrices should all be 0). We also want to systematically compare the results of the iterates of the high only operator to the results of nic's code.
addpath('../test_data/')
fname = '4x4x4_amphi_500';
load([fname '.mat'], ...
    'A_list', 'A_diag_list', 'b_list');

A = A_list{200};
A = (A + A')/2;

b = b_list{200};

A_func = @(x) A * x;
Ahat_func = A_func;

%4. Alternating low and high fidelity with no correction
%disable adaptive high fidelity switching
opts.outer.adaptive = 'none';
opts.inner.enabled = false;
opts.inner.max_iter = 1; %take 1 low fidelity step between high fidelities


%5. Alternating low and high fidelity with correction (SR1)
opts.outer.max_iter = 10;
opts.outer.correction = true;
opts.outer.store_updates = true;
opts.outer.low_update = 'bfgs';

[~, info] = multi_fidelity_solver(A_func, A_func, b, zeros(size(b)), opts);
matrix_errors = zeros(1, opts.outer.max_iter);
for i = 2:opts.outer.max_iter
    matrix_errors(i) = norm(info.outer.updates{i}, 'fro');
end

disp(size(matrix_errors));
disp(matrix_errors);
plot(matrix_errors);
