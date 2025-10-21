addpath('../inexact-or-multifidelity/utilities/');
addpath('../solvers/');
n = 100;
matrices = construct_test_matrices(n);
Amat = matrices.linear50.matrix;
b = 2*(randn(n, 1) - 1/2);

fg = LOCAL_create_fg(Amat, b);
opts.storeIts = true;
opts.solver = 'apgd';
opts.b = b;
opts.A = @(x) Amat * x;
[x, info] = accelerated_prox_grad(fg, zeros(n, 1), opts);

cvx begin
    variable x_cvx(n)
    minimize( (1/2)*quad_form(x_cvx, Amat) + dot(x_cvx, b))
    subject to
        x_cvx >= 0
cvx end



function fg = LOCAL_create_fg(A_init,b)
    %this function returns a function handle fg that fits into nic's code/template given a matrix evluation A and a vector b.
    if ~isa(A_init, 'function_handle')
        A = @(x) A_init*x;
    end
    fg = @temp_fg;
    function [f_k, grad_k] = temp_fg(x_k, Ax_k)

        %If we don't have a previous evaluation, compute it
        if nargin < 2 || isempty(Ax_k) 
            Ax_k = A(x_k);
        end

        grad_k = Ax_k + b;
        f_k = 0.5*x_k'*Ax_k + b'*x_k;
    end

end