function fg = create_fg(A_init, b)
    %this function returns a function handle fg that fits into nic's code/template given a matrix evluation A and a vector b.
    if ~isa(A_init, 'function_handle')
        A = @(x) A_init*x;
    end
    fg = @temp_fg;
    function [f_k, grad_k, Ax_k] = temp_fg(x_k, Ax_km1, Aq, eta, type)

        %Allow for the function to reuse previous evaluations to create a new evaluation.
        if nargin < 2 || isempty(Ax_km1) || isempty(eta) || isempty(Aq)
            Ax_k = A(x_k);

        elseif strcmp(type, 'acceleration')
            %in this example, we have Ax_km1 = Ax_k and Aq = Ax_km1. This could be done better.
            Ax_k = (1 + eta)*Ax_km1 - eta*Aq;
        else
            Ax_k = Ax_km1 + eta*Aq;
        end

        grad_k = Ax_k + b;
        f_k = 0.5*x_k'*Ax_k + b'*x_k;
    end
end