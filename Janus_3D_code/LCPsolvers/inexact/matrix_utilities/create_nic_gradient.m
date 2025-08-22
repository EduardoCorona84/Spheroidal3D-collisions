function myFunc = create_nic_gradient(matrix, b)
    myFunc = @myGrad;
    function [f_k, grad_k] = myGrad(x)
        f_k = 0.5 * x' * matrix * x + b' * x;
        grad_k = matrix * x + b;
    end

end