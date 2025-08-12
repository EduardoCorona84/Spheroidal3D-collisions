function [x, info] = warm_start(Ahat, b, x0, info)
    %This function simply evaluates the low fidelity operator for a fixed number of steps and then uses this as an initial guess for the high fidelity solver. TODO: Add support for a tolerance.

    for i = 1:opts.warm_start.steps 
        grad_0 = Ahat(x0) + b;
        if i == 1
            step = 1;
        else
            s_1 = x0 - x_m1;
            y_1 = grad_0 - grad_m1;
            step = (s_1'*s_1)/(s_1'*y_1);
        end
        x_1 = x0 - step * grad_0;
        %update variables
        x_m1 = x0;
        x0 = x_1;
        grad_m1 = grad_0;
    end
    x = x0;

end