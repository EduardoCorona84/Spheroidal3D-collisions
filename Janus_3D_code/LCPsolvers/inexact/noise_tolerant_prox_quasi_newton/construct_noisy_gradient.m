function noisy_gradient = construct_noisy_gradient(A_noisy, b)
    % Constructs a noisy gradient function based on the provided noisy matrix A_noisy.
    % The function returns a handle to compute the gradient at any point x.
    

    noisy_gradient = @n_gradient;
    function [function_value, gradient_value] = n_gradient(x)

        noisey_eval = A_noisy(x);
        gradient_value = noisey_eval + b;
        function_value = 0.5 * (x' * noisey_eval) + b' * x;

    end
end