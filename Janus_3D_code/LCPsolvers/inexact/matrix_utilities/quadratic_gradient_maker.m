function quadratic_gradient_function = quadratic_gradient_maker(A, b, noise_level)

    if class(A) ~= "function_handle"
        A_func = @(x) A * x; % Ensure A is a function handle
    else
        A_func = A;
    end

    quadratic_gradient_function = @quad_grad_f;

    function [objective, gradient] = quad_grad_f(input)

        %First we generate noise
        noise = sqrt(2)*randn(size(input));
        noise = noise + 1/2; 
        signs = sign(randn(size(input))); %random signed noise
        noise = noise .* signs; %apply the signs to the noise
        radius_noise = (rand(1))^(1/size(input, 1)) * noise_level; %uniformly distributed over the sphere of radius noise_level
        noise = noise_level * radius_noise * (noise / norm(noise)); %scale the noise to the correct radius
        A_eval = A_func(input) + noise; %add the noise to the function evaluation
        objective = 0.5 * input' * A_eval + b' * input; %quadratic objective
        gradient = A_eval + b; %gradient of the quadratic function

    end
end