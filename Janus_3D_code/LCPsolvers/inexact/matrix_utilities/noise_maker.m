function noisy_function = noise_maker(f)
    %this function returns a function evaluation with noise that is uniformly distributed over the sphere of radius noise_level centered at the true function evaluation.
    noisy_function = @noisy_f;
    function output = noisy_f(input, noise_level)
        %First we generate noise
        noise = sqrt(2)*randn(size(input));
        noise = noise + 1/2; 
        signs = sign(randn(size(input))); %random signed noise
        noise = noise .* signs; %apply the signs to the noise
        radius_noise = (rand(1))^(1/size(input, 1)) * noise_level; %uniformly distributed over the sphere of radius noise_level
        noise = noise_level * radius_noise * (noise / norm(noise)); %scale the noise to the correct radius
        output = f(input) + noise; %add the noise to the function evaluation
    end
end