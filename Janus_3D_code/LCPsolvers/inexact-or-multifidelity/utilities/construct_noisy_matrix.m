function [A_noisy, error_struct] = construct_noisy_matrix(A, noise_level)
    %This function takes in a matrix A and and returns a noisy version of A, A_noisy, such that the relative error in the spectral norm is between greater than noise_level, (we use a forward tracking line search, is one way of thinking about it). The function also returns a struct error_struct which contains the absolute and relative error in both the spectral and frobenius norms.
    
    iter = 0;
    noise_matrix = 2*rand(size(A)) - 1;
    noise_matrix = (noise_matrix + noise_matrix')/2;
    noise_matrix = noise_matrix/norm(noise_matrix, 'fro');
    curr_noise = 1e-1;

    while iter < 1000
        noise_matrix = noise_matrix*curr_noise;
        A_noisy = A + noise_matrix;
        error_struct.spectral.abs = norm(A - A_noisy);
        error_struct.spectral.rel = error_struct.spectral.abs/norm(A);
        
        %we check if we have enough noise
        if error_struct.spectral.rel >= noise_level
            break;
        end

        %if not, we continue, but now with more noise
        curr_noise = curr_noise * 1.1;
    end

    if iter == 1000
        A_noisy = [];
        error_struct = [];
        warning('Could not construct a positive definite noisy matrix after 1000 iterations. Returning empty matrixs.');
        return;

    else
        error_struct.frobenius.abs = norm(A - A_noisy, 'fro');
        error_struct.frobenius.rel = error_struct.frobenius.abs/norm(A, 'fro');
    end

end