function [A_noisy, error_struct] = construct_noisy_matrix(A, noise_level)
    %This function takes in a matrix A and a noise level and returns a noisy version of A. Importantly, we assume that A is symmetric positive definite and that the noisy matrix is also symmetric positive definite. This function will repeatedly add noise until the noisy matrix is positive definite.
    indefinite = true;
    
    iter = 0;
    while indefinite && iter < 100
        %Randomly Construct a Matrix with entries from U(-1, 1);
        noise_matrix = 2*rand(size(A)) - 1;
        %Make it symmetric
        noise_matrix = (noise_matrix + noise_matrix')/2;
        noise_matrix = noise_level * noise_matrix/norm(noise_matrix, 'fro');
        A_noisy = A + noise_matrix;
        %ensure positive definiteness
        min_eig = min(eig(A_noisy));
        if min_eig > 0
            indefinite = false;
        end
    end

    if iter == 100
        A_noisy = [];
        error_struct = [];
        warning('Could not construct a positive definite noisy matrix after 100 iterations. Returning empty matrixs.');
        return;

    else
        error_struct.spectral.abs = norm(A - A_noisy);
        error_struct.spectral.rel = error_struct.spectral.abs/norm(A);
        error_struct.frobenius.abs = norm(A - A_noisy, 'fro');
        error_struct.frobenius.rel = error_struct.frobenius.abs/norm(A, 'fro');
    end

end