function [A_noisy, error_struct] = construct_noisy_matrix(A, noise_level)
    %This function takes in a matrix A and and returns a noisy version of A, A_noisy, such that the relative error in the spectral norm is between noise_level and 1.1*noise_level, (we use a forward tracking line search, is one way of thinking about it). The function also returns a struct error_struct which contains the absolute and relative error in both the spectral and frobenius norms.
    
    iter = 0;
    noise_matrix = 2*rand(size(A)) - 1;
    noise_matrix = (noise_matrix + noise_matrix')/2;
    noise_matrix = noise_matrix/norm(noise_matrix, 'fro');
    curr_noise = 1e-1;
    noise_matrix = curr_noise*noise_matrix;
    line_search_scaling = 1.1;
    mode = 'increasing';
    while iter < 1000
        A_noisy = A + noise_matrix;
        error_struct.spectral.abs = norm(A - A_noisy);
        error_struct.spectral.rel = error_struct.spectral.abs/norm(A);
        
        %we check if we have enough noise
        if error_struct.spectral.rel >= noise_level 

            %anytime we switch modes, we decrease the line_search_scaling
            if strcmpi(mode, 'increasing')
                line_search_scaling = (line_search_scaling + 1)/2;
            end

            mode = 'decreasing';

            if error_struct.spectral.rel <= 1.1*noise_level
                break;
            else
                %do nothing here
            end
        else
            if strcmpi(mode, 'decreasing')
                line_search_scaling = (line_search_scaling + 1)/2;
            end
            mode = 'increasing';  
        end


        %increase or decrease the noise depending on the mode
        switch mode
            case 'increasing'
                noise_matrix = noise_matrix * line_search_scaling;
            case 'decreasing'
                noise_matrix = noise_matrix / line_search_scaling;
        end


    iter = iter + 1;    

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