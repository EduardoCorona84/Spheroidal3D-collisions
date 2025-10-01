function construct_noisy_matrices_from_data(data, noise_level)
    %this function takes in a noise level, loads in data from a .mat file, and creates a new .mat file with noisy versions of the matrices.
    %data is a filepath to a .mat file with A_list 
    %noise level is the level we are targetting for the relative error (spectral norm) of the matrix.
    addpath('../../utilities');
    load(data, 'A_list');
    A_noisy_list = cell(size(A_list));
    errors_list = cell(size(A_list));
    for i = 1:length(A_list)
        A = A_list{i};
        [A_noisy, errors] = construct_noisy_matrix(A, noise_level);
        A_noisy_list{i} = A_noisy;
        errors_list{i} = errors;
        noise_target_list{i} = noise_level;
    end
    save_filename = sprintf('../data/noisy_matrices_rel_error_%1.2e.mat', noise_level);
    save(save_filename, 'A_noisy_list', 'errors_list', 'noise_target_list');

end