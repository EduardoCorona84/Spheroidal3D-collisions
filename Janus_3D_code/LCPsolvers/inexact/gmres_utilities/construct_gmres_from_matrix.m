function my_gmres = construct_gmres_from_matrix(matrix_inverse, matrix_norm, matrix_inverse_norm, restart, max_iter)
    %This function constructs a function that performs a matrix matvec with the gmres solve.

    my_gmres = @myGmresFunction;

    function [output, iters, rel_residual, abs_residual] = myGmresFunction(input, tol, tol_type)

        if strcmp(tol_type, 'absolute')
            tol_rel = tol / (norm(input)*matrix_inverse_norm);
        elseif strcmp(tol_type, 'relative')
            tol_rel = tol / (matrix_norm*matrix_inverse_norm);
        else
            error('Unknown tolerance type: %s', tol_type);
        end

        [output, ~, residual, iterations] = gmres(matrix_inverse, input, restart, tol_rel, max_iter);

        rel_residual = matrix_norm * matrix_inverse_norm * residual;

        abs_residual = norm(input) * matrix_norm * residual;

        iters = prod(iterations);

    end

end