function my_gmres = construct_gmres_gradient(matrix_inverse, b, matrix_norm, matrix_inverse_norm, restart, max_iter)

    my_gmres = @myGmresFunction;

    function [objective, gradient, iters, rel_residual, abs_residual] = myGmresFunction(input, tol, tol_type)

        if strcmp(tol_type, 'absolute')
            tol_rel = tol / (norm(input) * matrix_inverse_norm);
        elseif strcmp(tol_type, 'relative')
            tol_rel = tol / (matrix_norm * matrix_inverse_norm);
        else
            error('Unknown tolerance type: %s', tol_type);
        end
        [eval, ~, residual, iterations] = gmres(matrix_inverse, input, restart, tol_rel, max_iter);
        gradient = eval + b;
        objective = 0.5 * input' * eval + b' * input; % Quadratic objective

        rel_residual = matrix_norm * matrix_inverse_norm * residual;

        abs_residual = norm(input) * matrix_norm * residual;

        iters = prod(iterations);

    end

end