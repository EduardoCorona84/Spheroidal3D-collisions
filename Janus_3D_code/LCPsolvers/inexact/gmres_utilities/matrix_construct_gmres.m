function my_gmres = matrix_construct_gmres(matrix, matrix_inverse_norm, restart, max_iter)

    my_gmres = @myGmresFunction;

    function [output, iters, residual] = myGmresFunction(input, tol_abs)

        tol_rel = tol_abs / norm(input);
        [output, ~, residual, iterations] = gmres(matrix, input, restart, tol_rel, max_iter);
        residual = matrix_inverse_norm * residual * norm(input);
        iters = prod(iterations);

    end

end