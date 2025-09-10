function [s_k, y_k] = projected_noise_control(x_km1, kappa, p, h0, U, V, opts)
    noisy = true;
    iter = 0;
    while noisy && iter < opts.max_noise_control_iter
        
        kappa = opts.noise_control_line_search_param * kappa;
        
        y = x_km1 + kappa * p;
        x_k = prox(y, h0, U, V, opts);
        [~, grad_k] = fcnGrad(x_k);
        s_k = x_k - x_km1;
        y_k = grad_k - grad_km1;

        noisy = y_k' * s_k < 2 * (1 + opts.noise_control_parameter) * opts.noise * norm(s_k);

        iter = iter + 1;
    end
end