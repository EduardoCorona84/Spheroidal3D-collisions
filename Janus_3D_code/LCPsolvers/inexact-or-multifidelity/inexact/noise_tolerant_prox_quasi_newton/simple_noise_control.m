function [s_k, y_k] = simple_noise_control(x_km1, s_k, grad_km1, y_k, fcnGrad, opts)

    noisy = true;
    iter = 0;
    beta = 1;
    q = s_k;

    while noisy && (iter < opts.noise_control_max_iter)

        beta = opts.noise_control_line_search_param * beta;

        x_k = x_km1 + beta * q;
        [~, grad_k] = fcnGrad(x_k);

        s_k = x_k - x_km1;
        y_k = grad_k - grad_km1;

        noisy = y_k' * s_k < 2 * (1 + opts.noise_control_parameter) * opts.noise * norm(s_k);

        iter = iter + 1;
    end
end
