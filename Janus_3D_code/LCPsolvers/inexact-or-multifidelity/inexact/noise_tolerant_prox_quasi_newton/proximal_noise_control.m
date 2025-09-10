function [s_k, y_k] = proximal_noise_control(x_km1, kappa, p, grad_km1, h0, U, V, fcnGrad, opts)
    noisy = true;
    iter = 0;
    while noisy && iter < opts.noise_control_max_iter
        
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

function xstar = prox(y, h0, U, V, opts)
if isempty(U) && isempty(V)
    % Project with respect to the identity
    xstar = max(0, y);
elseif size(U,2) + size(V,2) == 1
    % The sign on sigma is counter intuitive, but remember B = B0 + UU' - VV'
    if ~isempty(U)
        sigma = -1;
        w = U; 
    else 
        sigma = 1;
        w = V;
    end
    xstar = prox_rank1(y, h0, w, sigma, opts);
    return  
else
    xstar = prox_rankr(y, h0, U, V, opts);
end
end % prox