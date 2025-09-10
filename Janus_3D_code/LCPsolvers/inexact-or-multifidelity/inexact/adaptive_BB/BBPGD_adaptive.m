function [x, x_iters, gmres_iters, kkt_errors, mu_estimates, L_estimates, iters] = BBPGD_adaptive(A, b, x0, max_iter, opt_params)
% May 2018 Wen Yan
%This does not resemble the original version really at all but it feels weird getting rid of attribution.

%set parameters
adaptive = opt_params.adaptive;
opt_tol = opt_params.opt_tol;
noise_tol = opt_params.noise;
noise_control_param = opt_params.noise_control_param;
noise_shrink_param = opt_params.noise_shrink_param;
projection = opt_params.projection;


%initialize variables
x_iters = zeros(length(x0), max_iter);
gmres_iters = zeros(1, max_iter);
kkt_errors = zeros(1, max_iter);
mu_estimates = zeros(1, max_iter);
L_estimates = zeros(1, max_iter);
mu0 = inf; 
L0 = -inf;

%Compute first gradient
iter = 1;
[eval, curr_gmres_iters, ~] = A(x0, noise_tol);
grad0 = eval + b;
stepsize = 1;
gmres_iters(iter) = curr_gmres_iters; % count the first GMRES iteration

while (iter <= max_iter )

    %take step
    x1 = x0 - stepsize*grad0; 
    if projection 
        x1(x1 < 0) = 0;
    end
    x_iters(: , iter) = x1;

    %Convergences check
    
    %KKT condition
    phi = min(x1, grad0);
    kkt_errors(iter) = 0.5 * (phi' * phi);

    %relative error
    if norm(x1 - x0)/norm(x0) < opt_tol
        break;
    end
    %KKT condition


    %Compute the gradient at the new point
    [eval, curr_gmres_iters, ~] = A(x1, noise_tol);
    grad1 = eval + b;
    %Iterations are defined by the gradient evaluation, so store in the next entry of the array.
    gmres_iters(iter + 1) = gmres_iters(iter + 1) + curr_gmres_iters;

    if iter + 1 == max_iter
        break;
    end


    % calculate BB1 step size
    sk = x1 - x0;
    yk = grad1 - grad0;
    sksk = sk' * sk;
    inner_prodct = sk' * yk;
    stepsize = sksk / inner_prodct;

    %calculate spectrum estimates
    mu1 = inner_prodct / sksk + 2*noise_tol/norm(sk);
    L1 = inner_prodct / sksk - 2*noise_tol/norm(sk);
    mu_estimates(iter) = min(mu0, mu1);
    L_estimates(iter) = max(L0, L1);
    
    %Now check if the noise control condition is violated and then update the noise tolerance
    if adaptive
        if -(grad1 - grad0)'*grad0 < 2*(1 + noise_control_param) * noise_tol*norm(grad0)
            noise_tol = noise_tol * noise_shrink_param;
            fprintf('Noise tolerance updated to: %e\n', noise_tol);
        end
    end


    %update variables
    x0 = x1;
    grad0 = grad1;
    mu0 = mu_estimates(iter);
    L0 = L_estimates(iter);
    iter = iter + 1;
end
iters = iter;
x = x1;

end

