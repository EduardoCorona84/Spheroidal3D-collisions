function [x, x_iters, gradient_evals, kkt_errors, mu_estimates, L_estimates, iters] = BBPGD_noise_tolerant(A, b, x0, max_iter, opt_params)
% May 2018 Wen Yan
%This does not resemble the original version really at all but it feels weird getting rid of attribution.

% Add the folder and all its subdirectories
%set parameters
opt_tol = opt_params.opt_tol;
noise = opt_params.noise;
method = opt_params.method;
noise_control_param = opt_params.noise_control_param;
line_search_param = opt_params.line_search_param;
projection = opt_params.projection;

% For 'affine true' method
if strcmp(method, 'affine true')
    L = opt_params.L;
    mu = opt_params.mu;
end

%initialize variables
x_iters = zeros(length(x0), max_iter);
gradient_evals = zeros(1, max_iter);
kkt_errors = zeros(1, max_iter);
mu_estimates = zeros(1, max_iter);
L_estimates = zeros(1, max_iter);
mu0 = inf; 
L0 = -inf;

%Compute first gradient
iter = 1;
grad0 = A(x0, noise) + b;
stepsize = 1;
gradient_evals(iter) = 1; % count the first gradient evaluation

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
    grad1 = A(x1, noise) + b;
    %Iterations are defined by the gradient evaluation, so store in the next entry of the array.
    gradient_evals(iter + 1) = gradient_evals(iter) + 1;

    if iter + 1 == max_iter
        break;
    end


    % calculate BB1 step size
    sk = x1 - x0;
    yk = grad1 - grad0;
    sksk = sk' * sk;
    inner_prodct = sk' * yk;

    %calculate spectrum estimates
    mu1 = inner_prodct / sksk + 2*noise/norm(sk);
    L1 = inner_prodct / sksk - 2*noise/norm(sk);
    mu_estimates(iter) = min(mu0, mu1);
    L_estimates(iter) = max(L0, L1);


    %check which noise tolerant method we are using
    switch method
        case 'none'
            % do nothing, use the standard BB step size
            stepsize = sksk / inner_prodct;

        case 'affine true'
            %use the actual spectrum to map the step size into the true spectrum
            a = (sksk * (L - mu)) / ((sksk * (L - mu)) + 4*norm(sk)*noise);
            b = L*sksk - a*(L*sksk + 2*norm(sk)*noise);
            stepsize = sksk / (a * inner_prodct + b);

        case 'affine estimate'
            %check the noise control parameter
            if -(grad1 - grad0)'*grad0 < 2*(1 + noise_control_param) * noise*norm(grad0)
                %if violated, use the spectrum estimates to map the step size into the true spectrum
                a = (sksk * (L_estimates(iter) - mu_estimates(iter))) / ((sksk * (L_estimates(iter) - mu_estimates(iter))) + 4*norm(sk)*noise);
                b = L_estimates(iter)*sksk - a*(L_estimates(iter)*sksk + 2*norm(sk)*noise);
                stepsize = sksk / (a * inner_prodct + b);
            
            else
                %do nothing, use the standard BB step size
                stepsize = sksk / inner_prodct;

            end

        case 'lengthening'
            if -(grad1 - grad0)'*grad0 < 2*(1 + noise_control_param) * noise*norm(grad0)
                %if violated, begin a line search
                new_yk = grad1 - grad0;
                new_line_search_param = line_search_param;
                while (-(new_yk)'*grad0 < 2*(1 + noise_control_param) * noise*norm(grad0))
                    new_grad = A(x0 - new_line_search_param * stepsize * grad0, noise) + b;
                    gradient_evals(iter + 1) = gradient_evals(iter + 1) + 1; % count the new gradient evaluation
                    new_yk = new_grad - grad0;
                    new_line_search_param = new_line_search_param * line_search_param;
                end
                sk = -new_line_search_param * stepsize * grad0;
                sksk = sk' * sk;
                inner_prodct = sk' * new_yk;    
                stepsize = sksk / inner_prodct;
                
            else
                %do nothing, use the standard BB step size
                stepsize = sksk / inner_prodct;

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

