function [s_k, y_k, opts] = noise_control(x_km1, x_k, s_k, grad_km1, grad_k, y_k, q, Aq, opts, residual_km1, residual_k, residual_step)
    %Check Noise Control Condition
    if y_k'*s_k >= (1 + opts.noise_control.parameter)*(residual_km1 + residual_k)*norm(s_k)
        %do nothing
    else
        switch opts.noise_control.type

        case 'simple_adaptive'
        %If the noise control condition is violated, gradients errors are reduced by a factor of opts.noise_control.reduction_parameter. Updates are not skipped. Assume no optimal step size
        opts.noise_control.noise = opts.noise_control.noise * opts.noise_control.reduction_parameter;

        case 'skip_adaptive'
        %If the noise control condition is violated, gradients errors are reduced by a factor of opts.noise_control.reduction_parameter. Updates are skipped. Assume no optimal step size
        opts.noise_control.noise = opts.noise_control.noise * opts.noise_control.reduction_parameter;
        opts.noise_control.skip = true; 

        case 'eig_adaptive'
        %gradient noise is chosen such that the lower bound on the BB step length is opts.noise_control.eig_adaptive_scaling * opts.noise_control.min_eig (min eigenvalue).
        opts.noise_control.skip = true; 
        
        case 'perturbed_adaptive'
        %gradient errors follow a forcing sequence like that of inexact newton, BB/uniform step sizes are used, curvature updates are skipped if the curvature condition is not satisfied. Assume no optimal step size
        opts.noise_control.skip = true;

        case 'independent_perturbed_adaptive'
        %gradient errors follow a forcing sequence like that of inexact newton, curvature information is updated independently of the gradient errors using the optimal step size, the noise level for this can be very large or small. 
        %To ensure a positive evaluation of pAp, we need to make sure the noise comparitively small compared to the p. Not sure how to actually choose this.
        beta_k = (2*(1 + opts.noise_control.parameter)*(residual_step)*norm(q))/(q'*Aq);
        if beta_k < 0
            opts.noise_control.skip = true;
        else
            corrected_x_k = x_km1 + beta_k * q;
            s_k = corrected_x_k - x_km1;
            y_k = beta_k * Aq;
        end

        end
    end

