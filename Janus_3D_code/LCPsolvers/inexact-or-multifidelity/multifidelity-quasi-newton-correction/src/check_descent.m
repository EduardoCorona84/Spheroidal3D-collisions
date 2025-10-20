function descent = check_descent(grad_0, p, Q, opts)
    %This function computes an upper bound on the error in the low fidelity evaluation and checks the descent condition
    
    %We want to compute a bound on the error of (A - A_hat)*p
    %Recall that we have that A_hat p is exact if p is in the span of the columns of S
    %Thus we we want to find the norm of the projection of p onto the orthogonal complement of the span of S
    %Suppose we have some estimate of the norm of A - A_hat, call abs_err
    %Thus the error would be bounded by abs_err * || P_{S_perp} p ||_2


    if opts.adaptive.enabled == false
        descent = true;
        return
    else
        %compute the projection of p onto the colummn space of Q
        p_Q = Q * (Q' * p);
        p_perp = p - p_Q;
        %compute the norm
        proj_norm = norm(p_perp);
        %compute the error bound
        error= opts.adaptive.abs_err_estimate * proj_norm;
        %compute upper bound for error
        error_bound = opts.adaptive.descent_parameter*(grad_0'*p)/norm(p);

        if error < error_bound
            descent = true;
        else
            descent = false;
        end

    end