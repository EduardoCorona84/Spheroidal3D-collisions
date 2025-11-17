function [x, info] = alternating_high_and_low(fg, x0, opts)
    %high fidelity initialization
    [opts.outer, info] = defaultLCPOpts(opts.outer, x0);
    checkOpts(opts.outer)
    n = numel(x0);
    x_k = x0;
    if all(x0 == 0)
        Ax_k = zeros(n,1);
    else 
        Ax_k = opts.outer.A(x_k);
    end
    [f_k, grad_k] = fg(x_k, Ax_k);
    eta = 1;
    s = [];
    y = [];
    k = 0;
    while true

        [converged, info] = checkConvergence(k, f_k, x_k, ...
            grad_k, eta, info, opts.outer);
        if converged
            x = x_k;
            break
        end
        k = k + 1;
        x_km1 = x_k;
        grad_km1 = grad_k;
        Ax_km1 = Ax_k;
        [H, h0, U, V] = updateHk(k, s, y, opts.outer);
        % quasi-newton step direction
        p = -H(grad_k);
        % step size direction
        kappa = stepSize(k, p, x_km1, Ax_km1, opts.outer); 
        x_k_half = x_km1 + kappa * p;
        x_k_half = prox(x_k_half, h0, U, V, opts.outer);

        % Now we will take the inner low fidelity step
        %But only if we have two high fidelity gradient evaluations
        if opts.inner.enabled == true  
            [x_k, ~] = inner_step(k, x_k_half, x_km1, Ax_km1, kappa, H, h0, U, V, s, y, opts);
        else
            x_k = x_k_half;
        end
        % Possibly a step length update after the projection
        q = x_k - x_km1;
        [eta, Aq] = stepSize(-1, q, x_km1, Ax_km1, opts.outer);
        x_k = x_km1 + eta*q;
        Ax_k = Ax_km1 + eta*Aq;
        % Increment the number of iterations
        
        [f_k, grad_k] = fg(x_k, Ax_k);
        s = x_k - x_km1;
        y = grad_k - grad_km1;
    end

end % proxQuasiNewton

function [H, h0, U, V] = updateHk(k, s, y_k, opts)

    switch lower(opts.qnUpdate)
        case 'bfgs'
            [H, h0, U, V] = get_H_BFGS(k, s, y_k, opts);
        case 'sr1'
            
        otherwise
            error([opts.qnUpdate ' update not implement'])
    end
end % updateHk

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

function [x_k, info_low] = inner_step(k, x_k_half, x_km1, Ax_km1, kappa, H, h0, U, V, s, y, opts)

    grad_low = opts.inner.A(x_k_half - x_km1) + Ax_km1 + opts.outer.b;
    if k == 1 
        %On the first iteration we just take a constant step size gradient descent step.
        p = -grad_low;
        kappa = opts.inner.stepSize.init;
        x_k = x_k_half + opts.inner.damping*kappa * p;
        x_k = max(0, x_k);

    else
        switch lower(opts.inner.solver)
            case 'proxquasinewton'
                %reuse the hessian approximation (and step size) from the high fidelity
                p = -H(grad_low);
                %The opts.inner.damping parameter can be used to scale the size of the step size (as its low fidelity we may not trust it as much)
                x_k = x_k_half + opts.inner.damping*kappa * p;
                x_k = prox(x_k, h0, U, V, opts.outer);
                
            case 'bbpgd'
                p = -grad_low;
                %The opts.inner.damping parameter can be used to scale the size of the step size (as its low fidelity we may not trust it as much)
                kappa = (s'*s)/(s'*y);
                x_k = x_k_half + opts.inner.damping*kappa * p;
                x_k = max(0, x_k);

            case 'apgd'

            otherwise
                error(['Inner solver ' opts.inner.solver ' not implemented']);
        end
    end
    %First we need to compute the affine corrected low fidelity gradient
    %Recall Ax_k_half = A(x_km1) + A(x_k_half - x_km1)
    %We will compute A(x_k_half - x_km1) via the low fidelity
    grad_low = opts.inner.A(x_k_half - x_km1) + Ax_km1 + opts.outer.b;
    
    info_low = []; %placeholder for now
end

function checkOpts(opts)
    assert(strcmpi(opts.stepSize.kappa, 'opt') ...
        || strcmpi(opts.stepSize.kappa, 'uniform'),...
        ['Quasi Newtwon method should use a uniform step size of 1 because the '...
        'BB step size is baked into the hessian approximation (unconstrained) optimal size']);
    assert(0 < opts.tau_min, "tau_min must be positive");
    assert(opts.tau_min < opts.tau_max, "tau_max must be larger than tau_min");
    assert(0 < opts.gamma && opts.gamma < 1, "gamma must be in (0,1)");

end