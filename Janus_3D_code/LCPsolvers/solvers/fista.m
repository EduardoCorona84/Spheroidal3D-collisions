function [x, info] = fista(fg, x0, opts)

    % May 2018 Wen Yan
    % Edited 2025 Nic Rummel  
    % I used the projected gradient descent code as a template
    % But most of the actual implementation is from the FASTA paper
    %Largely based on the FASTA paper https://github.com/tomgoldstein/fasta-matlab
    %{
    Author = {Goldstein, Tom and Studer, Christoph and Baraniuk, Richard},
    Title = {A Field Guide to Forward-Backward Splitting with a {FASTA} Implementation},
    year = {2014},
    journal = {arXiv eprint},
    volume = {abs/1411.3406},
    url = {http://arxiv.org/abs/1411.3406},
    ee = {http://arxiv.org/abs/1411.3406}
    %}

    %The main difference now is that we have an acceleration step at the end
    % And we now have gradients corresponding to our feasible iterates (x_k) and the accelerated iterates (z_k).
    % We check convergence with x_k and grad_x_k, but compute the step with grad_z_k.
    [opts, info] = defaultLCPOpts(opts, x0);
    n = numel(x0); eta = 1; k = 0; tau = 1;
    x_k = x0; z_k = x_k; Ax_k = zeros(n,1); s = []; y = [];
    if any(x0 ~= 0) 
        Ax_k = opts.A(x_k); s = x_k; y = Ax_k;
    end
    [f_k, grad_x_k] = fg(x_k, Ax_k);
    % First step we have no acceleration
    Az_k = Ax_k;
    grad_z_k = grad_x_k;
    opts.acceleration.alpha_k = 1;
    while true
        %Check convergence (evaluate merit function) with grad_x_k
        [converged, info] = checkConvergence(k, f_k, x_k, ...
            grad_x_k, eta, info, opts);
        if converged
            x = x_k;
            return 
        end
        k = k + 1;
        x_km1 = x_k; Ax_km1 = Ax_k; f_km1 = f_k; grad_x_km1 = grad_x_k; z_km1 = z_k; Az_km1 = Az_k; grad_z_km1 = grad_z_k;
        % Gradient descent direction, we use grad_z_km1 here
        tau = stepSize(k,-grad_z_km1,x_km1,Ax_km1,opts,s,y);
        q = -tau * grad_z_km1;
        % Select step size
        prox_k = @(xtilde) max(xtilde,0);
        % forward backward takes a normal forward backward step and will return 
        % updated feasible iterates x_k and grad_x_k 
        step_k = @(t, opts) fistaFwdBwdstep(t, z_km1, Az_km1, x_km1, Ax_km1, q, prox_k, fg, opts);

        [x_k, Ax_k, f_k, grad_x_k] = linesearch(x_km1, Ax_km1, f_km1, grad_x_km1, ...
            step_k, opts);
        % These are in terms of the feasible iterates

        % Acceleration Step
        [z_k, Az_k, opts] = acceleration(x_k, x_km1, Ax_k, Ax_km1, opts);
        [~, grad_z_k] = fg(z_k, Az_k);

        switch lower(opts.acceleration.curvature)
            case 'accelerated'
                s = z_k - z_km1;
                y = grad_z_k - grad_z_km1;
            case 'feasible'
                s = x_k - x_km1;
                y = grad_x_k - grad_x_km1;
            case 'mixed'
                % z_k and x_k are the two most recent iterates, use one from each
                s = z_k - x_k;
                y = grad_z_k - grad_x_k;
        end
    end

end

function [x_k,Ax_k,f_k,grad_k, eta, p, Ap] = fistaFwdBwdstep(t, z_km1, Az_km1, x_km1, Ax_km1, q, prox, fg, opts)
    if t == 0
        x_k = x_km1;
        Ax_k = Ax_km1;
        [f_k, grad_k] = fg(x_k, Ax_k);
        eta = 0;
        p = q;
        Ap = zeros(numel(p),1);
        return 
    end
    ztilde = z_km1 + t * q;
    xhat = prox(ztilde);
    switch lower(opts.acceleration.direction)
        % Intuitively, feasible makes the most sense, I think. In practice, I'm not sure. 
        case 'feasible'
            p = xhat - x_km1;
            iterate_km1 = x_km1;
            evaluation_km1 = Ax_km1;
        case 'mixed'
            p = xhat - z_km1;
            iterate_km1 = z_km1;
            evaluation_km1 = Az_km1;
    end
    [eta, Ap] = stepSize(-1, p, iterate_km1, evaluation_km1, opts);
    x_k = iterate_km1 + eta*p;
    Ax_k = evaluation_km1 + eta*Ap;
    [f_k, grad_k] = fg(x_k, Ax_k);
end

function [z_k, Az_k, opts] = acceleration(x_k, x_km1, Ax_k, Ax_km1, opts)

    % In general, I don't think there is any guarantee that z_k is feasible.
    % It is not necessarily a convex combination
    switch opts.acceleration.method
        case 'known'
            beta_k = (sqrt(opts.acceleration.cond) - 1) / (sqrt(opts.acceleration.cond) + 1);
            z_k = x_k + beta_k*(x_k - x_km1);
            Az_k = Ax_k + beta_k*(Ax_k - Ax_km1);

            if opts.acceleration.restart
                % If we ascened, we restart the momentum and switch to adaptive.
                if dot(z_k - x_k, x_k - x_km1) > 0
                    opts.acceleration.alpha_k = 1;
                    opts.acceleration.method = 'adaptive';
                end
            end

        case 'adaptive'
            alpha_km1 = opts.acceleration.alpha_k;
            opts.acceleration.alpha_k = 1 + sqrt(1 + 4*alpha_km1^2)/2;
            beta_k = (alpha_km1 - 1)/opts.acceleration.alpha_k;
            s = x_k - x_km1;
            z_k = x_k + beta_k*s;
            Az_k = Ax_k + beta_k*(Ax_k - Ax_km1);
            if opts.acceleration.restart  && dot(z_k - x_k, s) > 0
                %If we ascened, we restart the momentum.
                opts.acceleration.alpha_k = 1;
            end

    end

end




