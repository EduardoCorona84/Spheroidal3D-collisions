function [x, info] = accelerated_prox_grad(fg, x0, opts)
    %This is an accelerated proximal gradient method. 
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

    [opts, info] = defaultLCPOpts(opts, x0);
    checkOpts(opts)
    n = numel(x0);
    eta = 1;
    x_k = x0;
    if all(x0 == 0)
        Ax_k = zeros(n,1);
    else 
        Ax_k = opts.A(x_k);
    end
    s = [];
    y = [];
    [f_k, grad_k] = fg(x_k, Ax_k);
    k = 0;
    alpha_k = 1;
    while true
        [converged, info] = checkConvergence(k, f_k, x_k, ...
            grad_k, eta, info, opts);
        if converged
            x = x_k;
            return 
        end
        k = k + 1;
        x_km1 = x_k; Ax_km1 = Ax_k; grad_km1 = grad_k;
        % Gradient descent direction
        p = -grad_k;
        % Select step size
        kappa = stepSize(k, p, x_km1, Ax_km1,  opts, s, y); 
        % Gradient Descent Step
        x_k = x_km1 + kappa * p;
        % Projection to positive orthant
        x_k = max(0, x_k);
        %add backing tracking logic here
        x_k = back_track(fg, x_k, x_km1, kappa, p, opts);
        % Possibly a step length update after the projection
        q = x_k - x_km1;
        [eta, Aq] = stepSize(-1, q, x_km1, Ax_km1, opts);
        Ax_k = Ax_km1 + eta*Aq;
        x_k = x_km1 + eta*q;
        %compute f_k
        f_k = (1/2)*x_k'*Ax_k + opts.b'*x_k;
        %acceleration step
        alpha_km1 = alpha_k;
        alpha_k = 1 + sqrt(1 + 4*alpha_km1^2)/2;
        %intermediate acceleration step
        z_k = x_k + (alpha_km1 - 1)/alpha_k * (x_k - x_km1);
        Az_k = Ax_k + (alpha_km1 - 1)/alpha_k * (Ax_k - Ax_km1);
        if opts.acceleration.restart && dot(z_k - x_k, x_k - x_km1) > 0
            alpha_k = 1;
        end
        %compute grad_k
        [~, grad_k] = fg(z_k, Az_k);
        s = x_k - x_km1;
        y = grad_k - grad_km1;
    end

end

function checkOpts(opts)
assert(strcmpi(opts.stepSize.kappa, 'opt') ...
    || contains(lower(opts.stepSize.kappa), 'bb'),...
    'Projected Gradient Descent should use the the BB step size or the (unconstrained) optimal size');
end
