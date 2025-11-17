function [x, info, opts] = multifidelityProxQuasiNewton(fg, x0, opts)

    %Warm Start with Low Fidelity/Need to look into this more.
    opts.sub.max_iter = 4;
    AHatCnt = @(x) Acounter(x, AHat, false);
    ACnt = opts.A;
    opts.A = AHatCnt;
    fgMid =  @(x, Ax) quadraticLoss(x, AHatCnt, b, Ax);
    [x0, ~, opts.sub] = proxQuasiNewton(fgMid, x0, opts.sub);

    % For the outer part of the struct, we will use default opts
    % The memory will be the same size as max_iter (we will set this to be small).
    opts.high.m = opts.high.max_iter;
    [opts.high, info] = defaultLCPOpts(opts.high, x0);
    %checkOpts(opts)
    n = numel(x0); k = 0;
    x_k = x0; Ax_k = zeros(n,1); s = []; y = [];
    if ~all(x0 == 0) 
        % We will treat the warm start as a big step
        Ax_k = opts.high.A(x_k);
        s = x_k;
        y = Ax_k;
        % For saving computation, we will store the low fidelity evaluation.
        Ahats = opts.sub.Ax_k;
    end
    [f_k, grad_k] = fg(x_k, Ax_k);
    while true
        [converged, info] = checkConvergence(k, f_k, x_k, ...
            grad_k, [], info, opts.high);
        if converged
            x = x_k; 
            break
        end
        % Increase k 
        k = k + 1;
        x_km1 = x_k; Ax_km1 = Ax_k; f_km1 = f_k; grad_km1 = grad_k; Ahatx_km1 = Ahatx_k;
        % Update the low fidelity Hessian with high fidelity data.
        % The implementation for prox quasi-Newton will not work. Need something here.
        [~, opts] = updateBk(s, y, Ahats, opts);
        % This solves the subproblem and saves the new low fidelity memory. 
        [x_k, Ahatx_k, opts, info] = solveSubProblem(grad_km1, opts, info);
        
        % Use the optimal step size 
        p = x_k - x_km1;
        [eta, Ap] = stepSize(-1, p, x_km1, Ax_km1, opts);
        x_k = x_km1 + eta*p;
        Ax_k = Ax_km1 + eta*Ap;
        [f_k, grad_k] = fg(x_k, Ax_k);
        % Update the low fidelity Ax_k with the optimal step size
        Ahatx_k = Ahatx_km1 + eta*(Ahatx_k - Ahatx_km1);
        opts.sub.Ax_k = Ahatx_k;
        % Save secant conditions
        % High
        s = x_k - x_km1;
        y = grad_k - grad_km1;
        % Low
        Ahats = Ahatx_k - Ahatx_km1;
    end

end 


function opts = updateBk(s, y, Ahats, opts)
    
    switch(opts.high.qn.update)
        case 'sr1'
            
        case 'bfgs'
            % Because there are so few iterations, we can store BFGS with the unrolled update and not use the compact representation. 
            % (Memory does not fall out of the window).
            % Will use the S and Y arrays for U and V instead
            U = opts.high.qn.S;
            V = opts.high.qn.Y;
            % Find the first empty (zeros) column.
            r = find(all(U == 0, 1), 1, 'first');
            % Compute the Ahat_k(s) matvec
            quantity = Ahats + U(1:r - 1)*(U(1:r - 1)'*s) - V(1:r - 1)*(V(1:r - 1)'*s);
            % Form and store the BFGS update (probably want to add a curvature check here)
            U(:, r) = y/sqrt(y'*s);
            V(:, r) = quantity/sqrt(s'*quantity);
            opts.high.qn.S = U;
            opts.high.qn.Y = V; 
            opts.high.qn.r = r;

        otherwise
            error(['Not implemented update ' opts.high.qn.update])
    end

end

function [x_k, Ahatx_k, opts, info] = solveSubProblem(x_km1, grad_km1, opts)

    %{
    We are storing the low fidelity secant conditions/memory in opts.sub. We are assuming these are of the form:
    Ahat_0S = Y
    That is, we are storing just the action of the low fidelity with no rank updates.
    %}

    % The high fidelity updates of the low are stored in S and Y (which we will denote U and V)
    U = opts.high.qn.S;
    V = opts.high.qn.Y;
    r = find(all(U == 0, 1), 1, 'first');
    U = U(:, 1:r - 1);
    V = V(:, 1:r - 1);
    
    % We are going to use the proxQuasiNewton function to solve the subproblem
    % Specifically with the high fidelity updated Ahat_0 (Ahat_k).
    % This means we need to update opts.sub.qn.Y such that it has the action of the rank updates
    opts.sub.qn.Y = opts.sub.qn.Y + U * (U' * opts.sub.qn.S) - V * (V' * opts.sub.qn.S);
    opts.sub.qn.rho = reshape(sum(opts.sub.qn.S.*opts.sub.qn.Y, 1).^-1, [], 1);

    % We also will update the last iterate, as we use it in constructing opts.sub.b
    opts.sub.Ax_k = opts.sub.Ax_k + U * (U' * opts.sub.Ax_k) - V * (V' * opts.sub.Ax_k);

    % Sub Problem Opts
    opts.sub.A(x) = opts.low.A(x) + opts.outer.qn.S(:, opts.outer.qn.r)*(opts.outer.qn.S(:, opts.outer.qn.r)'*x) - opts.outer.qn.Y(:, opts.outer.qn.r)*(opts.outer.qn.Y(:, opts.outer.qn.r)'*x);
    opts.sub.b = grad_km1 - opts.sub.Ax_k;

    % TODO: Set subproblem tolerance and max_iter
    opts.sub.max_iter = opts.low.max_iter;

    % Maybe use quadratic loss here?
    fg_sub = @(x, Ax) deal(1/2 * x' * Ax, Ax + opts.sub.b); 

    % Solve the subproblem with the sub.opts
    [x_k, info, opts.sub] = proxQuasiNewton(fg_sub, x_km1, opts.sub);

    % After solving the subproblem, we need to update the secant conditions again to remove the low rank updates (and the last iterate).
    opts.sub.qn.Y = opts.sub.qn.Y - U * (U' * opts.sub.qn.S) + V * (V' * opts.sub.qn.S);

    opts.sub.Ax_k = opts.sub.Ax_k - U * (U' * opts.sub.qn.S) + V * (V' * opts.sub.qn.S);
    Ahatx_k = opts.sub.Ax_k; %


end