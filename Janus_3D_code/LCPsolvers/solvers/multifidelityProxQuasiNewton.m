function [x, info, opts] = multifidelityProxQuasiNewton(fg, x0, opts)
    %Warm Start with Low Fidelity/Need to look into this more.
    fgMid =  @(x, Ax) quadraticLoss(x, opts.low.A, opts.low.b, Ax);
    opts.sub.b = opts.low.b;
    opts.sub.A = opts.low.A;
    opts.sub = defaultLCPOpts(opts.sub, x0);
    if opts.low.initWithLofi
        [x0, ~, opts.sub] = proxQuasiNewton(fgMid, x0, opts.sub);
        Ahatx_k = opts.sub.Ax_k;
    else
        Ahatx_k = opts.sub.A(x0);
        opts.sub.Ax_k = Ahatx_k;
    end
    % For the outer part of the struct, we will use default opts
    % The memory will be the same size as max_iter (we will set this to be small).
    opts.m = opts.max_iter;
    [opts, info] = defaultLCPOpts(opts, x0);
    %checkOpts(opts)
    n = numel(x0); k = 0;
    x_k = x0; Ax_k = opts.A(x_k);
    % the first secant condition can is free
    % s = (x0 - 0), y = A[x0] - A[0] = A[x0]
    s = x_k; y = Ax_k;
    % For efficiency also cache low fidelity evaluation of \hat{A}[s].
    Ahats = Ahatx_k;
    x_km1=[]; grad_km1=[];
    [f_k, grad_k] = fg(x_k, Ax_k);
    while true
        [converged, info] = checkConvergence(k, f_k, x_k, ...
            grad_k, [], info, opts, x_km1, grad_km1);
        if converged
            k
            x = x_k; 
            break
        end
        % Increase k 
        k = k + 1;
        x_km1 = x_k; Ax_km1 = Ax_k; f_km1 = f_k; grad_km1 = grad_k; Ahatx_km1 = Ahatx_k;
        % Update the low fidelity Hessian with high fidelity data.
        % The implementation for prox quasi-Newton will not work. Need something here.
        opts = updateBk(s, y, Ahats, opts);
        % This solves the subproblem and saves the new low fidelity memory. 
        [xhat_k, Ahatx_k, opts, ~] = solveSubProblem(x_km1, grad_km1, opts);
        
        % Use the optimal step size 
        p = xhat_k - x_km1;
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
    if isempty(s) || isempty(y) || isempty(Ahats) || all(s == 0)
        opts.qn.r = [];
        opts.sub.qn.r2 = [];
        return 
    end
    switch(opts.qn.update)
        case 'sr1'
            
        case 'bfgs'
            %curvature check 
            if y'*s < 1e-8 * norm(s)*norm(y)
                return;
            end
            % Because there are so few iterations, we can store BFGS with the unrolled update and not use the compact representation. 
            % (Memory does not fall out of the window).
            U = opts.qn.U;
            V = opts.qn.V;
            % Find the first empty (zeros) column.
            r = find(all(U == 0, 1), 1, 'first');
            % Compute the B^{(k)}(s) = (\hat{A} + UU^\top - VV^\top)s 
            if r == 1
                Bs = Ahats;
            else
                Bs = Ahats + U(:, 1:r - 1)*(U(:, 1:r - 1)'*s) - V(:, 1:r - 1)*(V(:, 1:r - 1)'*s);
            end

            % Form and store the BFGS update 
            rho = 1/dot(y,s);
            if rho <= 1e8
                U(:, r) = y*rho;
                V(:, r) = Bs/sqrt(s'*Bs);
            elseif r == 1
                r = [];
            else
                r = r-1;
            end
            opts.qn.U = U;
            opts.qn.V = V; 
            opts.qn.r = r;
        otherwise
            error(['Not implemented update ' opts.qn.update])
    end

end

function [x_k, Ahatx_k, opts, info] = solveSubProblem(x_km1, grad_km1, opts)
    %{
    We are storing the low fidelity secant conditions/memory in opts.sub. We are assuming these are of the form:
    Ahat_0S = Y
    That is, we are storing just the action of the low fidelity with no rank updates.
    %}

    % The low fidelity approximation of A can be updated to be \hat{A}^\hat{A} + UU^\top - VV^\top 
    U = opts.qn.U;
    V = opts.qn.V;
    r = opts.qn.r;
    U = U(:, 1:r);
    V = V(:, 1:r);
    
    % Using proxQuasiNewton to solve the subproblem benefits from passing 
    % secant conditions from one problem to the next.
    % \argmin_{x>0} 1/2 x^\top B^{(k)}x + x^\top c
    % The memory stored is secant conditions for \hat{A} : Y = \hat{A} S 
    % We want secant conditions for B : Y = B S
    S = opts.sub.qn.S;
    AhatS = opts.sub.qn.Y;
    r2 = find(~all(S == 0, 1), 1, 'last');     
    S = S(:,1:r2);
    AhatS = AhatS(:,1:r2);
    BS = AhatS + U * (U' * S) - V * (V' * S);
    rho = reshape(sum(S.*BS, 1).^-1, [], 1);
    
    % Fill the opts struct for the subproblem
    B = @(x) opts.low.A(x) + U*(U'*x) - V*(V'*x);
    Bx_km1 = B(x_km1);
    c = grad_km1 - Bx_km1;
    opts.sub.A = B;
    opts.sub.b = c;
    opts.sub.qn.Y(:,1:r2) = BS;
    opts.sub.qn.rho(1:r2) = rho;
    opts.sub.Ax_k = Bx_km1;
    fg_sub = @(x, Bx) quadraticLoss(x, B, c, Bx);


    % n = numel(c);
    % BB = B(eye(n));
    % BB = (BB + BB')/2;
    % Bx_km1 = BB*x_km1;
    % c = grad_km1 - Bx_km1;
    % cvx_begin quiet
    %     variable z(n)
    %     minimize 1/2*dot(z, BB*z) + dot(z,c)
    %     subject to
    %     0 <= z %#ok<NODEF,NOPRT>
    % cvx_end
    % cvx_begin quiet
    %     variable y(n)
    %     minimize 1/2*dot(y - x_km1, BB*(y-x_km1)) + dot(y-x_km1, grad_km1)
    %     subject to
    %     0 <= y %#ok<NODEF,NOPRT>
    % cvx_end
    % norm(y - z) / norm(z)

    % fprintf('relErr BS = %.6g\n', norm(BS - B(S)) / norm(B(S)))
    % Solve the subproblem 
    % n = numel(c);
    % opts.sub.qn.S = zeros(n,n);
    % opts.sub.qn.Y = zeros(n,n);
    [x_k, info, opts.sub] = proxQuasiNewton(fg_sub, x_km1, opts.sub);
    % norm(x_k - z) / norm(z)
    Bx_k = opts.sub.Ax_k;
    % After solving the subproblem, in order for bookkeeping to be simplified,
    % we need to store the secant conditions for \hat{A} (not B).
    % \hat{A}S = BS - U*U^\topS + V*V^\top*S
    S = opts.sub.qn.S;
    BS = opts.sub.qn.Y;
    r2 = find(~all(S == 0, 1), 1, 'last');
    S = S(:,1:r2);
    BS = BS(:,1:r2);
    AhatS = BS - U * (U' * S) + V * (V' * S);
    opts.sub.qn.Y(:, 1:r2) = AhatS;
    Ahatx_k = Bx_k - U * (U' * x_k) + V * (V' * x_k);
    %% idiot checks
    % fprintf('relErr BS = %.6g\n', norm(BS - B(S)) / norm(B(S)))
    % fprintf('relErr AhatS = %.6g\n', norm(AhatS - opts.low.A(S)) / norm(opts.low.A(S)))
    % fprintf('relErr Ahatx_k = %.6g\n', norm(Ahatx_k - opts.low.A(x_k)) / norm(opts.low.A(x_k)))
end