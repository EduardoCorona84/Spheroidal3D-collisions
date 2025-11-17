function [x, info] = multifidelity_quasi_newton(fg, x0, opts)

    %we don't really need this, but for now we will just use the same opts (really only care about tolerance, max iters, storage, and optimal step size)
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
        %Quasi-Newton update
        update_matrix = update_Bk(k, s, y, opts);
        %In this function, we want to solve the subproblem using the quasi-Newton updated low fidelity model.
        %During the solve, we will store evaluations of the low fidelity model
        %and construct a quasi-Newton approximation to it as well.
        %We will pass these evaluations between iterations as essentially enhanced memory of the low fidelity.
        x_k = solve_subproblem(x_km1, grad_km1, update_matrix, opts);

        % This can all stay the same
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

function update_matrix = update_Bk(k, s, y, opts)
    %For now, we will just use a dense quasi-Newton update (BFGS)
    persistent high_correction;
    if k == 1
        n = length(opts.outer.b);
        high_correction = zeros(n,n);
        %do nothing on first iteration besides initialization
    else
        quantity = opts.inner.A(s) + high_correction*s;
        high_correction = high_correction + (y*y')/(y'*s) - (quantity*quantity')/(s'*quantity);
    end

    update_matrix = high_correction;

end

function x_k = solve_subproblem(x_km1, grad_km1, high_correction, opts)
    %Solve the subproblem using a quasi newton updated low fidelity model
    %Specifically, solve min 0.5 z' B_k z + grad_km1'z s.t. z >= -x_km1

    %Create a function handle for the low fidelity model with quasi-Newton update

    % We are going to use a low fidelity matvec here, that we normally would not need to use, just for compatibility with proxQuasiNewton
    clear sub_opts;
    sub_opts.A = @(x) opts.inner.A(x) + high_correction*x;
    sub_opts.b = grad_km1 - sub_opts.A(x_km1);


    fg_sub = @(x, Ax) deal(0.5 * x' * Ax + sub_opts.b' * x, Ax + sub_opts.b);

    %set sub problem opts struct
    sub_opts.solver = 'proxquasinewton';
    sub_opts.stepSize.init = 'uniform';
    sub_opts.stepSize.kappa = 'uniform';
    sub_opts.stepSize.eta = 'opt';
    sub_opts.storeIts = false;
    sub_opts.max_iter = opts.inner.max_iter;
    sub_opts.qn_memory.warm = false;
    sub_opts.qn_memory.store = false;
    %Now we are going to solve the subproblem using proxQuasiNewton
    %TO DO: Add warm starting logic for proxQuasiNewton quasi-Newton approximation

    [x_k, ~] = proxQuasiNewton(fg_sub, x_km1, sub_opts);

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