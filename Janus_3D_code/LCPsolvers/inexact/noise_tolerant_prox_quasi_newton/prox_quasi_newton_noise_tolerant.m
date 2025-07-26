function [x, iter, errStruct] = prox_quasi_noise_tolerant(fcnGrad, x0, opts)

    opts = defaultOpts(opts, length(x0));
    noise = opts.noise;
    x_k = x0;
    [f_k, grad_k] = fcnGrad(x0);
    x_km1 = x_k;
    grad_km1 = grad_k;
    errStruct = struct( ...
            'err',zeros(opts.max_iter+1,1), ...
            'f',zeros(opts.max_iter+1,1), ...
            'xk',zeros(length(x_k), opts.max_iter+1), ...
            'reason','');
    errStruct.err(1) = NaN;
    errStruct.f(1) = f_k;
    errStruct.xk(:,1) = x_k;
    kappa = opts.kappa;
    % prox-grad descent for initial step
    p = - grad_km1;
    % For QP, this is the optimal step length (see page 56 of n&W)
    % kappa = -dot(p, grad_k)/ dot(p, opts.Q*p);
    y = x_km1 + kappa * p;
    if any(y < 0 )
        x_k = max(0, y);
        q = x_k - x0;
        Qqvec = opts.Q(q);
        eta = min(-dot(q, grad_k)/ dot(q, Qqvec),1);
        x_k = x0 + eta*q;
    else 
        x_k = y;
    end
    [f_k, grad_k] = fcnGrad(x_k);
    % quasi-newton iterations
    kappa = opts.kappa;
    for k = 1:opts.max_iter

        [converged, errStruct] = checkConvergence(k, f_k, x_k, ...
            grad_k, errStruct, opts);
        if converged
            x = x_k;
            iter = k;
            return
        end

        s_k = x_k - x_km1;
        y_k = grad_k - grad_km1;
        
        %Check Noise Control Condition
        if opts.noise_control == true
            switch opts.noise_control_lengthening
                case 'simple'
                    switch opts.noise_control_line_search
                        case 'line search'
                            if (grad_k - grad_km1)'*q >= 2*(1 + opts.noise_control_parameter)*opts.noise*norm(q)
                            %do nothing
                            else
                                [s_k, y_k] = grad_search(x_km1, q, eta, fcnGrad, opts);
                            end
    
                        case 'quadratic'
                            denom = q'*Qqvec - noise*norm(q);
                                if denom >= 0
                                    min_val = (2*opts.noise_control_parameter*noise)/denom;
                                    if eta >= min_val
                                        % do nothing
                                    else
                                        %eta is the step length along the direction q, from x_km1
                                        eta = min_val;
                                        s_k = eta*q;
                                        [grad_k_new, ~] = fcnGrad(x_km1 + s_k);
                                        y_k = grad_k_new - grad_km1;
                                    end
                                else
                                    min_val = (2*opts.noise_control_parameter*noise)/denom;
                                    if eta <= min_val
                                        % do nothing
                                    else
                                        eta = min_val;
                                        s_k = eta*q;
                                        [grad_k_new, ~] = fcnGrad(x_km1 + s_k);
                                        y_k = grad_k_new - grad_km1;
                                    end
                                end
                        
                    end
                case 'proximal'
                    switch opts.noise_control_line_search
                        case 'line search'
                            % TODO implement line search
                            if (grad_k - grad_km1)'*q >= 2*(1 + opts.noise_control_parameter)*opts.noise*norm(q)
                                %do nothing
                            else
                                [s_k, y_k] = prox_search(x_km1, grad_km1, kappa, h0, U, V, fcnGrad, opts);
                            end
                        case 'quadratic'
                            denom = q'*Qqvec - noise*norm(q);
                            if denom >= 0
                                min_val = (2*opts.noise_control_parameter*noise)/denom;
                                if eta >= min_val
                                    % do nothing
                                else
                                    [s_k, y_k] = prox_search(x_km1, grad_km1, kappa, h0, U, V, fcnGrad, opts);
                                end
                            else
                                min_val = (2*opts.noise_control_parameter*noise)/denom;
                                if eta <= min_val
                                    % do nothing
                                else
                                    [s_k, y_k] = prox_search(x_km1, grad_km1, kappa, h0, U, V, fcnGrad, opts);
                                end
                            end  
                    end
            end
        end


        %
        [h0, H_k, U, V, opts] = updateHk(k, s_k, y_k, opts);

        % quasi-newton step
        p = -H_k * grad_k;
        %     % For QP, this is the optimal step length (see page 56 of n&W)
        kappa = -dot(p, grad_k)/ dot(p, opts.Q(p));
        x_km1 = x_k;
        y = x_k + kappa * p;
        % if kappa <= opts.kappa / 2^5
        %     kappa = opts.kappa;
        % end
        % % armijo linesearch TODO remove
        % for t = 1:10 
        %     xbar = x_k + kappa* p;
        %     [fnew, ~] = fcnGrad(xbar);
        %     % fpred = f_k + dot(p, grad_k + (diag(1 ./ h0) + V * diag(Lambda) * V')* p);
        %     if fnew > f_k %|| fnew / fpred > 2
        %         kappa = kappa / 2;
        %     else 
        %         break
        %     end
        % end
        % proximal step
        x_k = prox(y, h0, U, V, opts);
        q = x_k - x_km1;
        Qqvec = opts.Q(q);
        eta = min(-dot(q, grad_k)/ dot(q, Qqvec),1);
        x_k = x_km1 + eta*q;
        grad_km1 = grad_k;
        f_km1 = f_k;
        [f_k, grad_k] = fcnGrad(x_k);
        % assert(f_k <= f_km1);   
    end

    warning('Maximum iterations reached')
    x = x_k;
    iter = opts.max_iter;

    end % sr1_custom

    function opts = defaultOpts(opts, n)

        if ~isfield(opts, 'max_iter')
            opts.max_iter = 100;
        end

        if ~isfield(opts, 'tol_rel')
            opts.tol_rel = 1e-4;
        end

        if ~isfield(opts, 'tol_abs')
            opts.tol_abs = 1e-6;
        end

        if ~isfield(opts, 'gamma')
            opts.gamma = 0.8;
        end

        if ~isfield(opts, 'tau')
            % warning("Tau should be set as the 1/Lhat the the best estimate of the lipschitz constant of the the A mat");
            opts.tau = 1;
        end

        if ~isfield(opts, 'tau_min')
            opts.tau_min = 1e-14;
        end

        if ~isfield(opts, 'tau_max')
            opts.tau_max = Inf;
        end

        if ~isfield(opts, 'kappa')
            % warning("Tau should be set as the 1/Lhat the the best estimate of the lipschitz constant of the the A mat");
            opts.kappa = 1;
        end

        if ~isfield(opts, 'r')
            opts.r = 1;
        else 
            assert(opts.r > 0, 'Memory/effective-rank or hessian must by positive');
            opts.r = min(opts.r, n);
        end
        opts.S = zeros(n, opts.r);
        opts.Y = zeros(n, opts.r);
end % defaultOpts

function [h0, H, U, V, opts] = updateHk(k, s_k, y_k, opts)

assert(0 < opts.tau_min, "tau_min must be positive");
assert(opts.tau_min < opts.tau_max, "tau_max must be larger than tau_min");
assert(0 < opts.gamma && opts.gamma < 1, "gamma must be in (0,1)");
n= length(s_k);
tau_bb2 = dot(s_k,y_k) / norm(y_k,2)^2;
tau_bb2 = clip(tau_bb2, opts.tau_min, opts.tau_max);
if tau_bb2 == opts.tau_min
    warning('Convexity of cost function is stagnating'); 
end
h0 = opts.gamma * tau_bb2;

r = min(k, opts.r);
% Handeling the case where we have skipped previous hessian updates
skipped_r = find(arrayfun(@(i)all(opts.S(:,i) == 0),1:r),1,'first');
if ~isempty(skipped_r) && r < k
    r = max(skipped_r-1, 1);
end

% loose memory regardless of skipping nonsense
if k > opts.r
    opts.S(:,1:r-1) = opts.S(:,2:r);
    opts.Y(:,1:r-1) = opts.Y(:,2:r);
    opts.S(:,r) = 0;
    opts.Y(:,r) = 0;
end
rho = 1 / dot(y_k,s_k);
if 1/rho >= 1e-8
    opts.S(:,r) = s_k;
    opts.Y(:,r) = y_k;
else
    r = r -1;
end
S = opts.S(:,1:r);
Y = opts.Y(:,1:r);
% ApD = zeros(r,r);
% for i = 1:r
%     for j = 1:r
%         ApD(i,j) = dot(S(:,i), Y(:,j));  
%     end
% end
% 
% BB = chol(S'*diag(1 ./ h0)*S - ApD); % for B
% HH = chol(ApD' - Y'*diag(h0)*Y); % for H
% V = (Y - diag(1 ./ h0)*S) / BB;
% U = (S - diag(h0)*Y) / HH;

%% BFGS 
U = zeros(n, r);
V = zeros(n, r);
H = diag(h0);
B = diag(1./h0);
for j = 1:r
    s = S(:,j); 
    y = Y(:,j); 
    rho = 1 / dot(y,s);
    if 1/rho < 1e-8
    % TODO something better than skip if the curvature condition is bad
        continue
    end
    UU = (eye(n) - rho*y*s');
    H = UU'*H*UU + rho*(s*s');

    v = B*s / sqrt(dot(s, B*s));
    u = y / sqrt(dot(y,s));
 
    B = B + u*u' - v*v' ;
    U(:, j) = u;
    V(:, j) = v;
end
% % IDIOT CHECK 
% assert(norm(B*H - eye(n)) < 1e-8)
%{
%% SR1
U = zeros(n, r);
V = zeros(n, r);
H = eye(n) * h0;
B = eye(n) / h0;
for j = 1:r 
    s = S(:,j); 
    y = Y(:,j);
    % update H (and mem)
    denomH = dot(s - H*y,y);
    sigma = sign(denomH);
    denomH = sqrt(sigma*denomH);
    assert(isreal(denomH));
    u = (s - H*y) / denomH;
    H = H + sigma*(u*u');
    % update H (and mem)
    denomB = dot(y - B*s,s);
    lambda = sign(denomB);
    denomB = sqrt(lambda*denomB);
    assert(isreal(denomB));
    v = (y - B*s) / denomB;
    B = B + lambda*(v*v');
    % Store for proof of concept
    if lambda > 0 
        U(:,j) = v;
    else
        V(:,j) = v;
    end
end
mask = arrayfun(@(i) all(U(:,i) == 0),1:r);
U = U(:,~mask);
V = V(:,mask);
assert(norm(B*H - eye(n))/norm(H)/norm(B) < 1e-8)
%}
end % updateHk

function xstar = prox(y, h0, U, V, opts)
if size(U,2) + size(V,2) == 1
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
end 
xstar = prox_rankr(y, h0, U, V, opts);
end
function [converged, errStruct] = checkConvergence(k, f_k, x_k, ...
    grad_k, errStruct, opts)

if k >= opts.max_iter 
    converged = true; 
    return
end
% kkt conditions / LCP being satisified is equivalen to
% the grad_k(i) = 0 perp  x_k(i) = 0
phi = min(grad_k,x_k);
old_err = errStruct.err(k);
err = 0.5*dot(phi, phi);
errStruct.err(k+1) = err;
errStruct.f(k+1) = f_k;
errStruct.xk(:,k+1) = x_k;
converged = false;

if (abs(err - old_err) / abs(old_err)) < opts.tol_rel  % Relative stopping criteria
    errStruct.reason = "relative";
     converged = true;
elseif err < opts.tol_abs   % Absolute stopping criteria
    errStruct.reason = "absolute";
    converged = true;
end
if converged 
    % include the stats for x0 and up to the kth iterate
    errStruct.f = errStruct.f(1:k+1);
    errStruct.err = errStruct.err(1:k+1);
    errStruct.xk = errStruct.xk(:,1:k+1);
end

end % checkConvergence

function [sk, yk] = prox_search(x0, grad0, kappa0, h0, U, V, fcnGrad, opts)

    noisey = true;
    iters = 0;
    while noisey && iters < 10000
        %increase step size inside the prox
        kappa0 = kappa0 * opts.noise_control_line_search_parameter;
        y = x0 + kappa0 * grad0;

        %compute the prox
        x1 = prox(y, h0, U, V, opts);
        [~, grad1] = fcnGrad(x1);

        %Check noise control with the new prox gradient and direction
        if ((grad1 - grad0)'*(x1 - x0)) > 2*(1 + opts.noise_control_parameter)*opts.noise*norm(x1 - x0) 
            noisey = false;
        end
        iters = iters + 1;
    end

    %update the curvature pair
    sk = x1 - x0;
    yk = grad1 - grad0;

end

function [sk, yk] = grad_search(x0, q, eta0, fcnGrad, opts)

    noisey = true;
    iters = 0;
    while noisey && iters < 10000
        %increase step size
        eta0 = eta0 * opts.noise_control_line_search_parameter;
        x1 = x0 + eta0 *q;

        %compute the gradient at the new point
        [~, grad1] = fcnGrad(y);

        %Check noise control with the new gradient and direction
        if ((grad1 - grad0)'*(q)) > 2*(1 + opts.noise_control_parameter)*opts.noise*norm(q) 
            noisey = false;
        end
        iters = iters + 1;
    end

    %update the curvature pair
    sk = x1 - x0;
    yk = grad1 - grad0;

end