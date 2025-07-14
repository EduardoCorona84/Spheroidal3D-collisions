function [x, info] = proxQuasiNewton(fcnGrad, x0, opts)
[opts, info] = defaultOpts(opts, x0);
checkOpts(opts)
n = numel(x0);
kappa = 1;
x_k = x0;
x_km1 = NaN*ones(n,1);
grad_km1 = NaN*ones(n,1);
k = 0;
while true
    [f_k, grad_k] = fcnGrad(x_k);
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, kappa, info, opts);
    if converged
        x = x_k;
        break
    end
    s_k = x_k - x_km1;
    y_k = grad_k - grad_km1;
    x_km1 = x_k;
    grad_km1 = grad_k;
    [h0, H, U, V, opts] = updateHk(k, s_k, y_k, opts);
    % quasi-newton step direction
    p = -H(grad_k);
    % step size direction
    kappa = stepSize(k, p, grad_k, opts); 
    x_k = prox(x_km1 + kappa * p, h0, U, V, opts);    
    % Possibly a step length update after the projection
    q = x_k - x_km1;
    eta = min(1, stepSize(-1, q, grad_k, opts));
    x_k = x_km1 + eta*q;
    % Increment the number of iterations
    k = k + 1;
end

end % proxQuasiNewton

function [h0, H, U, V, opts] = updateHk(k, s_k, y_k, opts)

if k == 0 
    % Initial step we do not compute any curvature information because only
    % information about x0 is known, thus no secant equaion could be
    % satisfied
    h0 = 1; 
    H = @(x) x; 
    U = []; 
    V = [];
    return 
end

assert(0 < opts.tau_min, "tau_min must be positive");
assert(opts.tau_min < opts.tau_max, "tau_max must be larger than tau_min");
assert(0 < opts.gamma && opts.gamma < 1, "gamma must be in (0,1)");
n = length(s_k);
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

% Let memory fall out of context window regardless of wether we skip or not
if k > opts.r
    opts.S(:,1:r-1) = opts.S(:,2:r);
    opts.Y(:,1:r-1) = opts.Y(:,2:r);
    opts.S(:,r) = 0;
    opts.Y(:,r) = 0;
end
% Curvature check
rho = 1 / dot(y_k,s_k);
if 1/rho >= 1e-8
    opts.S(:,r) = s_k;
    opts.Y(:,r) = y_k;
else
    r = r - 1;
end
S = opts.S(:,1:r);
Y = opts.Y(:,1:r);

switch lower(opts.qnUpdate)
    case 'bfgs'
        U = zeros(n, r);
        V = zeros(n, r);
        H = eye(n) * h0;
        B = eye(n) / h0;
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
        % IDIOT CHECK 
        assert(norm(B*H - eye(n)) < 1e-8)
        % TODO: Make this matrix free...
        H = @(x) H*x;
    case 'sr1'
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
        % TODO: Make this matrix free...
        H = @(x) H*x;
        % We may have an indefinite matrix from the SR1 update
        if ~isempty(U)  
            % TODO do this better with a shifted power method
            assert(min(eig(H)) > 0)
        end
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

function checkOpts(opts)
assert(strcmpi(opts.kappa.fwd, 'opt') ...
    || strcmpi(opts.kappa.fwd, 'uniform'),...
    ['Quasi Newtwon method should use a uniform step size of 1 because the '...
    'BB step size is baked into the hessian approximation (unconstrained) optimal size']);
end