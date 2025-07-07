function [x, iter, errStruct] = proxQuasiNewton(fcnGrad, x0, opts)
opts = defaultOpts(opts, length(x0));

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
    eta = min(-dot(q, grad_k)/ dot(q, opts.Q*q),1);
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
    [h0, H_k, U, V, opts] = updateHk(k, s_k, y_k, opts);

    % quasi-newton step
    p = -H_k * grad_k;
    % if isfield(opts, 'Q')
    %     % For QP, this is the optimal step length (see page 56 of n&W)
    % kappa = -dot(p, grad_k)/ dot(p, opts.Q*p);
        % c1 = 1e-4;
        % c2 = 0.9;
        % [f_test, grad_test] = fcnGrad(x_k + kappa * p);
        % assert( f_test <= f_k + c1*kappa*dot(grad_k, p), 'Sufficient decrease condition not satisfied');
        % assert( dot(grad_test, p) >= c2 *dot(grad_k, p), 'Curvature condition not satisfied');
    % end
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
    % if x_k ~= y
        % q = x_k - x_km1;
        % eta = min(-dot(q, grad_k)/ dot(q, opts.Q*q),1);
        % x_k = x_km1 + eta*q;
    % end
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
% U = zeros(n, r);
% V = zeros(n, r);
% H = diag(h0);
% B = diag(1./h0);
% for j = 1:r
%     s = S(:,j); 
%     y = Y(:,j); 
%     rho = 1 / dot(y,s);
%     if 1/rho < 1e-8
%         % TODO something better than skip if the curvature condition is bad
%         continue
%     end
%     UU = (eye(n) - rho*y*s');
%     H = UU'*H*UU + rho*(s*s');
% 
%     v = B*s / sqrt(dot(s, B*s));
%     u = y / sqrt(dot(y,s));
% 
%     B = B + u*u' - v*v' ;
%     U(:, j) = u;
%     V(:, j) = v;
% end
% % IDIOT CHECK 
% assert(norm(B*H - eye(n)) < 1e-8)
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