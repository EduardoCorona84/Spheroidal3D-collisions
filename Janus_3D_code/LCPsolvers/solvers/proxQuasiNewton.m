function [x, info] = proxQuasiNewton(fcnGrad, x0, opts)
opts = defaultOpts(opts, length(x0));
if ~isempty(opts.errFcn)
    if ishandle(opts.errFcn)
            info.errHist = zeros(opts.max_iter+1,1);
            info.errHist(1) = opts.errFcn(x0);
    elseif iscell(opts.errFcn)
        info.errHist = zeros(opts.max_iter+1,numel(opts.errFcn));
        for i = 1:numel(opts.errFcn)
            fcn = opts.errFcn{i};
            info.errHist(1,i) = fcn(x0);
        end
    end
end
x_k = x0;
[f_k, grad_k] = fcnGrad(x0);
[converged, info] = checkConvergence(0, f_k, x_k, ...
        grad_k, info, opts);
if converged
    x = x0;
    info.iter = 0;
    info.msg = "initial guess is sufficient";
    info.flag = -1;
    if ~isempty(opts.errFcn)
        info.errHist =  info.errHist(1,:);  
    end
    return 
end
x_km1 = x_k;
grad_km1 = grad_k;
kappa = opts.kappa;
% prox-grad descent for initial step
p = -grad_km1;
% For QP, this is the optimal step length (see page 56 of n&W)
% kappa = -dot(p, grad_k)/ dot(p, opts.Q*p);
y = x_km1 + kappa * p;
% if any(y < 0 )
%     x_k = max(0, y);
%     q = x_k - x0;
%     eta = min(-dot(q, grad_k)/ dot(q, opts.Q*q),1);
%     x_k = x0 + eta*q;
% else 
    x_k = y;
% end
[f_k, grad_k] = fcnGrad(x_k);
% quasi-newton iterations
kappa = opts.kappa;
for k = 1:opts.max_iter
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, info, opts);
    if converged
        x = x_k;
        break
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
    % f_km1 = f_k;
    [f_k, grad_k] = fcnGrad(x_k);
    % assert(f_k <= f_km1);   
end
info.iter = k;

if k == opts.max_iter
    warning('Maximum iterations reached')
    x = x_k;
    info.flag = 8; 
    info.msg = 'Maximum iterations';
end

if ~isempty(opts.errFcn)
   info.errHist =  info.errHist(1:k+1,:);  
end

end % sr1_custom

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
end
if ~isempty(U)  
    % TODO do this better with a shifted power method
    assert(min(eig(H)) > 0)
end

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

function [converged, info] = checkConvergence(k, f_k, x_k, ...
    grad_k, info, opts)
persistent old_kkt
if k >= opts.max_iter 
    converged = true; 
    return
end
% kkt conditions / LCP being satisified is equivalent to
% the grad_k(i) = 0 \perp x_k(i) = 0
phi = min(grad_k,x_k);
kkt = 0.5*dot(phi, phi);

converged = false;
if ~isempty(old_kkt) && (abs(kkt - old_kkt) / abs(kkt)) < opts.tol_rel  
    % Relative stopping criteria
    info.msg = "relative";
    info.flag = 3;
    converged = true;
elseif kkt < opts.tol_abs   
    % Absolute stopping criteria
    info.msg = "absolute";
    info.flag = 4;
    converged = true;
end
if ~isempty(opts.errFcn)
    if ishandle(opts.errFcn)
        info.errHist(k+1) = opts.errFcn(x_k);
    elseif iscell(opts.errFcn)
        for i = 1:numel(opts.errFcn)
            fcn = opts.errFcn{i};
            info.errHist(k+1,i) = fcn(x_k);
        end
    end
end
if converged 
    % include the stats for x0 and up to the kth iterate
    info.f = f_k;
    info.kkt = kkt;
end
old_kkt = kkt;

end % checkConvergence