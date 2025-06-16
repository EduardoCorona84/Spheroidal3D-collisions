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
% prox-grad descent for initial step
p0 = - opts.kappa * opts.tau * grad_km1;
x_k = max(0, x_km1 + p0);
[f_k, grad_k] = fcnGrad(x_k);

% quasi-newton iterations
kappa = opts.kappa;
for k = 1:opts.max_iter
    [converged, errStruct] = checkConvergence(k, f_k, x_k, ...
        grad_k, errStruct, opts);
    if converged
        x = x_k;
        iter = k-1;
        return
    end

    s_k = x_k - x_km1;
    y_k = grad_k - grad_km1;
    [d, H_k, V, Lambda, opts] = updateHk(k, s_k, y_k, opts);

    % quasi-newton step
    p = -H_k * grad_k;
    if isfield(opts, 'Q')
        % For QP, this is the optimal step length (see page 56 of N&W)
        kappa = -dot(p, grad_k)/ dot(p, opts.Q*p);
        c1 = 1e-4;
        c2 = 0.9;
        [f_test, grad_test] = fcnGrad(x_k + kappa * p);
        assert( f_test <= f_k + c1*kappa*dot(grad_k, p), 'Sufficient decrease condition not satisfied');
        assert( dot(grad_test, p) >= c2 *dot(grad_k, p), 'Curvature condition not satisfied');
    end
    xbar = x_k + kappa * p;
    % if kappa <= opts.kappa / 2^5
    %     kappa = opts.kappa;
    % end
    % % armijo linesearch TODO remove
    % for t = 1:10 
    %     xbar = x_k + kappa* p;
    %     [fnew, ~] = fcnGrad(xbar);
    %     % fpred = f_k + dot(p, grad_k + (diag(1 ./ d) + V * diag(Lambda) * V')* p);
    %     if fnew > f_k %|| fnew / fpred > 2
    %         kappa = kappa / 2;
    %     else 
    %         break
    %     end
    % end
    % proximal step
    x_km1 = x_k;
    x_k = prox(xbar, d, V, Lambda);
    grad_km1 = grad_k;
    [f_k, grad_k] = fcnGrad(x_k);
end

warning('Maximum iterations reached')
x = x_k;
iter = opts.max_iter;

end % sr1_custom

function opts = defaultOpts(opts, N)

    if ~isfield(opts, 'max_iter')
        opts.max_iter = 100;
    end

    if ~isfield(opts, 'tol_rel')
        opts.tol_rel = 1e-4;
    end

    if ~isfield(opts, 'tol_abs')
        opts.tol_abs = 1e-7;
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
    end
    opts.S = zeros(N, opts.r);
    opts.Y = zeros(N, opts.r);
end % defaultOpts

function [d, H, V, Lambda, opts] = updateHk(k, s_k, y_k, opts)

assert(0 < opts.tau_min, "tau_min must be positive");
assert(opts.tau_min < opts.tau_max, "tau_max must be larger than tau_min");
assert(0 < opts.gamma && opts.gamma < 1, "gamma must be in (0,1)");
N= length(s_k);
tau_bb2 = dot(s_k,y_k) / norm(y_k,2)^2;
tau_bb2 = clip(tau_bb2, opts.tau_min, opts.tau_max);
if tau_bb2 == opts.tau_min
    warning('Convexity of cost function is stagnating'); 
end
d = opts.gamma * tau_bb2 * ones(N,1);

r = min(k, opts.r);
if k > r 
    opts.S(:,1:r-1) = opts.S(:,2:r);
    opts.Y(:,1:r-1) = opts.Y(:,2:r);
end
opts.S(:,r) = s_k;
opts.Y(:,r) = y_k;

S = opts.S(:,1:r);
Y = opts.Y(:,1:r);
% ApD = zeros(r,r);
% for i = 1:r
%     for j = 1:r
%         ApD(i,j) = dot(S(:,i), Y(:,j));  
%     end
% end
% 
% BB = chol(S'*diag(1 ./ d)*S - ApD); % for B
% HH = chol(ApD' - Y'*diag(d)*Y); % for H
% V = (Y - diag(1 ./ d)*S) / BB;
% U = (S - diag(d)*Y) / HH;

%% BFGS 
% U = zeros(N, 2*r);
% Sigma = zeros(2*r,1);
V = zeros(N, 2*r);
Lambda = zeros(2*r,1);
H = diag(d);
B = diag(1./d);
for j = 1:r
    s = S(:,j); 
    y = Y(:,j); 
    rho = 1 / dot(y,s);
    assert(1/rho > 1e-8);
    U = (eye(N) - rho*y*s');
    H = U'*H*U + rho*s*s';

    a = B*s / dot(s, B*s);
    b = y / dot(y,s);
    
    B = B - a*a' + b*b';
    Lambda(2*(j-1)+1) = -1;
    Lambda(2*(j-1)+2) = 1;
    V(:, 2*(j-1)+1) = a;
    V(:, 2*(j-1)+2) = b;
end


%% SR1
% idiot check 2 
% Hfast = H; 
% Bfast = B;
% U = zeros(N, r);
% Sigma = zeros(r,1);
% V = zeros(N, r);
% Lambda = zeros(r,1);
% H = diag(d);
% B = diag(1./d);
% for j = 1:r 
%     s = S(:,j); 
%     y = Y(:,j);
%     % update H (and mem)
%     denomH = dot(s - H*y,y);
%     sigma = sign(denomH);
%     denomH = sqrt(sigma*denomH);
%     assert(isreal(denomH));
%     u = (s - H*y) / denomH;
%     H = H + sigma*(u*u');
%     % update H (and mem)
%     denomB = dot(y - B*s,s);
%     lambda = sign(denomB);
%     denomB = sqrt(lambda*denomB);
%     assert(isreal(denomB));
%     v = (y - B*s) / denomB;
%     B = B +lambda*(v*v');
%     % Store for proof of concept
%     U(:,j) = u;
%     Sigma(j) = sigma;
%     V(:,j) = v;
%     Lambda(j) = lambda;
% end
% assert(norm(B*H - eye(N)) < 1e-8)
% assert(norm(B - Bfast) < 1e-9*norm(B))
% assert(norm(H - Hfast) < 1e-9*norm(H))

% idiot check 1
% H = diag(d) + U*diag(Sigma)*U';
% B = diag(1 ./ d) + V*diag(Lambda)*V';
% N = length(s_k);
% assert(norm(B*H - eye(N)) < 1e-8)

end % updateHk

function x = prox(xbar, d, V, Lambda)
    if all(xbar >= 0)
        x = xbar;
        return 
    end
    N = length(xbar); %#ok<NASGU>
    B = diag(1./d)+V*diag(Lambda)*V';
    cvx_begin quiet
        variable z(N)
        minimize( dot(xbar-z,B*(xbar-z)) )
        subject to 
            z >= 0 %#ok<NOPRT>
    cvx_end
    x = z;

    % TODO: solve with grad descent for a couple of iters then newton.
    % cvx_begin
    %     variable a(r)
    %     minimize( sum_square( -V'* V*a * d(1) + a ) )
    %     subject to 
    %         0 <= xbar + V*a * d(1)
    % cvx_end 
    % x = max(0, xbar + (V*a) ./ d);
end% prox_rankr

function [converged, errStruct] = checkConvergence(k, f_k, x_k, ...
    grad_k, errStruct, opts)
% kkt conditions / LCP being satisified is equivalen to
% the grad_k(i) = 0 perp  x_k(i) = 0
phi = min(grad_k,x_k);
old_err = errStruct.err(k);
err = 0.5*(phi'*phi);
errStruct.err(k+1) = err;
errStruct.f(k+1) = f_k;
errStruct.xk(:,k+1) = x_k;
converged = false;
if k == 1
    return 
end

if (abs(err - old_err) / abs(old_err)) < opts.tol_rel  % Relative stopping criteria
    errStruct.reason = "relative";
     converged = true;
elseif err < opts.tol_abs   % Absolute stopping criteria
    errStruct.reason = "absolute";
    converged = true;
end
if converged 
     errStruct.f = errStruct.f(1:k+1);
    errStruct.err = errStruct.err(1:k+1);
    errStruct.xk = errStruct.xk(:,1:k+1);
end

end % checkConvergence