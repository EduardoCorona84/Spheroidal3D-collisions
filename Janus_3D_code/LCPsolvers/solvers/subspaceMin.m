function [x, info] = subspaceMin(fg, x0, opts)
[opts, info] = defaultLCPOpts(opts, x0);
% Check that options are set properly for this solver
checkOpts(opts);
% Reset memory to blank
updateSubSpace(); 
%% Do one step of Projected Gradient Descent to Get things going
x_k = x0;
if all(x0 == 0)
    n = numel(x0);
    Ax_k = zeros(n,1);
else 
    Ax_k = opts.A(x_k);
end
[f_k, grad_k] = fg(x0,Ax_k);
xprime = x0;
Axprime = Ax_k;
eta = 1;
s=[]; y=[]; z_k=[]; x_km1=[]; grad_km1=[];
xhat_k = x_k;
%% Start the subspace minimization
k = 0;
burnIn = opts.subspaceMin.kThresh;
m = opts.max_iter;
X = zeros(n, opts.max_iter);
AX = zeros(n, opts.max_iter);
Eta = cell(opts.max_iter,1);
P = cell(opts.max_iter,1);
while true
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts, x_km1, grad_km1);
    if converged
        x = x_k;
        return
    end
    %% The next search direction should be a descent direction from where we currently are...
    k = k + 1;
    x_km1 = x_k; Ax_km1 = Ax_k; grad_km1 = grad_k; z_km1 = z_k;
    if opts.subspaceMin.useOneStepIter
        x_km1 = xhat_k;
    end
    %% Choose next step direction
    [H, h0, U, V] = updateHk(k, s, y, opts);
    % quasi-newton step direction
    q = -H(grad_km1); 
    % step size direction
    kappa = stepSize(k, q, x_km1, Ax_km1, opts); % Should always be 1... 
    y_k = x_km1 + kappa * q;
    xhat_k = prox(y_k, h0, U, V, opts);
    p = xhat_k - x_km1;
    [~, Ap] = stepSize(-k, p, x_km1, Ax_km1, opts);
    % if k <= 1
    %     Eta_k = eta; 
    %     P_k = p;
    %     x_k = x_km1 + eta*p;
    %     Ax_k = Ax_km1 + eta*Ap;
    %     xprime = x_k;
    %     Axprime = Ax_k;
    % else
        % if k - burnIn > m 
        %     i = k-burnIn-m+1;
        %     xprime = X(:, i);
        %     Axprime = AX(:, i);
        % end
        % Update subspace matrix
        [P_k, AP_k] = updateSubSpace(k, p, Ap, opts);
        % Solve subproblem 
        Eta_k = solveSubProblem(k, xprime, Axprime, P_k, AP_k, opts);
        % update x1
        x_k = xprime+P_k*Eta_k;
        Ax_k = Axprime + AP_k*Eta_k;
    % end
    %% Check for convergence
    [f_k, grad_k] = fg(x_k, Ax_k);
    f_k
    s = x_k - x_km1;
    y = grad_k - grad_km1;
    X(:,k) = x_k;
    AX(:,k) = Ax_k;
    Eta{k} = Eta_k;
    P{k} = P_k;
end

end % subspaceMin

function [H, h0, U, V] = updateHk(k, s, y_k, opts)
switch lower(opts.qnUpdate)
    case 'bfgs'
        [H, h0, U, V] = get_H_BFGS(k, s, y_k, opts);
    % case 'sr1'
        
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

function [Pk, APk] = updateSubSpace(k, p, Ap, opts)
persistent P AP 

if ~exist('k','var') 
    % reset memory
    P = [];
    AP = [];
    Pk = [];
    APk = [];
    return 
end
n = numel(p);
m = opts.max_iter;
burnIn = opts.subspaceMin.kThresh;
i = min(m, k - burnIn);
if isempty(P) || isempty(AP)
    assert(k <= 2, 'This should only happen on the first iteration');
    P = zeros(n,m);
    AP = zeros(n,m);
elseif k-burnIn > m
    P = [P(:, 2:m) zeros(n,1)];
    AP = [AP(:, 2:m) zeros(n,1)];
end
P(:,i) = p;
AP(:,i) = Ap;
[Pk,R] = qr(P(:,1:i),0);%'econ');
APk = AP(:,1:i) / R;

B = Pk'*APk; 
Bt = B';
if norm(B-Bt) / norm(B) > 1e-6 || cond(B) > 1e12 || any(isnan(B(:)))
    warning('We are losing condition of the problem.')
    APk = opts.AA*Pk;
end

end

function z = solveSubProblem(k, xprime, Axprime, P, AP, opts)
k = k - opts.subspaceMin.kThresh;
if k == 1 
    % The one dimentional optimal solution
    x = xprime;
    b = opts.b;
    p = P;
    Ax = Axprime;
    Ap = AP;
    z_unc = -(dot(p, Ax+b)) / dot(p, Ap);
    if z_unc < 0
        mask = p > 0;
    else 
        mask = p < 0;
    end
    z = min([z_unc; - x(mask) ./ p(mask)]);
    return 
end
innersolvername = lower(opts.subspaceMin.innerSolver);
b = opts.b;
switch innersolvername
    case 'cvx'
        % try 
            z = solveSubProblemCVX(xprime, Axprime, P, AP, b);
        % catch 
        %     z = z_old;
        % end
    case 'pgd'
        z = solveSubProblemPGD(xprime, Axprime, P, AP, b, z_old);
    otherwise
        error([innersolvername ' not recognized as an inner solver'])
end
end

function z = solveSubProblemPGD(xprime, Axprime, P, AP, b, z_old, debug) %#ok<INUSD>
if ~exist('debug','var') || isempty(debug)
    debug = false;
end
%% Precompute some stuff
B = P'*AP; 
B = 0.5*(B+B');
dB = decomposition(B);
c = P'*(b + Axprime); 
d = -xprime;
n = size(P,1); 
C = P*(dB\P');
r = P*(dB\c) + d; 
%% 
% maximize dot(-0.5*C*u,u) + dot(r, u)
% minimize dot(0.5*C*u, u) + dot(-r, u)
fg = @(u, Cu) quadraticLoss(u, @(u) C*u, -r, Cu); 
u0 = zeros(n,1); % Perhaps initialize better?
% Use the default opts
dualToPrimalMap = @(u) B \ (P'*u-c);
primalLoss = @(z) dot(0.5*B*z + c, z);
subOpts = struct( ...
    'max_iter',1000, ...
    'kkt_rel',[], ...
    'kkt_abs',[], ...
    'arg_rel',1e-12, ...
    'arg_abs',1e-12, ...
    'step_abs',1e-12, ...
    'stepSize',struct('init','uniform',...
        'kappa','uniform',...
        'eta','opt'),...
    'A', @(u) C*u, ...
    'b', -r ...
);
[subOpts, ~] = defaultLCPOpts(subOpts,u0);
dualLoss = @(u) dot(0.5*C*u - r, u);
if debug
    subOpts.errFcn = {
        @(u) rel_iter(u), ...
        @(u) primalLoss(dualToPrimalMap(u)) - dot(u, P*dualToPrimalMap(u) - d),...
        @(u) -dualLoss(u)...
    };
    subOpts.errFcn{1}([]); % reset
end
% Solve with PGD
[u, info] = proxQuasiNewton(fg, u0, subOpts );
% apply dual-primal map
z = dualToPrimalMap(u);
if debug
    k = size(P,2); %#ok<NASGU>
    cvx_begin
        variable zCVX(k)
        dual variable w
        minimize dot(0.5*B*zCVX + c, zCVX)
        subject to
        uCVX : d <= P*zCVX  %#ok<NOPRT>
    cvx_end
    disp("primal loss pgd")
    primalLoss(z)
    disp("primal loss cvx")
    primalLoss(zCVX)
    if primalLoss(zCVX) < primalLoss(z)
        disp("CVX is better")
    else
        disp("PGD is better")
    end
    disp('relPrimal Loss Gap')
    abs( (primalLoss(zCVX) - primalLoss(z)) / primalLoss(z))
    disp("primal rel argmin gap ")
    norm(zCVX - z) / norm(zCVX) %#ok<NOPRT>
    figure()
    hold on 
    subfigure(3,1,1)
    semilogy(info.errHist(1,:))
    ylabel('Relative Iterate Error')
    subfigure(3,1,2)
    semilogy(info.errHist(1,:))
    ylabel('Primal Loss')
    subfigure(3,1,3)
    semilogy(info.errHist(1,:))
    ylabel('Dual Loss')
end
end % solveSubProblemPGD

function e = rel_iter(u_k)
persistent u_km1
e = 1e6;
if ~exist('u_k', 'var') || isempty(u_k)
    u_km1 = [];
    return 
elseif isempty(u_km1) 
    u_km1 = u_k;
    return
end
e = norm(u_k-u_km1) / max(norm(u_k) , norm(u_km1));
u_km1 = u_k;
end % rel_iter

function eta = solveSubProblemCVX(xprime, Axprime, P, AP, b, debug) %#ok<STOUT>
if ~exist('debug','var') || isempty(debug)
    debug = false;
end
D = P'*AP; 
D = 0.5*(D+D'); %#ok<NASGU>
d = P'*(b + Axprime); 
k = size(P, 2); %#ok<NASGU>
if cond(D) > 1e12
    warning('Conditioning of D is bad')
elseif cond(P) > 1e12 
    warning('Conditioning of P is bad')
end
cvx_begin quiet
    % cvx_precision best
    variable eta(k)
    dual variable u
    minimize dot(0.5*D*eta + d, eta)
    subject to 
    u : 0 <= P*eta + xprime  %#ok<NODEF,NOPRT>
cvx_end 

if debug
    dualToPrimalMap = @(u) B \ (P'*u-c);
    disp('Checking the relErr of the dualToPrimalMap')
    norm(dualToPrimalMap(w) - z) / norm(z) %#ok<NOPRT>
    n = numel(xprime); %#ok<NASGU>
    C = P*(B\P'); %#ok<NASGU>
    r = P*(B\c) + d; %#ok<NASGU>

    cvx_begin
        variable u(n)
        minimize dot(0.5*C*u-r, u)
        subject to
        0 <= u  %#ok<NOPRT>
    cvx_end
    primalLoss = @(z) dot(0.5*B*z + c, z);
    disp("primal loss via cvx solving the dual problem ")
    primalLoss(dualToPrimalMap(u))
    disp("primal loss via cvx solving the primal problem ")
    primalLoss(z)
    if primalLoss(dualToPrimalMap(u)) < primalLoss(z)
        disp("dual is better")
    else
        disp('primal is better')
    end
    disp('relPrimal Loss Gap')
    abs( (primalLoss(dualToPrimalMap(u)) - primalLoss(dualToPrimalMap(w))) / primalLoss(z))
    disp("primal rel argmin gap ")
    norm(dualToPrimalMap(u) - z) / norm(z) %#ok<NOPRT>
end
end % solveSubProblemCVX

function checkOpts(opts)
    % assert(isfield(opts, 'subspaceMin'), 'Necessary to have specification for the subspaceMin algo')
    % assert(isfield(opts.subspaceMin, 'innerSolver'), 'Necessary to have specify for the innerSolver')
end
