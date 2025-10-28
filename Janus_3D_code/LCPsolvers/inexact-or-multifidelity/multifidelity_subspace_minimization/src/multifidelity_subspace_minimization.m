function [x, info] = multifidelity_subspace_minimization(fg, x0, opts)
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
    %{
    if opts.subSpaceMin.useOneStepIter
        x_km1 = xhat_k;
    end
    %}
    %take one high-fidelity step and record the pre-projection step direction
    [x_k_half, eta, p, Ap] = high_step_algo(k, x_km1, Ax_km1, grad_km1, s, y, opts);
    %take one low-fidelity step and record the high-fidelity free low fidelity step direction
    [] = low_step_algo(k, x_k, Ax_k, )
    %I don't understand the point of burn in
    %{
    if k <= opts.subSpaceMin.kThresh
        % burnIn
        Ax_k = Ax_km1 + eta*Ap;
        xprime = x_k;
        Axprime = Ax_k;
    %}
    
        % Update subspace matrix
    [P, AP] = updateSubSpace(k, p, Ap, opts);
    % Solve subproblem 
    z_old = [z_km1; eta];
    z_k = solveSubProblem(k, xprime, Axprime, P, AP, z_old, opts);
    % update x1
    x_k = xprime+P*z_k;
    Ax_k = Axprime + AP*z_k;
    %% Check for convergence
    [f_k, grad_k] = fg(x_k, Ax_k);
    s = x_k - x_km1;
    y = grad_k - grad_km1;
end

end % subspaceMin

function [x_k_half, search_direction] = high_step_algo(k, x_km1, Ax_km1, grad_km1, s, y, opts)
innerStep = lower(opts.subSpaceMin.innerStepSelection);
switch innerStep
    case 'pgd'
        q = -grad_km1; % use gradient descent?
        [kappa, ~] = stepSize(k, q, x_km1, Ax_km1, opts, s, y);
        x_k_half = x_km1 + kappa * q;
        % Projection to positive orthant
        x_k_half = max(0, x_k_half);
    case 'pqn'
        [H, h0, U, V] = updateHk(k, s, y, opts);
        % quasi-newton step direction
        q = -H(grad_km1);
        % step size direction
        kappa = stepSize(k, q, x_km1, Ax_km1, opts);
        x_k = x_km1 + kappa * q;
        x_k = prox(x_k, h0, U, V, opts);
    otherwise
        error([innerStep ' not implemented'])
end
% Possibly a step length update after the projection
p = x_k - x_km1;
p = p / norm(p);
[eta, Ap] = stepSize(-1, p, x_km1, Ax_km1, opts);
x_k = x_km1 + eta*p; 
end % oneStepAlgo

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
k = k - opts.subSpaceMin.kThresh;
if isempty(P) || isempty(AP)
    assert(k <= 2, 'This should only happen on the first iteration');
    n = length(p);
    P = zeros(n,n);
    AP = zeros(n,n);
end
orthoMethod = lower(opts.subSpaceMin.orthoMethod);
switch orthoMethod
case 'cgs'
    % Make orthogonal to all other vectors in the subspace (build orthonormal P)
    % Here we are using graham-schmidt, which may be UNSTABLE, so we should update this later...
    p = p / norm(p);
    for j = 1:k-1
        pj = P(:,j);
        p = p - dot(pj, p) / dot(pj,pj) * pj; 
        p = p / norm(p);
    end
    P(:,k) = p;
    AP(:,k) = Ap;
    %% Form inner optimization problem 
    Pk = P(:,1:k);
    APk = AP(:,1:k);
case 'qr'
    P(:,k) = p;
    AP(:,k) = Ap;
    [Pk,R] = qr(P(:,1:k),0);%'econ');
    APk = AP(:,1:k) / R;
end


% B = Pk'*APk; 
% Bt = B';
% if norm(B-Bt) / norm(B) > 1e-6
%     warning('We are losing condition of the problem.')
%     APk = opts.AA*Pk;
% end

end

function z = solveSubProblem(k, xprime, Axprime, P, AP, z_old, opts)
k = k - opts.subSpaceMin.kThresh;
if k == 1 
    % The one dimentional optimal solution is found by the linesearch algo
    z = z_old;
    return 
end
innersolvername = lower(opts.subSpaceMin.innerSolver);
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

function z = solveSubProblemCVX(xprime, Axprime, P, AP, b, debug) %#ok<STOUT>
if ~exist('debug','var') || isempty(debug)
    debug = false;
end
B = P'*AP; 
B = 0.5*(B+B');
c = P'*(b + Axprime); 
d = -xprime;
k = size(P, 2); %#ok<NASGU>
cvx_begin quiet
    % cvx_precision best
    variable z(k)
    dual variable u
    minimize dot(0.5*B*z + c, z)
    subject to 
    u : d <= P*z  %#ok<NODEF,NOPRT>
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
    assert(isfield(opts, 'subSpaceMin'), 'Necessary to have specification for the subSpaceMin algo')
    assert(isfield(opts.subSpaceMin, 'innerSolver'), 'Necessary to have specify for the innerSolver')
end
