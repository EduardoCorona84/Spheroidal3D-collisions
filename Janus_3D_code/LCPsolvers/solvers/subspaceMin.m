function [x, info] = subspaceMin(fg, x0, opts, debug)
if ~exist('debug','var') || isempty(debug)
    debug = false;
end
[opts, info] = defaultLCPOpts(opts, x0);
% Check that options are set properly for this solver
checkOpts(opts);
% Reset memory to blank
updateSubSpace([],[],opts); 
%% Do one step of Projected Gradient Descent to Get things going
x_k = x0;
n = numel(x0);
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
s=[]; y=[]; x_km1=[]; grad_km1=[];
%% Start the subspace minimization
k = 0;
X = zeros(opts.n, opts.subspaceMin.m);
AX = zeros(opts.n, opts.subspaceMin.m);
while true
    if debug 
        fprintf('f_k = %.4g', f_k);
    end
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, eta, info, opts, x_km1, grad_km1);
    if converged
        x = x_k;
        return
    end
    %% The next search direction should be a descent direction from where we currently are...
    k = k + 1;
    x_km1 = x_k; Ax_km1 = Ax_k; f_km1 = f_k; grad_km1 = grad_k; 
    %% Choose next step direction
    [H, opts] = updateHk(s, y, opts);
    % Define this step 
    q = -H(grad_km1);
    prox_k = @(xtilde) prox(xtilde, opts);
    step_k = @(t, opts) fwdBwdstep(t, x_km1, Ax_km1, q, prox_k, fg, opts);
    % Linesearch 
    [x_k, Ax_k, ~, ~, eta, p, Ap] = linesearch(x_km1, Ax_km1, f_km1, grad_km1, ...
        step_k, opts);
    [P_k, AP_k] = updateSubSpace(p, Ap, opts);
    i = size(P_k,2);
    if k > opts.subspaceMin.m
        xprime = X(:,1);
        Axprime = AX(:,1);
        X(:,1:end-1) = X(:,2:end);
        AX(:,1:end-1) = AX(:,2:end);
    end
    if k > 1 
        % Solve subproblem 
        eta = solveSubProblem(k, xprime, Axprime, P_k, AP_k, opts);
        x_k = xprime+P_k*eta;
        Ax_k = Axprime + AP_k*eta;
    end
    %% Check for convergence
    [f_k, grad_k] = fg(x_k, Ax_k);
    s = x_k - x_km1;
    y = grad_k - grad_km1;
    X(:,i) = x_k;
    AX(:,i) = Ax_k;
end

end % subspaceMin

function [Pk, APk] = updateSubSpace(p, Ap, opts)
persistent P AP 
Pk = [];
APk = [];
n = opts.n;
m = opts.subspaceMin.m;
if ~exist('p','var') || isempty(p)
    % reset memory
    P = zeros(n,m);
    AP = zeros(n,m);
    return 
end

i = find(arrayfun(@(i) all(P(:,i) == 0), 1:m), 1, 'first');
if isempty(i)
    P(:, 1:m-1) = P(:,2:m);
    AP(:,1:m-1) = AP(:,2:m);
    i = m;
end
P(:,i) = p;
AP(:,i) = Ap;
[Pk,R] = qr(P(:,1:i),0); % 'econ' mode in new syntax
APk = AP(:,1:i) / R;

B = Pk'*APk; 
Bt = B';
if norm(B-Bt) / norm(B) > 1e-6 || cond(B) > 1e12 || any(isnan(B(:)))
    warning('We are losing condition of the problem.')
    APk = opts.AA*Pk;
end

end

function z = solveSubProblem(k, xprime, Axprime, P, AP, opts)
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
d = -xprime;
c = P'*(b + Axprime); 
k = size(P, 2); %#ok<NASGU>
if cond(B) > 1e12
    warning('Conditioning of B is bad')
elseif cond(P) > 1e12 
    warning('Conditioning of P is bad')
end
cvx_begin quiet
    % cvx_precision best
    variable z(k)
    dual variable u
    minimize dot(0.5*B*z + c, z)
    subject to 
    u : d <= P*z   %#ok<NODEF,NOPRT>
cvx_end 

if debug
    dualToPrimalMap = @(u) B \ (P'*u-c);
    disp('Checking the relErr of the dualToPrimalMap')
    disp(norm(dualToPrimalMap(u) - z) / norm(z))
    n = numel(xprime); %#ok<NASGU>
    C = P*(B\P'); %#ok<NASGU>
    r = P*(B\c) + d; %#ok<NASGU>
    cvx_begin quiet
        variable u1(n)
        minimize dot(0.5*C*u1-r, u1)
        subject to
        0 <= u1  %#ok<NOPRT>
    cvx_end;
    primalLoss = @(z) dot(0.5*B*z + c, z);
    disp("primal loss via cvx solving the dual problem ")
    disp(primalLoss(dualToPrimalMap(u1)))
    disp("primal loss via cvx solving the primal problem ")
    disp(primalLoss(z))
    if primalLoss(dualToPrimalMap(u1)) < primalLoss(z)
        disp("dual is better")
    else
        disp('primal is better')
    end
    disp('relPrimal Loss Gap')
    disp(abs( (primalLoss(dualToPrimalMap(u)) - primalLoss(dualToPrimalMap(u1))) / primalLoss(z)))
    disp("primal rel argmin gap ")
    disp(norm(dualToPrimalMap(u1) - z) / norm(z))
end
end % solveSubProblemCVX

function checkOpts(opts) %#ok<INUSD>
    % assert(isfield(opts, 'subspaceMin'), 'Necessary to have specification for the subspaceMin algo')
    % assert(isfield(opts.subspaceMin, 'innerSolver'), 'Necessary to have specify for the innerSolver')
end
