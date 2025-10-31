function [x,info] = manifoldIdentificationExperiment
n = 10;
max_iter = 1000;
tol = 1e-6;
condNum = 1e2;
% rng(1)
%% Set paths
setPaths();
%% Try to initialize cvx if on the cluster
try %#ok<TRYNC>
    run(fullfile('/projects', getenv('USER'), 'cvx', 'cvx_startup.m'))
end
%% Preallocate structures
opts = struct( ...
    'max_iter',max_iter, ...
    'kkt_rel',tol, ...
    'kkt_abs',tol, ...
    'storeIts', true, ...
    'subspaceMin', struct( ...
    'innerStepSelection', 'pqn', ...
    'innerSolver','cvx', ...
    'orthoMethod','qr', ...
    'kThresh', 0, ...
    'useOneStepIter', false) ...
);
%% Generate random problem
B = randn(n,n);
vv = [1; rand(n-2,1)*condNum;condNum];
[Q, ~] = qr(B);
A = Q*diag(vv)*Q';
[vecs,vals] = eig(A);
vals = max(1,diag(vals));
A = vecs*diag(vals)*vecs';
x_unconstrained = randn(n,1);
while ~any(x_unconstrained<0) || all(x_unconstrained < 1)
    x_unconstrained = randn(n,1);
end
b = -A*x_unconstrained;

%% Build cost function
Acnt = @(x) Acounter(x,A, false);
x0 = zeros(n,1);
fg = @(x, Ax) quadraticLoss(x, Acnt, b, Ax);
%% Fill the opts with problem specific information
opts.A = Acnt;
opts.b = b;
xlb = A \ -b;
opts.fstar = 1/2*dot(xlb, A*xlb) + dot(xlb, b);
evals = eig(A);
opts.L = max(evals);
opts.mu = min(evals);
opts.AA = A;

[xstar,~] = callCVX([], x0, opts);
opts.metricNames = {'abs_kkt', 'rel_kkt', 'MVP', 'rel_iter', 'abs_iter', 'obj'};
opts.errFcn = {
    @(x) abs_kkt(x, A, b);
    @(x) rel_kkt(x, A, b);
    @(x) Acnt('cnt');
    @(x) rel_iter(x,xstar);
    @(x) abs_iter(x,xstar);
    @(x) dot(x,0.5*A*x+b);
    };
opts.solver = 'subspaceMin';
[x, info] = subspaceMin(fg, x0, opts);
if ~all(x > -tol)
    disp('Constraints are being violated too much')
end

if fg(x, []) < fg(xstar, []) 
    xstar = x;
end
tol = 1e-3;
K = size(info.iterHist,1);
I_star = (abs(xstar) < tol)';
foundActiveSet = zeros(K,1);
for k = 1:K
    xk = info.iterHist(k,:);
    I_Ak = abs(xk) < tol;
    if all(I_star == I_Ak)
        foundActiveSet(k) = 2;
    elseif all(I_Ak(I_star))
        foundActiveSet(k) = 1;
    end
end 
gcf;
clf; 
plot(1:k, foundActiveSet, '-','Marker','.', 'MarkerSize',30)
ylim([-1,3])
yticks([0,1,2])
xticks(1:k)
xlabel('Iteration k')
yticklabels({'$\mathcal{A}(x^*) \not\subseteq \mathcal{A}(x^{(k)})$', '$\mathcal{A}(x^*) \subset \mathcal{A}(x^{(k)})$', '$\mathcal{A}(x^*) = \mathcal{A}(x^{(k)})$'})
set(gca, 'TickLabelInterpreter', 'latex')
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
title('Manifold Identification For subspaceMin','FontSize',30)
disp('Iterates [x_1, ..., x_K, x^*]^T')
disp([info.iterHist;
xstar'])
end

function Ax = Acounter(x, A, transpose)
persistent matVecCnt

if isempty(matVecCnt)
    matVecCnt = 0;
end
if ~exist('transpose','var') || isempty(transpose)
    transpose = false;
end

if ischar(x)
    if strcmpi(x, 'cnt')
        Ax = matVecCnt;
        return
    elseif strcmpi(x, 'reset')
        Ax = matVecCnt;
        matVecCnt = 0;
        return
    else
        assert(false, ['Option: ' x ' not recognized'])
    end
end
if transpose
    Ax = A'*x;
else
    Ax = A*x;
end
matVecCnt = matVecCnt + 1;
end

function [x, info] = callCVX(~, x0, opts) %#ok<STOUT>
n = length(x0); %#ok<NASGU>
% f = @(x) 1/2*sum_square(opts.A(x)+opts.b);
f = @(x) 1/2*dot(x,opts.A(x)) + dot(opts.b,x); %#ok<NASGU>
cvx_begin quiet
variable x(n)
minimize f(x)
subject to
0 <= x %#ok<NODEF,NOPRT>
cvx_end

info.iter = NaN;
info.kkt = NaN;
info.errHist = NaN;
info.iterHist = NaN;
end

function e = abs_kkt(x, A, b)
phi = min(x,A*x + b);
e = 1/2*dot(phi, phi);
end % abs_kkt

function diff = rel_kkt(x, A, b)
persistent olde
if ischar(x) && strcmpi(x,'reset')
    olde = [];
    diff = NaN;
    return
end
phi = min(x,A*x + b);
e = 1/2*dot(phi, phi);
if isempty(olde)
    diff = NaN;
else
    diff = abs(e - olde) / abs(e);
end
olde = e;
end % rel_kkt

function e = abs_iter(x, xstar)
e = norm(x - xstar);
end % abs_iter

function e = rel_iter(x, xstar)
e = norm(x-xstar) / norm(xstar);
end % rel_iter
