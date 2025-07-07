function xstar = prox_rankr(y, d, U, V, opts)
if all(y >= 0) 
    xstar = y;
    return
elseif isempty(U) && isempty(V)
    xstar = max(y,0);
    return 
end
% TODO better handling of defaults    
if ~exist('opts','var') || isempty(opts) || ~isfield(opts, 'prox')
    opts = struct('prox', struct( ...
        'maxiter', 100, ...
        'tol', 1e-12, ...
        'newton', true, ...
        'gradDescentWarmStartIter', 0, ...
        'c1', 1e-4, ...
        'c2', .9, ...
        'debug', false));
end
% parameters
maxiter = opts.prox.maxiter;
tol = opts.prox.tol;
newton = opts.prox.newton;
gradDescentWarmStartIter = opts.prox.gradDescentWarmStartIter;
c1 = opts.prox.c1;
c2 = opts.prox.c2;
debug = opts.prox.debug;
n = length(y);
r1 = size(U,2);
r2 = size(V,2);
r = r1 + r2;
% 
B0 = @(x) x ./ d;
H0 = @(x) x .* d;
Dsqrt = @(x) sqrt(d) .* x ;
% Disqrt = @(x) x ./ sqrt(d);
B1 = @(x) B0(x) + U*(U'*x);
B1_ = B0(eye(n,n)) + U*U';
% assert(min(eig(B1_)) > 0);
% B = B1_ - V*V'; 
% eigVals = eig(B);
% assert(min(eigVals) > 0);
R = chol(B1_);
B1inv = @(x) R\(R'\x);
B1sqrt = @(x) (R'*x);
% B1sqrti = @(x) (R'\x);
% Dtilde = @(x) cat(1, B0(x(1:r1,:)), B1(x(r1+1:end,:)));
% Dtildei = @(x) cat(1, H0(x(1:r1,:)), B1inv(x(r1+1:end,:)));
% Dtildesqrt = @(x) cat(1, Dsqrt(x(1:r1,:)), B1sqrt(x(r1+1:end,:)));
% Dtildesqrti = @(x) cat(1, Dsqrti(x(1:r1,:)), B1sqrti(x(r1+1:end,:)));
Atilde = @(a) cat(1, U' *B0(U*a(1:r1,:)), V'*B1(V*a(r1+1:end,:)));

%% Define functions
Utilde = cat(2, -H0(U) , B1inv(V));
xa = @(a) max(0, y + Utilde * a);
g = @(a) 1/2*dot(a, a + Atilde(a)) ...
    - 1/2*sum(cat(1, ...
    Dsqrt(xa(a)) - y - B1inv(V*a(r1+1:end,:)), ...
    B1sqrt(xa(a)) - y).^2);
grad_g = @(a) cat(1, ...
    U' * (y + B1inv(V*a(r1+1:end,:)) - xa(a)), ...
    V' * (y - xa(a))...
    ) + a;
Lambda = @(a) diag(sign(xa(a)));
hess_g = @(a) cat(1, ...
    cat(2, U'*Lambda(a)*H0(U), U'*(eye(n) - Lambda(a))*B1inv(V)), ...
    cat(2, V'*Lambda(a)*H0(U), -V'*Lambda(a)*B1inv(V)) ...
    ) + eye(r);
% 
a0 = zeros(r,1);
ak = a0;
akp1 = a0;
for iter = 1:maxiter
    ak = akp1; 
    gradk = grad_g(ak);
    hessk = hess_g(ak);
    if newton && iter >= gradDescentWarmStartIter
        eta = 1;
        pk = -hessk\gradk;
    else
        pk = -gradk;
        eta = 1;%pk'*gradk / (pk'*hessk*pk);
    end
    akp1 = ak + eta*pk;
    gradkp1 = grad_g(akp1);
    while ( g(akp1) > g(ak) + c1 * eta * dot(pk, gradk)... % sufficient decrease 
        && dot(gradkp1,pk) < c2*dot(gradk, pk) )
        eta = eta*.5;
        akp1 = ak + eta*pk;
        gradkp1 = grad_g(akp1);
    end
    relErr = norm(ak - akp1) / norm(ak);
    if relErr < tol 
        if debug 
            disp('Converged')
        end
        xstar = xa(akp1);
        break
    end 
end

if iter == maxiter
    if debug 
        disp('Maximum Iterations hit... did not converge')
        disp(['  Relative Error is ' relErr])
    end
    xstar = xa(akp1);
end



if debug 
    B = diag(1./d) + U*U' - V*V';
    % cvx
    tic
    f = @(x)1/2*dot(x-y, B*(x-y)) ;
    cvx_begin quiet
            variable xRef(n)
            minimize f(xRef)
            subject to 
            0 <= xRef
    cvx_end 
    err = norm(xRef - xstar) / norm(xRef);
    assert(err < 1e-4)
end

end % prox