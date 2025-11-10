%% Set up problem 
rng(1);
n = 1000; % size of the original problem 
r = 40; % rank of the update
d = 5;
D = @(x) d .* x;
Dinv = @(x) x ./ d;
Dsqrt = @(x) sqrt(d) .* x ;
Disqrt = @(x) x./ sqrt(d);
U = rand(n,r);
B = D(eye(n,n)) + U*U'; 
y = randn(n,1);
%% Define functions
xa = @(a) max(0, y  - Dinv(U * a));
g = @(a) -1/2* sum((Dsqrt(xa(a)) - Dsqrt(y) + Disqrt(U*a)).^2) ...
    + 1/2*(a'*a + Dinv(dot(U*a, U*a)));
grad_g = @(a) a + U'*(y - xa(a));
I_r = eye(r,r);
hess_g = @(a) I_r + Dinv(U'*diag(sign(xa(a)))*U);
%% Test it out
a0 = zeros(r,1);
ak = a0;
akp1 = a0;
maxiter = 100;
tol = 1e-6;
newton = true;
c1 = 1e-4;
c2 = .9;
tic;
for iter = 1:maxiter
    ak = akp1; 
    gradk = grad_g(ak);
    hessk = hess_g(ak);
    if newton && iter > 5
        eta = 1;
        pk = -hessk\gradk;
    else
        pk = -gradk;
        eta = pk'*gradk / (pk'*hessk*pk);
    end
    akp1 = ak + eta*pk;
    gradkp1 = grad_g(akp1);
    while (g(akp1) > g(ak) + c1 * eta * dot(pk, gradk)... % sufficient decrease 
        && dot(gradkp1,pk) < c2*dot(gradk, pk)...
        )
        % disp('here')
        eta = eta*.5;
        akp1 = ak + eta*pk;
        gradkp1 = grad_g(akp1);
    end
    if iter > 5
        relErr = norm(ak - akp1) / norm(ak);
        if relErr < tol 
            disp('Converged')
            break 
        end 
    end
end
if iter == maxiter 
    disp('Maximum Iterations hit... did not converge')
    disp(['  Relative Error is ' relErr])
end

runTime = toc;
xstar = xa(ak);
tic
cvx_begin quiet
        variable xRef(n)
        minimize 1/2*dot(xRef-y, B*(xRef-y)) 
        subject to 
        0 <= xRef
cvx_end 
cvxTime = toc;
disp(['Problem Size: n = ' num2str(n) ', r = ' num2str(r)])
if newton
    disp('Semi-smooth newton')
else 
    disp('Gradient Descent')
end
disp(['  iter = '  num2str(iter) ', run time ' num2str(runTime) 's'])
disp('  xstar')
disp(['    ' num2str(xstar')])
disp('CVX')
disp(['  iter = Not Avalailable, run time ' num2str(cvxTime) 's'])
disp('  xRef')
disp(['    ' num2str(xRef')])
disp('Relative Error')
disp(norm(xRef - xstar) / norm(xRef))