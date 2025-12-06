function xhat = prox_rankr(xtilde, opts)
%% Return simple answer in the trivial cases
if all(xtilde >= 0) 
    xhat = xtilde;
    return
end
%% Unpack parameters
maxiter = opts.prox.maxiter;
res_abstol = opts.prox.res_abstol;
res_reltol = opts.prox.res_reltol;
alp_abstol = opts.prox.alp_abstol;
alp_reltol = opts.prox.alp_reltol;
prox_B0 = opts.prox.prox_B0;
B0 = opts.prox.B0;
H0 = opts.prox.H0;
U = opts.prox.U;
V = opts.prox.V;
verbose = opts.prox.verbose;
runCVX = opts.prox.runCVX;
%% Base case of no low rank update
if isempty(U) && isempty(V)
    xhat = prox_B0(xtilde);
    return 
end
n = length(xtilde);
r1 = size(U,2);
r2 = size(V,2);
r = r1 + r2;
%% Make as efficient as possible by using matrix free implementations 
if isempty(V)
    Utilde = -H0(U);
    xa = @(a) prox_B0(xtilde + Utilde * a);
    L = @(a,xa) U' * (xtilde - xa) + a;
    % This is a bad approximation in the non-diagonal B0 case...
    Lambda = @(xa) diag(abs(xa) >= 1e-8);
    J_L = @(a,xa) U'*Lambda(xa)*H0(U) + eye(r);
else
    C = B0(eye(n,n)) + U*U';
    R = chol(C);
    Cinv = @(x) R\(R'\x);
    Utilde = cat(2, -H0(U) , Cinv(V));
    xa = @(a) prox_B0(xtilde + Utilde * a);
    L = @(a,xa) cat(1, ...
        U' * (xtilde + Cinv(V*a(r1+1:end,:)) - xa), ...
        V' * (xtilde - xa)...
        ) + a;
    % This is a bad approximation in the non-diagonal B0 case...
    Lambda = @(xa) diag(abs(xa) >= 1e-8);
    J_L = @(a,xa) cat(1, ...
        cat(2, U'*Lambda(xa)*H0(U), U'*(eye(n) - Lambda(xa))*Cinv(V)), ...
        cat(2, V'*Lambda(xa)*H0(U), -V'*Lambda(xa)*Cinv(V)) ...
        ) + eye(r);
end
%% Start semi-smooth newton iterations
k = 0;
a_km1 = zeros(r,1); % a0
xa_km1 = xa(a_km1);
L_km1 = L(a_km1, xa_km1);
while true
    k = k + 1;
    J_L_km1 = J_L(a_km1, xa_km1);
    p = -J_L_km1\L_km1;
    a_k = a_km1 + p;
    alp_abserr = norm(a_k - a_km1);
    alp_relerr = alp_abserr / norm(a_k);
    if alp_abserr < alp_abstol 
        if verbose 
            disp('Iterate absolute error is below tolerence')
        end
        break
    elseif alp_relerr < alp_reltol 
        if verbose 
            disp('Iterate relative error is below tolerence')
        end
        break 
    end
    xa_k = xa(a_k);
    L_k = L(a_k,xa_k);
    res_abserr = norm(L_k);
    res_relerr =  norm(L_k - L_km1) / norm(L_k);
    if res_abserr < res_abstol 
        if verbose 
            disp('Residiual absolute error is below tolerence')
        end
        break 
    elseif res_relerr < res_reltol 
        if verbose 
            disp('Residiual relative error is below tolerence')
        end
        break 
    elseif k == maxiter
        if verbose 
            disp('Maximum Iterations hit... did not converge')
            disp(['  Iterate Relative error : ' alp_relerr])
        end
        break
    end
    a_km1 = a_k;
    xa_km1 = xa_k;
    L_km1 = L_k;
end
%% Now that we have the root of L, we apply use it to obtain x(alphastar)
xhat = xa(a_k);

if runCVX 
    %% Debug with CVX if necessary
    B = B0(eye(n)) + U*U' - V*V';
    % cvx
    tic
    f = @(x)1/2*dot(x-xtilde, B*(x-xtilde)) ; %#ok<NASGU>
    cvx_begin quiet
            variable xRef(n)
            minimize f(xRef)
            subject to 
            0 <= xRef %#ok<NOPRT>
    cvx_end 
    cvxErr = norm(xRef - xhat) / norm(xRef);
    assert(cvxErr < 1e-4, 'Compared to CVX we have a bad answer')
end

end % prox