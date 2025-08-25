function xstar = prox_rankr(y, h0, U, V, opts)
%% Return simple answer in the trivial cases
if all(y >= 0) 
    xstar = y;
    return
elseif isempty(U) && isempty(V)
    xstar = max(y,0);
    return 
end
%% Unpack parameters
maxiter = opts.prox.maxiter;
res_abstol = opts.prox.res_abstol;
res_reltol = opts.prox.res_reltol;
alp_abstol = opts.prox.alp_abstol;
alp_reltol = opts.prox.alp_reltol;
verbose = opts.prox.verbose;
runCVX = opts.prox.runCVX;
n = length(y);
r1 = size(U,2);
r2 = size(V,2);
r = r1 + r2;
%% Make as efficient as possible by using matrix free implementations where
% we can
B0 = @(x) x ./ h0;
H0 = @(x) x .* h0;
C = B0(eye(n,n)) + U*U';
R = chol(C);
B1inv = @(x) R\(R'\x);
Utilde = cat(2, -H0(U) , B1inv(V));
xa = @(a) max(0, y + Utilde * a);
L = @(a) cat(1, ...
    U' * (y + B1inv(V*a(r1+1:end,:)) - xa(a)), ...
    V' * (y - xa(a))...
    ) + a;
Lambda = @(a) diag(sign(xa(a)));
J_L = @(a) cat(1, ...
    cat(2, U'*Lambda(a)*H0(U), U'*(eye(n) - Lambda(a))*B1inv(V)), ...
    cat(2, V'*Lambda(a)*H0(U), -V'*Lambda(a)*B1inv(V)) ...
    ) + eye(r);
%% Start semi-smooth newton iterations
k = 0;
a_km1 = zeros(r,1); % a0
L_km1 = L(a_km1);
while true
    k = k + 1;
    J_L_km1 = J_L(a_km1);
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
    L_k = L(a_k);
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
    L_km1 = L_k;
end
%% Now that we have the root of L, we apply use it to obtain x(alphastar)
xstar = xa(a_k);

if runCVX 
    %% Debug with CVX if necessary
    B = 1/h0 * eye(n) + U*U' - V*V';
    % cvx
    tic
    f = @(x)1/2*dot(x-y, B*(x-y)) ;
    cvx_begin quiet
            variable xRef(n)
            minimize f(xRef)
            subject to 
            0 <= xRef
    cvx_end 
    cvxErr = norm(xRef - xstar) / norm(xRef);
    assert(cvxErr < 1e-4, 'Compared to CVX we have a bad answer')
end

end % prox