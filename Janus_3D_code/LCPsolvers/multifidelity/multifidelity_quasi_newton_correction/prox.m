function xstar = prox(y, opts)
    if isempty(opts.U) && isempty(opts.V)
        % Project with respect to the identity
        xstar = max(0, y);
    elseif size(opts.U,2) + size(opts.V,2) == 1
        % The sign on sigma is counter intuitive, but remember B = B0 + UU' - VV'
        if ~isempty(opts.U)
            sigma = -1;
            w = opts.U; 
        else 
            sigma = 1;
            w = opts.V;
        end
        xstar = prox_rank1(y, opts.h0, w, sigma, opts);
        return  
    else
        xstar = prox_rankr(y, opts.h0, opts.U, opts.V, opts);
    end
end % prox



function xstar = prox_rank1(y, h0, w, sigma, opts)

    if all(y >= 0)
        xstar = y;
        return 
    elseif isempty(w)
        xstar = max(y, 0);
        return 
    end

    if ~exist('opts','var') || isempty(opts) || ~isfield(opts, 'prox')
        opts = struct('prox', struct('debug', false));
    end
    %% Input checks
    assert(size(y,2) == 1, 'y should be a column vector');
    assert(size(h0,1) == 1 && size(h0,2) ==1, 'h0 should be a scalar');
    assert(size(w,2) == 1, 'w should be a column vector');
    assert(sigma == 1 || sigma == -1, 'Sigma should be the sign of the rank 1 update')
    %% Find points of discontinuituy in the piecewise linear function
    alphas = - sigma * y ./ (h0 .* w);
    alphas = sort(alphas); % perhaps one day do something smarter than sorting...
    %% Specify the values of the piecewise linear function $L$ at each value of
    % alpha
    la = @(a) y + sigma * a * h0 .* w;
    N = length(y);
    Ls = zeros(N,1);
    for n = 1:N
        a = alphas(n);
        l = la(a);
        mask = l > 0;
        Ls(n) = a + dot(w, y) - dot(w(mask), l(mask));
    end
    %% Specify the interval where the root occurs
    [ix, ~] = find(Ls < 0, 1, 'last');
    if isempty(ix) % This means that all the constraints are active.
        a = alphas(1) - 10; % arbitrary
    else
        a = alphas(ix);
    end
    l = la(a);
    mask = l > 10*eps() | (abs(l) < 10*eps() & sigma * h0 * w > 0);
    alphastar = (dot(w,y) - dot(w(mask), y(mask))) ...
        / (sigma * dot(w(mask), h0 * w(mask)) - 1 );

    xstar = max(y + sigma * alphastar*h0.*w, 0);

    if opts.prox.runCVX 
        %% idiot check
        % H = h0 + sigma * u * u'
        % B = 1 / h0 * I - sigma * w * w';
        n = length(y);
        B = eye(n) ./ h0 - sigma * (w * w');
        try %#ok<TRYNC>
            cvx_begin quiet
                variable z(N)
                minimize( dot(y-z,B*(y-z)) )
                subject to 
                    z >= 0
            cvx_end
            cvxErr = norm(z - xstar) / norm(xstar);
            assert(cvxErr < 1e-4, 'Our solution does not match CVX')
        end 
    end 

end % prox_rank1


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