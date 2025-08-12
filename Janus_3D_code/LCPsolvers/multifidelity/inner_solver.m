function [x_inner_1, info_inner] = inner_solver(x_outer_1, grad_outer_1, grad_inner, s_1, y_1, outer_iter, opts)
    %This is a wrapper function that takes one initial high fidelity step (using whatever method is prescribed by opts.step) and then takes a set amount of steps or til a tolerance of the corrected low fidelity operator.
    switch lower(opts.outer.solver)
        case 'bb'
            x_inner_0 = x_outer_1 - ((s_1'*s_1)/(s_1'*y_1))*grad_outer_1;
        case 'prox'
            opts.outer.prox = updateH(outer_iter, s_1, y_1, opts.outer.prox);
            p = -opts.outer.prox.H(grad_outer_1);
            x_inner_0 = prox(x_outer_1 + p, opts.outer.prox);
    end

    %check if we want to take any inner steps or just high fidelity
    if opts.inner.enabled == false
        x_inner_1 = x_inner_0;
        info_inner = [];
        return
    else
        switch lower(opts.inner.solver)
            case 'outer prox qn'
                opts.inner.prox = opts.outer.prox;
                [x_inner_1, info_inner] = outer_prox_qn(x_inner_0, grad_inner, opts);

        end

    end

end

function opts = updateH(outer_iter, s_1, y_1, opts)

    %{
    if outer_iter == 0
        % Initial step we do not compute any curvature information because only
        % information about x0 is known, thus no secant equaion could be
        % satisfied
        h0 = 1; 
        H = @(x) x; 
        U = []; 
        V = [];
        return 
    end
    %}

    assert(0 < opts.tau_min, "tau_min must be positive");
    assert(opts.tau_min < opts.tau_max, "tau_max must be larger than tau_min");
    assert(0 < opts.gamma && opts.gamma < 1, "gamma must be in (0,1)");
    n = length(s_1);
    tau_bb2 = dot(s_1,y_1) / norm(y_1,2)^2;
    tau_bb2 = clip(tau_bb2, opts.tau_min, opts.tau_max);
    if tau_bb2 == opts.tau_min
        warning('Convexity of cost function is stagnating'); 
    end
    opts.h0 = opts.gamma * tau_bb2;

    r = min(outer_iter - 1, opts.r);
    % Handeling the case where we have skipped previous hessian updates
    skipped_r = find(arrayfun(@(i)all(opts.S(:,i) == 0),1:r),1,'first');
    if ~isempty(skipped_r) && r < outer_iter
        r = max(skipped_r-1, 1);
    end

    % Let memory fall out of context window regardless of wether we skip or not
    if outer_iter - 1 > opts.r
        opts.S(:,1:r-1) = opts.S(:,2:r);
        opts.Y(:,1:r-1) = opts.Y(:,2:r);
        opts.S(:,r) = 0;
        opts.Y(:,r) = 0;
    end
    % Curvature check
    rho = 1 / dot(y_1,s_1);
    if 1/rho >= 1e-8
        opts.S(:,r) = s_1;
        opts.Y(:,r) = y_1;
    else
        r = r - 1;
    end
    S = opts.S(:,1:r);
    Y = opts.Y(:,1:r);

    switch lower(opts.qnUpdate)
        case 'bfgs'
            U = zeros(n, r);
            V = zeros(n, r);
            H = eye(n) * opts.h0;
            B = eye(n) / opts.h0;
            for j = 1:r
                s = S(:,j); 
                y = Y(:,j); 
                rho = 1 / dot(y,s);
                if 1/rho < 1e-8
                    % TODO something better than skip if the curvature condition is bad
                    continue
                end
                UU = (eye(n) - rho*y*s');
                H = UU'*H*UU + rho*(s*s');
            
                v = B*s / sqrt(dot(s, B*s));
                u = y / sqrt(dot(y,s));
            
                B = B + u*u' - v*v' ;
                U(:, j) = u;
                V(:, j) = v;
            end
            % IDIOT CHECK 
            %assert(norm(B*H - eye(n)) < 1e-8)
            % TODO: Make this matrix free...
            opts.U = U;
            opts.V = V;
            opts.H = @(x) H*x;
            opts.B = B;
        %{
        case 'sr1'
            U = zeros(n, r);
            V = zeros(n, r);
            H = eye(n) * h0;
            B = eye(n) / h0;
            for j = 1:r 
                s = S(:,j); 
                y = Y(:,j);
                % update H (and mem)
                denomH = dot(s - H*y,y);
                sigma = sign(denomH);
                denomH = sqrt(sigma*denomH);
                assert(isreal(denomH));
                u = (s - H*y) / denomH;
                H = H + sigma*(u*u');
                % update H (and mem)
                denomB = dot(y - B*s,s);
                lambda = sign(denomB);
                denomB = sqrt(lambda*denomB);
                assert(isreal(denomB));
                v = (y - B*s) / denomB;
                B = B + lambda*(v*v');
                % Store for proof of concept
                if lambda > 0 
                    U(:,j) = v;
                else
                    V(:,j) = v;
                end
            end
            mask = arrayfun(@(i) all(U(:,i) == 0),1:r);
            U = U(:,~mask);
            V = V(:,mask);
            assert(norm(B*H - eye(n))/norm(H)/norm(B) < 1e-8)
            % TODO: Make this matrix free...
            H = @(x) H*x;
            % We may have an indefinite matrix from the SR1 update
            if ~isempty(U)  
                % TODO do this better with a shifted power method
                assert(min(eig(H)) > 0)
            end
        otherwise
            error([opts.qnUpdate ' update not implement'])
        %}
    end
end 

