function [x, flag, relres, iter, info] = SpheroidalMS_augmented_GMRES(A, b, opts)
%{
Augmented GMRES wrapper with a fixed supplied basis.

The basis is meant to target the current hard modes of the operator
(for example, active-contact body modes). The preconditioner is absorbed
into the operator so the coarse correction is built for the same
left-preconditioned system that GMRES actually sees.
%}

if nargin < 3 || isempty(opts)
    opts = struct();
end
restart = opts.rst;
tol = opts.tol;
maxit = opts.maxit;
prec = opts.prec;
U = LOCAL_modified_gram_schmidt(opts.deflate_basis);
base_solver = LOCAL_normalize_krylov_solver(opts.krylov_solver);

A_op = LOCAL_get_matvec_func(A);
Abar = @(x) LOCAL_apply_preconditioner(prec, A_op(x));
bbar = LOCAL_apply_preconditioner(prec, b);

info = struct( ...
    'basis_dim', 0, ...
    'coarse_relres', 1.0, ...
    'base_solver', base_solver, ...
    'krylov_cycles', 0, ...
    'krylov_total_iters', 0, ...
    'krylov_recycle_dim', 0 ...
);

% No augmented basis, so run usual solver on preconditioned system
if isempty(U)
    fprintf('No augmentation: running the usual GMRES(m).');
    [x, flag, relres, iter, current_info] = LOCAL_run_krylov_solver(base_solver, Abar, bbar, restart, tol, maxit, []);
    info = LOCAL_update_krylov_info(info, current_info);
    return;
end

AU = LOCAL_apply_operator(Abar, U);
[U, AU] = LOCAL_prune_basis(U, AU);

% Supplied basis does not give a meaningful augmented base, so fallback to the usual solver
if isempty(U)
    [x, flag, relres, iter, current_info] = LOCAL_run_krylov_solver(base_solver, Abar, bbar, restart, tol, maxit, []);
    info = LOCAL_update_krylov_info(info, current_info);
    return;
end

H = U' * AU; % U^TAU
rhs_small = U' * bbar; % U^Tb
alpha = LOCAL_augmented_solve(H, rhs_small);
x0 = U * alpha; % Lift from coarse space to full space
r_coarse = bbar - AU * alpha;

info.basis_dim = size(U, 2);
info.coarse_relres = LOCAL_relative_norm(r_coarse, bbar);

% Converged, so return
if info.coarse_relres <= tol
    x = x0;
    flag = 0;
    relres = info.coarse_relres;
    iter = [0 0];
    info.krylov_recycle_dim = size(U, 2);
    return;
end

[x, flag, relres, iter, current_info] = LOCAL_run_krylov_solver(base_solver, Abar, bbar, restart, tol, maxit, x0);
info = LOCAL_update_krylov_info(info, current_info);
end %% END MAIN FUNCTION

function solver = LOCAL_normalize_krylov_solver(base_solver)
    solver = lower(char(base_solver));
    solver_key = regexprep(solver, '[^a-z0-9]', '');
    switch solver_key
        case 'gmres'
            solver = 'gmres';
        otherwise
            error('Unsupported augmented Krylov solver "%s". Use ''gmres''.', solver);
    end
end

function A_op = LOCAL_get_matvec_func(A)
    if isnumeric(A)
        A_op = @(x) A * x;
    else
        A_op = A;
    end
end

function Y = LOCAL_apply_operator(A_op, X)
    if isempty(X)
        Y = X;
        return;
    end

    if size(X, 2) == 1
        Y = A_op(X);
        return;
    end

    Y = zeros(size(X));
    for col_idx = 1:size(X, 2)
        Y(:, col_idx) = A_op(X(:, col_idx));
    end
end

function Y = LOCAL_apply_preconditioner(prec, X)
    if isempty(prec)
        Y = X;
        return;
    end

    if isnumeric(prec)
        Y = prec \ X;
    else
        Y = prec(X);
    end
end

function [U, AU] = LOCAL_prune_basis(U, AU)
    if isempty(U)
        return;
    end

    % Drop basis vectors whose image under the preconditioned operator is
    % numerically zero; they cannot improve the coarse correction.
    norms_AU = vecnorm(AU);
    keep = norms_AU > 1e-12 * max(1, max(norms_AU));
    U = U(:, keep);
    AU = AU(:, keep);
    if isempty(U)
        return;
    end

    % Keep only a numerically independent subset of Abar(U). The matching
    % columns of U define the reduced coarse space used by the augmentation.
    [~, R, perm] = qr(AU, 'vector');
    diagR = abs(diag(R));
    if isempty(diagR)
        U = [];
        AU = [];
        return;
    end

    rank_AU = nnz(diagR > 1e-10 * diagR(1));
    rank_AU = min(rank_AU, numel(perm));
    if rank_AU == 0
        U = [];
        AU = [];
        return;
    end

    keep_cols = perm(1:rank_AU);
    U = U(:, keep_cols);
    AU = AU(:, keep_cols);
end

function x = LOCAL_augmented_solve(H, rhs)
    if isempty(H)
        x = zeros(0, size(rhs, 2));
        return;
    end

    if size(H, 1) == 1
        if abs(H) < 1e-14
            x = zeros(size(rhs));
        else
            x = rhs / H;
        end
        return;
    end

    if rcond(H) > 1e-12
        x = H \ rhs;
    else
        x = pinv(H) * rhs;
    end
end

function rel = LOCAL_relative_norm(r, b)
    denom = norm(b);
    if denom == 0
        rel = norm(r);
    else
        rel = norm(r) / denom;
    end
end

function [x, flag, relres, iter, info] = LOCAL_run_krylov_solver(base_solver, Abar, bbar, restart, tol, maxit, x0)
    if nargin < 7
        x0 = [];
    end

    info = struct('cycles', 0, 'total_iters', 0, 'recycle_dim', 0);
    switch base_solver
        case 'gmres'
            if isempty(x0)
                [x, flag, relres, iter] = gmres(Abar, bbar, restart, tol, maxit);
            else
                [x, flag, relres, iter] = gmres(Abar, bbar, restart, tol, maxit, [], [], x0);
            end
            info.cycles = iter(1);
            info.total_iters = LOCAL_total_gmres_iters(iter, restart);
        otherwise
            error('Unsupported augmented Krylov solver "%s". Use ''gmres''.', base_solver);
    end
end

function info = LOCAL_update_krylov_info(info, current_info)
    info.krylov_cycles = current_info.cycles;
    info.krylov_total_iters = current_info.total_iters;
    info.krylov_recycle_dim = current_info.recycle_dim;
end

function total_iters = LOCAL_total_gmres_iters(it, restart)
    if numel(it) ~= 2
        total_iters = NaN;
        return;
    end
    if it(1) <= 0
        total_iters = it(2);
    else
        total_iters = (it(1) - 1) * restart + it(2);
    end
end

function Q = LOCAL_modified_gram_schmidt(X)
    % Process columns of X
    if isempty(X)
        Q = [];
        return;
    end
    keep = all(isfinite(X), 1) & (vecnorm(X) > 1e-12);
    X = X(:, keep);

    n = size(X, 1);
    Q = zeros(n, 0);
    for col_idx = 1:size(X, 2)
        v = X(:, col_idx);
        for q_idx = 1:size(Q, 2)
            v = v - Q(:, q_idx) * (Q(:, q_idx)' * v);
        end

        nv = norm(v);
        if nv > 1e-10
            Q(:, end+1) = v / nv;
        end
    end
end
