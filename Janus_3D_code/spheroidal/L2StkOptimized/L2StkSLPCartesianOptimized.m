function [res_x, res_y, res_z] = L2StkSLPCartesianOptimized(X_eval, params, sigma_x, sigma_y, sigma_z, source_Gmatrix, target_precomp)
%{
Cartesian-family evaluation for the Stokes SLP.

This applies the SLP through
    sigma = rho_+ ehat_+ + rho_- ehat_- + rho_0 e_z,
where
    ehat_+ = (e_x + i e_y)/sqrt(2),
    ehat_- = (e_x - i e_y)/sqrt(2),
    rho_+ = (sigma_x - i sigma_y)/sqrt(2),
    rho_- = (sigma_x + i sigma_y)/sqrt(2),
    rho_0 = sigma_z.

For each Cartesian family it evaluates
    S[W] = 0.5 * (S^(I) + S^(y) - S^(x)),

Self evaluation returns the continuous surface trace.
%}

p = params.p;
u0 = params.u0;
a = params.a;
oblate = params.oblate;
if nargin < 7
    target_precomp = [];
end

[np, nf, ns] = size(sigma_x);

if isscalar(u0)
    u0 = u0 .* ones(1, ns);
end
if isscalar(a)
    a = a .* ones(1, ns);
end
if isscalar(oblate)
    oblate = oblate .* ones(1, ns);
end

if isempty(X_eval)
    X_eval = LOCAL_source_geometry_all(p, u0, a, oblate, ns);
    allow_surface = true;
else
    allow_surface = false;
end

[rho_plus, rho_minus, rho_zero] = LOCAL_build_cartesian_family_densities(sigma_x, sigma_y, sigma_z);
[shc_plus, shc_minus, shc_zero] = LOCAL_shAna_batched(rho_plus, rho_minus, rho_zero);
[Gshc_plus, Gshc_minus, Gshc_zero] = LOCAL_apply_source_Gmatrix(source_Gmatrix, shc_plus, shc_minus, shc_zero);

if isempty(target_precomp)
    target_precomp = spheroidal_cartesian_target_precompute(p, u0, a, oblate, X_eval, false, allow_surface);
end

Gshc_main = cat(2, Gshc_plus, Gshc_minus, Gshc_zero);
SL_main = LOCAL_scalar_single_layer(target_precomp, Gshc_main);

Gshc_weighted = LOCAL_apply_source_coordinate_coupling_matrices(p, u0, a, oblate, source_Gmatrix, Gshc_main, nf);

Gshc_derivative = cat(2, Gshc_main, Gshc_weighted);
[dplus_all, dminus_all, dzero_all] = cartesian_family_directional_derivatives(p, u0, a, oblate, Gshc_derivative, X_eval, target_precomp);

idx_main = 1:(3 * nf);
idx_weighted = (3 * nf) + (1:(9 * nf));

dplus_main = slice_density(dplus_all, idx_main);
dminus_main = slice_density(dminus_all, idx_main);
dzero_main = slice_density(dzero_all, idx_main);

dplus_weighted = slice_density(dplus_all, idx_weighted);
dminus_weighted = slice_density(dminus_all, idx_weighted);
dzero_weighted = slice_density(dzero_all, idx_weighted);

[res_x, res_y, res_z] = LOCAL_combine_family_sl_terms( ...
    X_eval, SL_main, dplus_main, dminus_main, dzero_main, ...
    dplus_weighted, dminus_weighted, dzero_weighted, nf);

if params.isReal && isreal(sigma_x) && isreal(sigma_y) && isreal(sigma_z)
    res_x = real(res_x);
    res_y = real(res_y);
    res_z = real(res_z);
end
end

function [rho_plus, rho_minus, rho_zero] = LOCAL_build_cartesian_family_densities(sigma_x, sigma_y, sigma_z)
    rho_plus = (sigma_x - 1i .* sigma_y) ./ sqrt(2);
    rho_minus = (sigma_x + 1i .* sigma_y) ./ sqrt(2);
    rho_zero = sigma_z;
end

function varargout = LOCAL_shAna_batched(varargin)
    varargout = cell(1, nargin);
    for j = 1:nargin
        density = varargin{j};
        [~, nf, ns] = size(density);
        shc = zeros(size(shAna(density(:, :, 1)), 1), nf, ns);
        for k = 1:ns
            shc(:,:,k) = shAna(density(:,:,k));
        end
        varargout{j} = shc;
    end
end

function varargout = LOCAL_apply_source_Gmatrix(source_Gmatrix, varargin)
    varargout = cell(1, nargin - 1);
    for j = 1:nargin - 1
        shc = varargin{j};
        [sp, nf, ns] = size(shc);
        Gshc = zeros(sp, nf, ns);

        for k = 1:ns
            if isempty(source_Gmatrix)
                Gshc(:,:,k) = shc(:,:,k);
            elseif iscell(source_Gmatrix)
                if isempty(source_Gmatrix{k})
                    Gshc(:,:,k) = shc(:,:,k);
                else
                    Gshc(:,:,k) = source_Gmatrix{k} \ shc(:,:,k);
                end
            else
                Gshc(:,:,k) = source_Gmatrix \ shc(:,:,k);
            end
        end
        varargout{j} = Gshc;
    end
end

function Gshc_weighted = LOCAL_apply_source_coordinate_coupling_matrices(p, u0, a, oblate, source_Gmatrix, Gshc_main, nf)
    [sp, ~, ns] = size(Gshc_main);
    Gshc_weighted = zeros(sp, 9 * nf, ns);

    for k = 1:ns
        source_coupling = ...
            LOCAL_source_coordinate_coupling_matrices(p, u0(k), a(k), oblate(k), LOCAL_get_Gmatrix_k(source_Gmatrix, k));

        Gmain_k = Gshc_main(:,:,k);
        Gshc_weighted(:,:,k) = [ ...
            source_coupling.yplus * Gmain_k(:, 1:nf), ...
            source_coupling.yminus * Gmain_k(:, 1:nf), ...
            source_coupling.yzero * Gmain_k(:, 1:nf), ...
            source_coupling.yplus * Gmain_k(:, nf + (1:nf)), ...
            source_coupling.yminus * Gmain_k(:, nf + (1:nf)), ...
            source_coupling.yzero * Gmain_k(:, nf + (1:nf)), ...
            source_coupling.yplus * Gmain_k(:, 2 * nf + (1:nf)), ...
            source_coupling.yminus * Gmain_k(:, 2 * nf + (1:nf)), ...
            source_coupling.yzero * Gmain_k(:, 2 * nf + (1:nf))];
    end
end

function source_coupling = LOCAL_source_coordinate_coupling_matrices(p, u0, a, oblate, source_Gmatrix)
    persistent key_cache source_coupling_cache
    if isempty(key_cache)
        key_cache = {};
        source_coupling_cache = {};
    end

    key = sprintf('%d|%.16g|%.16g|%d', p, u0, a, double(oblate));
    cache_id = find(strcmp(key_cache, key), 1);
    if ~isempty(cache_id)
        source_coupling = source_coupling_cache{cache_id};
        return
    end

    xy_alpha = a * sqrt(u0^2 + 2 * double(logical(oblate)) - 1) * sqrt(4 * pi / 3);
    z_alpha = a * u0 * sqrt(4 * pi / 3);

    source_coupling = struct();
    source_coupling.yplus = LOCAL_apply_Gmatrix_similarity( ...
        source_Gmatrix, LOCAL_source_coordinate_coupling_matrix(p, 1, xy_alpha) ...
    );
    source_coupling.yminus = LOCAL_apply_Gmatrix_similarity( ...
        source_Gmatrix, LOCAL_source_coordinate_coupling_matrix(p, -1, -xy_alpha) ...
    );
    source_coupling.yzero = LOCAL_apply_Gmatrix_similarity( ...
        source_Gmatrix, LOCAL_source_coordinate_coupling_matrix(p, 0, z_alpha) ...
    );

    key_cache{end + 1} = key;
    source_coupling_cache{end + 1} = source_coupling;
end

function transformed_matrix = LOCAL_apply_Gmatrix_similarity(source_Gmatrix, coupling_matrix)
    % Derive and cache the change-of-basis needed for this.
    % The y_k coupling matrix M_y is derived to apply to shc.
    % Thus, we actually need shc_{transformed} = M_y shc_{original}
    % We use Gshc = G^{-1} shc_{original}, so applying y_k directly to Gshc needs
    % the transformation G^{-1} M_y G.
    if isempty(source_Gmatrix)
        transformed_matrix = coupling_matrix;
    else
        transformed_matrix = full(source_Gmatrix \ (coupling_matrix * source_Gmatrix));
    end
end

function source_Gmatrix_k = LOCAL_get_Gmatrix_k(source_Gmatrix, k)
    if isempty(source_Gmatrix)
        source_Gmatrix_k = [];
    elseif iscell(source_Gmatrix)
        source_Gmatrix_k = source_Gmatrix{k};
    else
        source_Gmatrix_k = source_Gmatrix;
    end
end

function source_coupling_matrix = LOCAL_source_coordinate_coupling_matrix(p, k, alpha_k)
    % Build the sparse mapping for y_k * Y_n^m. Multiplication by the
    % degree-1 position vector y only couples one input mode (n,m) to
    % output modes (n-1,m-k) and (n+1,m-k), with coefficients given by the
    % triple-product/Wigner-3j formula.
    sp = (p + 1)^2;
    max_nnz = 2 * sp;

    row_idx = zeros(max_nnz, 1);
    col_idx = zeros(max_nnz, 1);
    vals = complex(zeros(max_nnz, 1));
    nnz_count = 0;

    for n = 0:p
        for m = -n:n
            in_idx = n^2 + n + m + 1;
            ell_m = m - k;

            for ell = [n - 1, n + 1]
                if ell < 0 || ell > p || abs(ell_m) > ell
                    continue
                end

                nnz_count = nnz_count + 1;
                row_idx(nnz_count) = ell^2 + ell + ell_m + 1;
                col_idx(nnz_count) = in_idx;
                vals(nnz_count) = alpha_k * ((-1) .^ (m - k)) .* ...
                                    sqrt(3 * (2 * n + 1) * (2 * ell + 1) / (4 * pi)) .* ...
                                    wigner3j(1, n, ell, 0, 0, 0) .* ...
                                    wigner3j(1, n, ell, -k, m, k - m);
            end
        end
    end

    source_coupling_matrix = sparse(row_idx(1:nnz_count), col_idx(1:nnz_count), vals(1:nnz_count), sp, sp);
end

function SL = LOCAL_scalar_single_layer(target_precomp, Gshc)
    [~, nf, ns] = size(Gshc);
    nt = target_precomp.sources(1).nt;
    SL = zeros(nt, nf, ns);

    for k = 1:ns
        source_precomp = target_precomp.sources(k);
        if source_precomp.nt == 0
            continue;
        end

        SLk = zeros(source_precomp.nt, nf);
        for r = 1:numel(source_precomp.regions)
            region = source_precomp.regions(r);
            idx = region.idx;
            SLk(idx, :) = ...
                ((region.Fr .* region.Yr0) .* (region.a .* region.source_radial_prefactors)) * Gshc(:,:,k);
        end

        SL(:,:,k) = SLk;
    end
end

function [res_x, res_y, res_z] = LOCAL_combine_family_sl_terms( ...
        X_eval, SL_main, dplus_main, dminus_main, dzero_main, ...
        dplus_weighted, dminus_weighted, dzero_weighted, nf)
    % Assemble final result (i.e. the SLP Cartesian formula)
    idx_plus = 1:nf;
    idx_minus = nf + (1:nf);
    idx_zero = 2*nf + (1:nf);

    idx_yplus_plus = 1:nf;
    idx_yminus_plus = nf + (1:nf);
    idx_yzero_plus = 2*nf + (1:nf);
    idx_yplus_minus = 3*nf + (1:nf);
    idx_yminus_minus = 4*nf + (1:nf);
    idx_yzero_minus = 5*nf + (1:nf);
    idx_yplus_zero = 6*nf + (1:nf);
    idx_yminus_zero = 7*nf + (1:nf);
    idx_yzero_zero = 8*nf + (1:nf);

    ehat_plus = [1, 1i, 0] ./ sqrt(2);
    ehat_minus = [1, -1i, 0] ./ sqrt(2);

    SL_plus = slice_density(SL_main, idx_plus);
    SL_minus = slice_density(SL_main, idx_minus);
    SL_zero = slice_density(SL_main, idx_zero);

    SI_x = (SL_plus + SL_minus) ./ sqrt(2);
    SI_y = 1i .* (SL_plus - SL_minus) ./ sqrt(2);
    SI_z = SL_zero;

    scalar_x_term = ...
        slice_density(dplus_main, idx_plus) ./ sqrt(2) + ...
        slice_density(dminus_main, idx_minus) ./ sqrt(2) + ...
        slice_density(dzero_main, idx_zero);

    x_eval = X_eval(:, 1, :);
    y_eval = X_eval(:, 2, :);
    z_eval = X_eval(:, 3, :);

    SX_x = x_eval .* scalar_x_term;
    SX_y = y_eval .* scalar_x_term;
    SX_z = z_eval .* scalar_x_term;

    Sy_plus = ...
        slice_density(dplus_weighted, idx_yplus_plus) ./ sqrt(2) + ...
        slice_density(dminus_weighted, idx_yplus_minus) ./ sqrt(2) + ...
        slice_density(dzero_weighted, idx_yplus_zero);
    Sy_minus = ...
        slice_density(dplus_weighted, idx_yminus_plus) ./ sqrt(2) + ...
        slice_density(dminus_weighted, idx_yminus_minus) ./ sqrt(2) + ...
        slice_density(dzero_weighted, idx_yminus_zero);
    Sy_zero = ...
        slice_density(dplus_weighted, idx_yzero_plus) ./ sqrt(2) + ...
        slice_density(dminus_weighted, idx_yzero_minus) ./ sqrt(2) + ...
        slice_density(dzero_weighted, idx_yzero_zero);

    SY_x = ehat_plus(1) .* Sy_plus + ehat_minus(1) .* Sy_minus;
    SY_y = ehat_plus(2) .* Sy_plus + ehat_minus(2) .* Sy_minus;
    SY_z = Sy_zero;

    res_x = 0.5 .* (SI_x + SY_x - SX_x);
    res_y = 0.5 .* (SI_y + SY_y - SX_y);
    res_z = 0.5 .* (SI_z + SY_z - SX_z);
end

function X_src = LOCAL_source_geometry_all(p, u0, a, oblate, ns)
    np = 2*p*(p+1);
    X_src = zeros(np, 3, ns);

    for k = 1:ns
        if oblate(k)
            X_src(:,:,k) = oblate_spheroid_shape(p, u0(k), a(k));
        else
            X_src(:,:,k) = prolate_spheroid_shape(p, u0(k), a(k));
        end
    end
end
