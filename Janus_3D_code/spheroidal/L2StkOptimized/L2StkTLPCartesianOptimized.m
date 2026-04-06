function [res_x, res_y, res_z] = L2StkTLPCartesianOptimized(X_eval, nu_eval, params, sigma_x, sigma_y, sigma_z, source_Gmatrix, dealiasing_flag, dealiasing_pad, target_precomp)
%{
Cartesian-family evaluation for the TSL.

This applies the TSL through the Cartesian family decomposition
    sigma = rho_+ ehat_+ + rho_- ehat_- + rho_0 e_z,
where
    ehat_+ = (e_x + i e_y)/sqrt(2),
    ehat_- = (e_x - i e_y)/sqrt(2),
    rho_+ = (sigma_x - i sigma_y)/sqrt(2),
    rho_- = (sigma_x + i sigma_y)/sqrt(2),
    rho_0 = sigma_z.

It returns the principal-value for self-evaluation, matching the other spectral layer potential codes.
%}

if nargin < 8 || isempty(dealiasing_flag)
    dealiasing_flag = false;
end
if nargin < 9 || isempty(dealiasing_pad)
    dealiasing_pad = 4;
end
if nargin < 10
    target_precomp = [];
end

[~, nf, ns] = size(sigma_x);

p = params.p;
u0 = params.u0;
a = params.a;
oblate = params.oblate;

if isscalar(u0)
    u0 = u0 .* ones(1, ns);
end
if isscalar(a)
    a = a .* ones(1, ns);
end
if isscalar(oblate)
    oblate = oblate .* ones(1, ns);
end

if dealiasing_flag
    fine_p = p + dealiasing_pad;
else
    fine_p = p;
end

[~, rho_plus, rho_minus, rho_zero, rho_y_plus, rho_y_minus, rho_y_zero] = ...
    LOCAL_build_cartesian_family_densities(p, fine_p, u0, a, oblate, sigma_x, sigma_y, sigma_z, dealiasing_flag);

if isempty(X_eval)
    [res_x, res_y, res_z] = ...
        LOCAL_TSL_surface_evaluation(p, fine_p, u0, a, oblate, sigma_x, sigma_y, sigma_z, source_Gmatrix, dealiasing_flag);
    if params.isReal
        res_x = real(res_x);
        res_y = real(res_y);
        res_z = real(res_z);
    end
    return;
end

[shc_plus, shc_minus, shc_zero, shc_y_plus, shc_y_minus, shc_y_zero] = ...
    LOCAL_shAna_batched(rho_plus, rho_minus, rho_zero, rho_y_plus, rho_y_minus, rho_y_zero);
[Gshc_plus, Gshc_minus, Gshc_zero, Gshc_y_plus, Gshc_y_minus, Gshc_y_zero] = ...
    LOCAL_apply_source_Gmatrix(source_Gmatrix, shc_plus, shc_minus, shc_zero, shc_y_plus, shc_y_minus, shc_y_zero);

Gshc_main = cat(2, Gshc_plus, Gshc_minus, Gshc_zero);
Gshc_weighted = cat(2, Gshc_y_plus, Gshc_y_minus, Gshc_y_zero);
if isempty(target_precomp)
    target_precomp = spheroidal_cartesian_target_precompute(p, u0, a, oblate, X_eval, true);
end

[dplus_all, dminus_all, dzero_all] = ...
    cartesian_family_directional_derivatives(p, u0, a, oblate, Gshc_main, X_eval, target_precomp);

Gshc_hessian = cat(2, Gshc_main, Gshc_weighted);
[Hn_x_all, Hn_y_all, Hn_z_all] = LOCAL_hessian_normal( ...
    nu_eval, target_precomp, Gshc_hessian);

idx_plus = 1:nf;
idx_minus = nf + (1:nf);
idx_zero = 2*nf + (1:nf);
idx_y_plus = 3*nf + idx_plus;
idx_y_minus = 3*nf + idx_minus;
idx_y_zero = 3*nf + idx_zero;

[res_x, res_y, res_z] = LOCAL_combine_family_traction_terms( ...
    X_eval, nu_eval, dplus_all, dminus_all, dzero_all, ...
    Hn_x_all, Hn_y_all, Hn_z_all, ...
    idx_plus, idx_minus, idx_zero, idx_y_plus, idx_y_minus, idx_y_zero);

if params.isReal && isreal(sigma_x) && isreal(sigma_y) && isreal(sigma_z)
    res_x = real(res_x);
    res_y = real(res_y);
    res_z = real(res_z);
end
end %% END MAIN FUNCTION

function [res_x, res_y, res_z] = LOCAL_TSL_surface_evaluation( ...
        p, fine_p, u0, a, oblate, sigma_x, sigma_y, sigma_z, source_Gmatrix, dealiasing_flag)
    [np, nf, ns] = size(sigma_x);
    res_x = zeros(np, nf, ns);
    res_y = zeros(np, nf, ns);
    res_z = zeros(np, nf, ns);

    for k = 1:ns
        Tself = LOCAL_get_averaged_backend_self_operator( ...
            p, fine_p, u0(k), a(k), oblate(k), LOCAL_get_source_Gmatrix_k(source_Gmatrix, k), dealiasing_flag);
        out = Tself * [sigma_x(:,:,k); sigma_y(:,:,k); sigma_z(:,:,k)];
        [res_x(:,:,k), res_y(:,:,k), res_z(:,:,k)] = LOCAL_unstack_cartesian_density(out, np);
    end
end

function source_Gmatrix_k = LOCAL_get_source_Gmatrix_k(source_Gmatrix, k)
    if iscell(source_Gmatrix)
        source_Gmatrix_k = source_Gmatrix{k};
    else
        source_Gmatrix_k = source_Gmatrix;
    end
end

function Tself = LOCAL_get_averaged_backend_self_operator(p, fine_p, u0, a, oblate, source_Gmatrix, dealiasing_flag)
    % Active self implementation for both prolates and oblates: build the
    % one-sided surface traces with the backend itself, then cache their
    % principal-value average. The exact formula-cache route is still
    % retained below as an inactive alternative.
    persistent keys values
    if isempty(keys)
        keys = strings(0, 1);
        values = {};
    end

    key = sprintf('backend-self|p=%d|u0=%.17g|a=%.17g|oblate=%d', ...
        p, u0, a, double(oblate));
    idx = find(keys == key, 1);
    if ~isempty(idx)
        Tself = values{idx};
        return;
    end

    Text = LOCAL_build_dense_surface_interaction_matrix(p, fine_p, u0, a, oblate, source_Gmatrix, dealiasing_flag, false);
    Tint = LOCAL_build_dense_surface_interaction_matrix(p, fine_p, u0, a, oblate, source_Gmatrix, dealiasing_flag, true);
    Tself = 0.5 .* (Text + Tint);

    keys(end + 1, 1) = string(key);
    values{end + 1, 1} = Tself;
end

function Ttrace = LOCAL_build_dense_surface_interaction_matrix(p, fine_p, u0, a, oblate, source_Gmatrix, dealiasing_flag, is_interior)
    np = 2*p*(p+1);
    X_src = LOCAL_source_geometry(p, u0, a, oblate);
    nu_eval = get_norm_vecs(p, u0, oblate);
    I = eye(np);
    Z = zeros(np, np);

    [tx, ty, tz] = LOCAL_get_surface_trace( ...
        p, fine_p, u0, a, oblate, I, Z, Z, source_Gmatrix, X_src, nu_eval, dealiasing_flag, is_interior);
    T_x = [tx; ty; tz];

    [tx, ty, tz] = LOCAL_get_surface_trace( ...
        p, fine_p, u0, a, oblate, Z, I, Z, source_Gmatrix, X_src, nu_eval, dealiasing_flag, is_interior);
    T_y = [tx; ty; tz];

    [tx, ty, tz] = LOCAL_get_surface_trace( ...
        p, fine_p, u0, a, oblate, Z, Z, I, source_Gmatrix, X_src, nu_eval, dealiasing_flag, is_interior);
    T_z = [tx; ty; tz];

    Ttrace = [T_x, T_y, T_z];
end

function [res_x, res_y, res_z] = LOCAL_get_surface_trace( ...
        p, fine_p, u0, a, oblate, sigma_x, sigma_y, sigma_z, source_Gmatrix, X_src, nu_eval, dealiasing_flag, is_interior)
    [np, nf] = size(sigma_x);
    sigma_x_3d = reshape(sigma_x, np, nf, 1);
    sigma_y_3d = reshape(sigma_y, np, nf, 1);
    sigma_z_3d = reshape(sigma_z, np, nf, 1);
    X_src_3d = reshape(X_src, np, 3, 1);
    nu_eval_3d = reshape(nu_eval, np, 3, 1);

    [~, rho_plus, rho_minus, rho_zero, rho_y_plus, rho_y_minus, rho_y_zero] = ...
        LOCAL_build_cartesian_family_densities(p, fine_p, u0, a, oblate, sigma_x_3d, sigma_y_3d, sigma_z_3d, dealiasing_flag);

    [shc_plus, shc_minus, shc_zero, shc_y_plus, shc_y_minus, shc_y_zero] = ...
        LOCAL_shAna_batched(rho_plus, rho_minus, rho_zero, rho_y_plus, rho_y_minus, rho_y_zero);
    [Gshc_plus, Gshc_minus, Gshc_zero, Gshc_y_plus, Gshc_y_minus, Gshc_y_zero] = ...
        LOCAL_apply_source_Gmatrix(source_Gmatrix, shc_plus, shc_minus, shc_zero, shc_y_plus, shc_y_minus, shc_y_zero);

    Gshc_main = cat(2, Gshc_plus, Gshc_minus, Gshc_zero);
    Gshc_weighted = cat(2, Gshc_y_plus, Gshc_y_minus, Gshc_y_zero);
    target_precomp = spheroidal_cartesian_target_precompute(p, u0, a, oblate, X_src_3d, true, true, is_interior);

    [dplus_all, dminus_all, dzero_all] = cartesian_family_directional_derivatives( ...
        p, u0, a, oblate, Gshc_main, X_src_3d, target_precomp);

    Gshc_hessian = cat(2, Gshc_main, Gshc_weighted);
    [Hn_x_all, Hn_y_all, Hn_z_all] = LOCAL_hessian_normal(nu_eval_3d, target_precomp, Gshc_hessian);

    idx_plus = 1:nf;
    idx_minus = nf + (1:nf);
    idx_zero = 2*nf + (1:nf);
    idx_y_plus = 3*nf + idx_plus;
    idx_y_minus = 3*nf + idx_minus;
    idx_y_zero = 3*nf + idx_zero;

    [res_x, res_y, res_z] = LOCAL_combine_family_traction_terms( ...
        X_src_3d, nu_eval_3d, dplus_all, dminus_all, dzero_all, ...
        Hn_x_all, Hn_y_all, Hn_z_all, ...
        idx_plus, idx_minus, idx_zero, idx_y_plus, idx_y_minus, idx_y_zero);

    res_x = res_x(:, :, 1);
    res_y = res_y(:, :, 1);
    res_z = res_z(:, :, 1);
end

function [sigma_x, sigma_y, sigma_z] = LOCAL_unstack_cartesian_density(rhs, np)
    sigma_x = rhs(1:np,:);
    sigma_y = rhs(np + 1:2*np,:);
    sigma_z = rhs(2*np + 1:3*np,:);
end

function [X_src, rho_plus, rho_minus, rho_zero, rho_y_plus, rho_y_minus, rho_y_zero] = ...
        LOCAL_build_cartesian_family_densities(p, fine_p, u0, a, oblate, sigma_x, sigma_y, sigma_z, dealiasing_flag)
    [np, nf, ns] = size(sigma_x);

    X_src = zeros(np, 3, ns);
    rho_plus = (sigma_x - 1i .* sigma_y) ./ sqrt(2);
    rho_minus = (sigma_x + 1i .* sigma_y) ./ sqrt(2);
    rho_zero = sigma_z;

    rho_y_plus = zeros(np, nf, ns);
    rho_y_minus = zeros(np, nf, ns);
    rho_y_zero = zeros(np, nf, ns);

    if ~dealiasing_flag
        [theta, phi] = gl_grid(p);
        theta = reshape(theta, np, 1, 1);
        phi = reshape(phi, np, 1, 1);

        xy_scale = reshape(a(:) .* sqrt(u0(:).^2 + 2 .* double(oblate(:)) - 1), 1, 1, ns);
        z_scale = reshape(a(:) .* u0(:), 1, 1, ns);

        x_src = xy_scale .* sin(theta) .* cos(phi);
        y_src = xy_scale .* sin(theta) .* sin(phi);
        z_src = z_scale .* cos(theta);

        X_src(:, 1, :) = x_src;
        X_src(:, 2, :) = y_src;
        X_src(:, 3, :) = z_src;

        rho_y_plus = (x_src + 1i .* y_src) .* rho_plus ./ sqrt(2);
        rho_y_minus = (x_src - 1i .* y_src) .* rho_minus./ sqrt(2);
        rho_y_zero = z_src .* rho_zero;
        return;
    end

    for i = 1:ns
        X_src_i = LOCAL_source_geometry(p, u0(i), a(i), oblate(i));
        X_src(:, :, i) = X_src_i;

        X_fine_i = LOCAL_source_geometry(fine_p, u0(i), a(i), oblate(i));

        sigma_x_fine = LOCAL_transfer_to_finer_grid(fine_p, p, sigma_x(:, :, i));
        sigma_y_fine = LOCAL_transfer_to_finer_grid(fine_p, p, sigma_y(:, :, i));
        sigma_z_fine = LOCAL_transfer_to_finer_grid(fine_p, p, sigma_z(:, :, i));

        rho_plus_fine = (sigma_x_fine - 1i .* sigma_y_fine) ./ sqrt(2);
        rho_minus_fine = (sigma_x_fine + 1i .* sigma_y_fine) ./ sqrt(2);
        rho_zero_fine = sigma_z_fine;

        rho_y_plus(:, :, i) = LOCAL_transfer_to_coarser_grid( ...
            p, fine_p, (X_fine_i(:, 1) + 1i .* X_fine_i(:, 2)) .* rho_plus_fine ./ sqrt(2));
        rho_y_minus(:, :, i) = LOCAL_transfer_to_coarser_grid( ...
            p, fine_p, (X_fine_i(:, 1) - 1i .* X_fine_i(:, 2)) .* rho_minus_fine ./ sqrt(2));
        rho_y_zero(:, :, i) = LOCAL_transfer_to_coarser_grid( ...
            p, fine_p, X_fine_i(:, 3) .* rho_zero_fine);
    end
end

function X_src_i = LOCAL_source_geometry(p, u0, a, oblate)
    if oblate
        X_src_i = oblate_spheroid_shape(p, u0, a);
    else
        X_src_i = prolate_spheroid_shape(p, u0, a);
    end
end

function varargout = LOCAL_shAna_batched(varargin)
    varargout = cell(1, nargin);
    for j = 1:nargin
        density = varargin{j};
        [~, nf, ns] = size(density);
        shc = zeros(size(shAna(density(:,:,1)), 1), nf, ns);
        for i = 1:ns
            shc(:,:,i) = shAna(density(:,:,i));
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

        for i = 1:ns
            if isempty(source_Gmatrix)
                Gshc(:,:,i) = shc(:,:,i);
            elseif iscell(source_Gmatrix)
                if isempty(source_Gmatrix{i})
                    Gshc(:,:,i) = shc(:,:,i);
                else
                    Gshc(:,:,i) = source_Gmatrix{i} \ shc(:,:,i);
                end
            else
                Gshc(:,:,i) = source_Gmatrix \ shc(:,:,i);
            end
        end

        varargout{j} = Gshc;
    end
end

function [Hn_x, Hn_y, Hn_z] = LOCAL_hessian_normal(nu_eval, target_precomp, solid_harmonic_Gshc)
    [~, nf_all, ns] = size(solid_harmonic_Gshc);

    Hn_x = zeros(size(nu_eval, 1), nf_all, ns);
    Hn_y = zeros(size(nu_eval, 1), nf_all, ns);
    Hn_z = zeros(size(nu_eval, 1), nf_all, ns);

    for k = 1:ns
        source_precomp = target_precomp.sources(k);
        ntk = source_precomp.nt;
        if ntk == 0, continue; end

        Hxk = zeros(ntk, nf_all);
        Hyk = zeros(ntk, nf_all);
        Hzk = zeros(ntk, nf_all);

        for r = 1:numel(source_precomp.regions)
            region = source_precomp.regions(r);
            idx = region.idx;
            nu_r = nu_eval(idx, :, k);

            [HnU, HnV, HnPHI] = ...
                LOCAL_normal_of_spheroidal_hessian(nu_r, region, solid_harmonic_Gshc(:, :, k));

            if source_precomp.oblate
                [Hxk(idx,:), Hyk(idx,:), Hzk(idx,:)] = convert_from_oblate_basis(HnU, HnV, HnPHI, region.u, region.v, region.phi);
            else
                [Hxk(idx,:), Hyk(idx,:), Hzk(idx,:)] = convert_from_prolate_basis(HnU, HnV, HnPHI, region.u, region.v, region.phi);
            end
        end

        Hn_x(:, :, k) = Hxk;
        Hn_y(:, :, k) = Hyk;
        Hn_z(:, :, k) = Hzk;
    end
end

function [HnU, HnV, HnPHI] = LOCAL_normal_of_spheroidal_hessian(nu_r, region, Gshc)
    [Umat, Vmat, PHImat] = build_directional_hessian_mats(region, nu_r, region.a, region.oblate);

    HnU = Umat * Gshc;
    HnV = Vmat * Gshc;
    HnPHI = PHImat * Gshc;
end

function [res_x, res_y, res_z] = LOCAL_combine_family_traction_terms( ...
        X_eval, nu_eval, dplus_all, dminus_all, dzero_all, ...
        Hn_x_all, Hn_y_all, Hn_z_all, ...
        idx_plus, idx_minus, idx_zero, idx_y_plus, idx_y_minus, idx_y_zero)
    dplus = slice_density(dplus_all, idx_plus) ./ sqrt(2);
    dminus = slice_density(dminus_all, idx_minus) ./ sqrt(2);
    dzero = slice_density(dzero_all, idx_zero);

    xplus = (X_eval(:, 1, :) + 1i .* X_eval(:, 2, :)) ./ sqrt(2);
    xminus = (X_eval(:, 1, :) - 1i .* X_eval(:, 2, :)) ./ sqrt(2);
    xzero = X_eval(:, 3, :);

    nx = nu_eval(:, 1, :);
    ny = nu_eval(:, 2, :);
    nz = nu_eval(:, 3, :);

    res_x = nx .* dplus - xplus .* slice_density(Hn_x_all, idx_plus) + slice_density(Hn_x_all, idx_y_plus);
    res_y = ny .* dplus - xplus .* slice_density(Hn_y_all, idx_plus) + slice_density(Hn_y_all, idx_y_plus);
    res_z = nz .* dplus - xplus .* slice_density(Hn_z_all, idx_plus) + slice_density(Hn_z_all, idx_y_plus);

    res_x = res_x + nx .* dminus - xminus .* slice_density(Hn_x_all, idx_minus) + slice_density(Hn_x_all, idx_y_minus);
    res_y = res_y + ny .* dminus - xminus .* slice_density(Hn_y_all, idx_minus) + slice_density(Hn_y_all, idx_y_minus);
    res_z = res_z + nz .* dminus - xminus .* slice_density(Hn_z_all, idx_minus) + slice_density(Hn_z_all, idx_y_minus);

    res_x = res_x + nx .* dzero - xzero .* slice_density(Hn_x_all, idx_zero) + slice_density(Hn_x_all, idx_y_zero);
    res_y = res_y + ny .* dzero - xzero .* slice_density(Hn_y_all, idx_zero) + slice_density(Hn_y_all, idx_y_zero);
    res_z = res_z + nz .* dzero - xzero .* slice_density(Hn_z_all, idx_zero) + slice_density(Hn_z_all, idx_y_zero);
end

function fine_grid_vec = LOCAL_transfer_to_finer_grid(fine_p, coarse_p, coarse_vec)
    coarse_sp = (coarse_p + 1)^2;
    fine_sp = (fine_p + 1)^2;
    coarse_shc = shAna(coarse_vec);
    padded_shc = [coarse_shc; zeros(fine_sp - coarse_sp, size(coarse_vec, 2))];
    fine_grid_vec = shSyn(padded_shc, isreal(coarse_vec));
end

function coarse_grid_vec = LOCAL_transfer_to_coarser_grid(coarse_p, fine_p, fine_vec)
    coarse_sp = (coarse_p + 1)^2;
    fine_shc = shAna(fine_vec);
    truncated_shc = fine_shc(1:coarse_sp, :);
    coarse_grid_vec = shSyn(truncated_shc, isreal(fine_vec));
end
