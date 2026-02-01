function [res_x, res_y, res_z] = L2StkTLPOptimized(X_trg, nu_trg, params, sigma_x, sigma_y, sigma_z, source_Gmatrix, alternate_spp_flag)
%{
Optimized version of L2StkTLP for the Stokes traction layer potential.
%}

if nargin < 8
    alternate_spp_flag = false;
end

[np, nf, ns] = size(sigma_x);

% Fallback to the reference implementation for the alternate formula.
if alternate_spp_flag
    if isempty(X_trg)
        X_cell = [];
        if isempty(nu_trg)
            nu_cell = cell(1, ns);
            for i = 1:ns
                nu_cell{i} = get_norm_vecs(params.p, params.u0(i), params.oblate(i));
            end
        else
            if ~isequal(size(nu_trg), [np 3 ns])
                error("nu_trg must be np x 3 x ns for self-eval mode.");
            end
            nu_cell = squeeze(num2cell(nu_trg, [1, 2])).';
        end
    else
        X_cell = squeeze(num2cell(X_trg, [1, 2])).';
        nu_cell = squeeze(num2cell(nu_trg, [1, 2])).';
    end

    [cx, cy, cz] = L2StkTLP(X_cell, nu_cell, params, sigma_x, sigma_y, sigma_z, ns, true);
    if iscell(cx)
        res_x = cat(3, cx{:});
        res_y = cat(3, cy{:});
        res_z = cat(3, cz{:});
    else
        res_x = cx;
        res_y = cy;
        res_z = cz;
    end
    return;
end

p = params.p;
u0 = params.u0;
a = params.a;
isReal = params.isReal;
oblate = params.oblate;

if isempty(X_trg)
    if ~isempty(nu_trg)
        if ~isequal(size(nu_trg), [np 3 ns])
            error("nu_trg must be np x 3 x ns for self-eval mode.");
        end
    end
else
    if ~isnumeric(X_trg) || size(X_trg, 2) ~= 3 || size(X_trg, 3) ~= ns
        error("X_trg must be nt x 3 x ns for numeric-only mode.");
    end
    if isempty(nu_trg) || ~isnumeric(nu_trg) || size(nu_trg, 2) ~= 3 || size(nu_trg, 1) ~= size(X_trg, 1) || size(nu_trg, 3) ~= ns
        error("nu_trg must be nt x 3 x ns to match X_trg.");
    end
end

% Compute y_j * sigma on the source surface (local frame)
Xloc_all = zeros(np, 3, ns);
y1x = zeros(np, nf, ns);
y1y = zeros(np, nf, ns);
y1z = zeros(np, nf, ns);
y2x = zeros(np, nf, ns);
y2y = zeros(np, nf, ns);
y2z = zeros(np, nf, ns);
y3x = zeros(np, nf, ns);
y3y = zeros(np, nf, ns);
y3z = zeros(np, nf, ns);

for i = 1:ns
    if params.oblate(i)
        Xloc_i = oblate_spheroid_shape(p, params.u0(i), params.a(i));
    else
        Xloc_i = prolate_spheroid_shape(p, params.u0(i), params.a(i));
    end
    Xloc_all(:,:,i) = Xloc_i;

    y1x(:,:,i) = sigma_x(:,:,i) .* Xloc_i(:, 1);
    y1y(:,:,i) = sigma_y(:,:,i) .* Xloc_i(:, 1);
    y1z(:,:,i) = sigma_z(:,:,i) .* Xloc_i(:, 1);

    y2x(:,:,i) = sigma_x(:,:,i) .* Xloc_i(:, 2);
    y2y(:,:,i) = sigma_y(:,:,i) .* Xloc_i(:, 2);
    y2z(:,:,i) = sigma_z(:,:,i) .* Xloc_i(:, 2);

    y3x(:,:,i) = sigma_x(:,:,i) .* Xloc_i(:, 3);
    y3y(:,:,i) = sigma_y(:,:,i) .* Xloc_i(:, 3);
    y3z(:,:,i) = sigma_z(:,:,i) .* Xloc_i(:, 3);
end

% Build normals if not provided
if isempty(X_trg)
    if isempty(nu_trg)
        nu_trg_all = zeros(np, 3, ns);
        for i = 1:ns
            nu_trg_all(:,:,i) = get_norm_vecs(p, u0(i), oblate(i));
        end
    else
        nu_trg_all = nu_trg;
    end
else
    nu_trg_all = nu_trg;
end

% Compute SH coefficients and Gmatrix transforms for sigma only
[shc_x, shc_y, shc_z] = generate_shc(sigma_x, sigma_y, sigma_z);
[Gshc_x, Gshc_y, Gshc_z] = generate_mod_shc(source_Gmatrix, shc_x, shc_y, shc_z);
Gshc_all = cat(2, Gshc_x, Gshc_y, Gshc_z);

% Compute SH coefficients and Gmatrix transforms for y_j * sigma
[shc_y1x, shc_y1y, shc_y1z, shc_y2x, shc_y2y, shc_y2z, shc_y3x, shc_y3y, shc_y3z] = generate_shc( ...
    y1x, y1y, y1z, y2x, y2y, y2z, y3x, y3y, y3z ...
);
[Gshc_y1x, Gshc_y1y, Gshc_y1z, Gshc_y2x, Gshc_y2y, Gshc_y2z, Gshc_y3x, Gshc_y3y, Gshc_y3z] = generate_mod_shc( ...
    source_Gmatrix, shc_y1x, shc_y1y, shc_y1z, shc_y2x, shc_y2y, shc_y2z, shc_y3x, shc_y3y, shc_y3z ...
);

% Build unit vectors for SP
[nu_x, nu_y, nu_z] = build_unit_nu(X_trg, np, ns);

% Evaluate SP in one pass
[SPx_all, SPy_all, SPz_all] = spheroidalSPOptimized(p, u0, a, oblate, Gshc_all, isReal, nu_x, nu_y, nu_z, X_trg);

idx1 = 1:nf;
idx2 = nf + (1:nf);
idx3 = 2*nf + (1:nf);

SP_sigmax_X = slice_density(SPx_all, idx1);
SP_sigmay_X = slice_density(SPx_all, idx2);
SP_sigmaz_X = slice_density(SPx_all, idx3);

SP_sigmax_Y = slice_density(SPy_all, idx1);
SP_sigmay_Y = slice_density(SPy_all, idx2);
SP_sigmaz_Y = slice_density(SPy_all, idx3);

SP_sigmax_Z = slice_density(SPz_all, idx1);
SP_sigmay_Z = slice_density(SPz_all, idx2);
SP_sigmaz_Z = slice_density(SPz_all, idx3);

% Compute grad-div SL terms (S'')
[gd_sig_U, gd_sig_V, gd_sig_PHI] = spheroidalgraddivSLOptimized(p, u0, a, oblate, Gshc_x, Gshc_y, Gshc_z, isReal, X_trg);
[gd_y1_U, gd_y1_V, gd_y1_PHI] = spheroidalgraddivSLOptimized(p, u0, a, oblate, Gshc_y1x, Gshc_y1y, Gshc_y1z, isReal, X_trg);
[gd_y2_U, gd_y2_V, gd_y2_PHI] = spheroidalgraddivSLOptimized(p, u0, a, oblate, Gshc_y2x, Gshc_y2y, Gshc_y2z, isReal, X_trg);
[gd_y3_U, gd_y3_V, gd_y3_PHI] = spheroidalgraddivSLOptimized(p, u0, a, oblate, Gshc_y3x, Gshc_y3y, Gshc_y3z, isReal, X_trg);

if isempty(X_trg)
    nt = np;
else
    nt = size(X_trg, 1);
end

gd_sig_x = zeros(nt, nf, ns);
gd_sig_y = zeros(nt, nf, ns);
gd_sig_z = zeros(nt, nf, ns);

gd_y1_x = zeros(nt, nf, ns);
gd_y1_y = zeros(nt, nf, ns);
gd_y1_z = zeros(nt, nf, ns);

gd_y2_x = zeros(nt, nf, ns);
gd_y2_y = zeros(nt, nf, ns);
gd_y2_z = zeros(nt, nf, ns);

gd_y3_x = zeros(nt, nf, ns);
gd_y3_y = zeros(nt, nf, ns);
gd_y3_z = zeros(nt, nf, ns);

for k = 1:ns
    if isempty(X_trg)
        Xk = Xloc_all(:,:,k);
    else
        Xk = X_trg(:,:,k);
    end

    S_eval_k = cart2spheroidal(Xk, a(k), oblate(k));
    [nu_x_sph, ~] = cartNu2spheroidal(nu_x(:,:,k), S_eval_k, a(k), oblate(k));
    [nu_y_sph, ~] = cartNu2spheroidal(nu_y(:,:,k), S_eval_k, a(k), oblate(k));
    [nu_z_sph, ~] = cartNu2spheroidal(nu_z(:,:,k), S_eval_k, a(k), oblate(k));

    gd_sig_x(:,:,k) = nu_x_sph(:, 1) .* gd_sig_U(:,:,k) + nu_x_sph(:, 2) .* gd_sig_V(:,:,k) + nu_x_sph(:, 3) .* gd_sig_PHI(:,:,k);
    gd_sig_y(:,:,k) = nu_y_sph(:, 1) .* gd_sig_U(:,:,k) + nu_y_sph(:, 2) .* gd_sig_V(:,:,k) + nu_y_sph(:, 3) .* gd_sig_PHI(:,:,k);
    gd_sig_z(:,:,k) = nu_z_sph(:, 1) .* gd_sig_U(:,:,k) + nu_z_sph(:, 2) .* gd_sig_V(:,:,k) + nu_z_sph(:, 3) .* gd_sig_PHI(:,:,k);

    gd_y1_x(:,:,k) = nu_x_sph(:, 1) .* gd_y1_U(:,:,k) + nu_x_sph(:, 2) .* gd_y1_V(:,:,k) + nu_x_sph(:, 3) .* gd_y1_PHI(:,:,k);
    gd_y1_y(:,:,k) = nu_y_sph(:, 1) .* gd_y1_U(:,:,k) + nu_y_sph(:, 2) .* gd_y1_V(:,:,k) + nu_y_sph(:, 3) .* gd_y1_PHI(:,:,k);
    gd_y1_z(:,:,k) = nu_z_sph(:, 1) .* gd_y1_U(:,:,k) + nu_z_sph(:, 2) .* gd_y1_V(:,:,k) + nu_z_sph(:, 3) .* gd_y1_PHI(:,:,k);

    gd_y2_x(:,:,k) = nu_x_sph(:, 1) .* gd_y2_U(:,:,k) + nu_x_sph(:, 2) .* gd_y2_V(:,:,k) + nu_x_sph(:, 3) .* gd_y2_PHI(:,:,k);
    gd_y2_y(:,:,k) = nu_y_sph(:, 1) .* gd_y2_U(:,:,k) + nu_y_sph(:, 2) .* gd_y2_V(:,:,k) + nu_y_sph(:, 3) .* gd_y2_PHI(:,:,k);
    gd_y2_z(:,:,k) = nu_z_sph(:, 1) .* gd_y2_U(:,:,k) + nu_z_sph(:, 2) .* gd_y2_V(:,:,k) + nu_z_sph(:, 3) .* gd_y2_PHI(:,:,k);

    gd_y3_x(:,:,k) = nu_x_sph(:, 1) .* gd_y3_U(:,:,k) + nu_x_sph(:, 2) .* gd_y3_V(:,:,k) + nu_x_sph(:, 3) .* gd_y3_PHI(:,:,k);
    gd_y3_y(:,:,k) = nu_y_sph(:, 1) .* gd_y3_U(:,:,k) + nu_y_sph(:, 2) .* gd_y3_V(:,:,k) + nu_y_sph(:, 3) .* gd_y3_PHI(:,:,k);
    gd_y3_z(:,:,k) = nu_z_sph(:, 1) .* gd_y3_U(:,:,k) + nu_z_sph(:, 2) .* gd_y3_V(:,:,k) + nu_z_sph(:, 3) .* gd_y3_PHI(:,:,k);
end

res_x = zeros(nt, nf, ns);
res_y = zeros(nt, nf, ns);
res_z = zeros(nt, nf, ns);

for k = 1:ns
    if isempty(X_trg)
        Xk = Xloc_all(:,:,k);
    else
        Xk = X_trg(:,:,k);
    end

    nu_k = nu_trg_all(:,:,k);
    nx_trg = nu_k(:, 1);
    ny_trg = nu_k(:, 2);
    nz_trg = nu_k(:, 3);
    n_dot_x = dot(Xk, nu_k, 2);

    SP_sum_term_x = nx_trg .* SP_sigmax_X(:,:,k) + ny_trg .* SP_sigmax_Y(:,:,k) + nz_trg .* SP_sigmax_Z(:,:,k);
    SP_sum_term_y = nx_trg .* SP_sigmay_X(:,:,k) + ny_trg .* SP_sigmay_Y(:,:,k) + nz_trg .* SP_sigmay_Z(:,:,k);
    SP_sum_term_z = nx_trg .* SP_sigmaz_X(:,:,k) + ny_trg .* SP_sigmaz_Y(:,:,k) + nz_trg .* SP_sigmaz_Z(:,:,k);

    graddivSL_sum_term_X = nx_trg .* gd_y1_x(:,:,k) + ny_trg .* gd_y2_x(:,:,k) + nz_trg .* gd_y3_x(:,:,k);
    graddivSL_sum_term_Y = nx_trg .* gd_y1_y(:,:,k) + ny_trg .* gd_y2_y(:,:,k) + nz_trg .* gd_y3_y(:,:,k);
    graddivSL_sum_term_Z = nx_trg .* gd_y1_z(:,:,k) + ny_trg .* gd_y2_z(:,:,k) + nz_trg .* gd_y3_z(:,:,k);

    res_x(:,:,k) = SP_sum_term_x + graddivSL_sum_term_X - n_dot_x .* gd_sig_x(:,:,k);
    res_y(:,:,k) = SP_sum_term_y + graddivSL_sum_term_Y - n_dot_x .* gd_sig_y(:,:,k);
    res_z(:,:,k) = SP_sum_term_z + graddivSL_sum_term_Z - n_dot_x .* gd_sig_z(:,:,k);
end

if isReal
    res_x = real(res_x);
    res_y = real(res_y);
    res_z = real(res_z);
end
end %% END MAIN FUNCTION