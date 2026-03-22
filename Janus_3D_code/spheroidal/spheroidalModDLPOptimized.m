function modDL = spheroidalModDLPOptimized(params, lambda, X_trg, mex_opts, radial_ws)
if nargin < 4
    mex_opts = struct();
end
if nargin < 5
    radial_ws = [];
end
mex_opts = modified_laplace_get_mex_options(mex_opts);
if isempty(radial_ws)
    radial_ws = spheroidal_modified_radial_workspace(params, lambda, X_trg, mex_opts);
end

p = params.p;
u0 = params.u0;
a = params.a;
isReal = params.isReal;
oblate = params.oblate;
shc = params.sigma_coefficients;
[sp, nf, ns] = size(shc);

if sp ~= (p + 1)^2
    error("shc first dimension must be (p+1)^2.");
end

if isscalar(u0)
    u0 = u0 .* ones(1, ns);
end
if isscalar(a)
    a = a .* ones(1, ns);
end
if isscalar(oblate)
    oblate = oblate .* ones(1, ns);
end

if isempty(X_trg)
    [theta2, ~] = gl_grid(p);
    nt = numel(theta2);
    modDL = zeros(nt, nf, ns);

    for k = 1:ns
        c_ang_k = radial_ws.c_ang(k);
        swfc_k = shc_to_swfc(shc(:, :, k), p, c_ang_k);
        Snm_k = LOCAL_compute_Snm(p, c_ang_k);

        target_meta = radial_ws.targets(k);
        group = radial_ws.groups(target_meta.group_index);
        surface_group_u_index = target_meta.surface_group_u_index;

        spectra_surf_k = LOCAL_get_modDLP_surface_spectrum(u0(k), oblate(k), radial_ws.c_rad(k), group, surface_group_u_index);
        spectra_matrix = repmat(spectra_surf_k(:), 1, nf);
        modDL(:, :, k) = Snm_k * (spectra_matrix .* swfc_k);
    end

    if isReal
        modDL = real(modDL);
    end
    return;
end

nt = size(X_trg, 1);
modDL = zeros(nt, nf, ns);

for k = 1:ns
    c_ang_k = radial_ws.c_ang(k);
    swfc_k = shc_to_swfc(shc(:, :, k), p, c_ang_k);

    target_meta = radial_ws.targets(k);
    ntk = target_meta.nt;
    if ntk == 0
        continue;
    end

        group = radial_ws.groups(target_meta.group_index);
        surface_group_u_index = target_meta.surface_group_u_index;

        [spectra_int_k, spectra_surf_k, spectra_ext_k] = LOCAL_get_modDLP_spectra(u0(k), oblate(k), radial_ws.c_rad(k), group, surface_group_u_index);
        spectra_regions = {spectra_int_k(:), spectra_surf_k(:), spectra_ext_k(:)};
        regions = {target_meta.indices_interior, target_meta.indices_surface, target_meta.indices_exterior};

    modDLk = zeros(ntk, nf);
    for r = 1:3
        idx = regions{r};
        if ~any(idx)
            continue;
        end

        v_x_r = target_meta.v_x(idx);
        phi_x_r = target_meta.phi_x(idx);
        target_group_u_indices_r = target_meta.target_group_u_indices(idx);
        nt_r = nnz(idx);

        Fr = ones(nt_r, sp);
        switch r
            case 1
                Fr = group.R1(target_group_u_indices_r, :);
            case 3
                Fr = group.R3(target_group_u_indices_r, :);
        end

        Sr = zeros(nt_r, sp);
        for n = 0:p
            An = ASWFnm(n, [], v_x_r, phi_x_r, c_ang_k, p, 0);
            Sr(:, n^2+1:(n+1)^2) = An;
        end

        SYr = Fr .* Sr;
        spectra_matrix = repmat(spectra_regions{r}, 1, nf);
        modDLcoefs_r = spectra_matrix .* swfc_k;
        modDLk(idx, :) = SYr * modDLcoefs_r;
    end

    if isReal
        modDLk = real(modDLk);
    end
    modDL(:, :, k) = modDLk;
end
end

function [lambda_int, lambda_surf, lambda_ext] = LOCAL_get_modDLP_spectra(u0, oblate, c_rad, group, sr)
    if oblate
        cnm = -c_rad * (u0^2 + 1);
    else
        cnm = 1j * c_rad * (u0^2 - 1);
    end

    R1_u0 = group.R1(sr, :);
    dR1_u0 = group.dR1(sr, :);
    R3_u0 = group.R3(sr, :);
    dR3_u0 = group.dR3(sr, :);

    lambda_int = cnm .* dR3_u0;
    lambda_surf = cnm .* (dR1_u0 .* R3_u0 + dR3_u0 .* R1_u0) / 2;
    lambda_ext = cnm .* dR1_u0;
end

function lambda_surf = LOCAL_get_modDLP_surface_spectrum(u0, oblate, c_rad, group, sr)
    [~, lambda_surf, ~] = LOCAL_get_modDLP_spectra(u0, oblate, c_rad, group, sr);
end

function S = LOCAL_compute_Snm(p, c)
    [theta2, phi_k] = gl_grid(p);
    v_k = cos(theta2);
    v_k = v_k(:);
    phi_k = phi_k(:);
    nt_r = length(v_k);
    sp = (p + 1)^2;

    S = zeros(nt_r, sp);
    v_row = real(v_k);
    for n = 0:p
        Yn = ASWFnm(n, [], v_row, phi_k, c, p, 0);
        S(:, n^2+1:(n+1)^2) = Yn;
    end
end
