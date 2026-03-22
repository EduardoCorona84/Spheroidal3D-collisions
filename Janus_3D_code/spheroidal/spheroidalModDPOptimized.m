function modDP = spheroidalModDPOptimized(params, lambda, X_trg, mex_opts, radial_ws)
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
    v_surf = cos(theta2);
    nt = numel(theta2);
    modDP = zeros(nt, nf, ns);

    for k = 1:ns
        c_ang_k = radial_ws.c_ang(k);
        swfc_k = shc_to_swfc(shc(:, :, k), p, c_ang_k);
        Snm_k = LOCAL_compute_Snm(p, c_ang_k);

        target_meta = radial_ws.targets(k);
        group = radial_ws.groups(target_meta.group_index);
        surface_group_u_index = target_meta.surface_group_u_index;
        u_surf = u0(k) .* ones(nt, 1);
        spectra_surf = LOCAL_get_modDP_surface_spectrum(p, u0(k), a(k), oblate(k), radial_ws.c_rad(k), group, surface_group_u_index, u_surf, v_surf(:));
        modDP(:, :, k) = (Snm_k .* spectra_surf) * swfc_k;
    end

    if isReal
        modDP = real(modDP);
    end
    return;
end

nt = size(X_trg, 1);
modDP = zeros(nt, nf, ns);

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

        [spectra_int_k, spectra_surf_k, spectra_ext_k] = LOCAL_get_modDP_spectra( ...
            p, u0(k), a(k), oblate(k), radial_ws.c_rad(k), group, surface_group_u_index, ...
            target_meta.u_x, target_meta.v_x, target_meta.target_group_u_indices, ...
            target_meta.indices_interior, target_meta.indices_surface, target_meta.indices_exterior);
    spectra_regions = {spectra_int_k, spectra_surf_k, spectra_ext_k};
    regions = {target_meta.indices_interior, target_meta.indices_surface, target_meta.indices_exterior};

    modDPk = zeros(ntk, nf);
    for r = 1:3
        idx = regions{r};
        if ~any(idx)
            continue;
        end

        v_x_r = target_meta.v_x(idx);
        phi_x_r = target_meta.phi_x(idx);
        spectra_r = spectra_regions{r};
        spectra_r = spectra_r(idx, :);
        nt_r = nnz(idx);

        Sr = zeros(nt_r, sp);
        for n = 0:p
            An = ASWFnm(n, [], v_x_r, phi_x_r, c_ang_k, p, 0);
            Sr(:, n^2+1:(n+1)^2) = An;
        end

        modDPk(idx, :) = (Sr .* spectra_r) * swfc_k;
    end

    if isReal
        modDPk = real(modDPk);
    end
    modDP(:, :, k) = modDPk;
end
end %% END MAIN FUNCTION

function spectra_surf = LOCAL_get_modDP_surface_spectrum(p, u0, a, oblate, c, group, surface_group_u_index, u_x, v_x)
    [~, spectra_surf, ~] = LOCAL_get_modDP_spectra(p, u0, a, oblate, c, group, surface_group_u_index, u_x, v_x, repmat(surface_group_u_index, numel(u_x), 1), false(size(u_x)), true(size(u_x)), false(size(u_x)));
end

function [spectra_int, spectra_surf, spectra_ext] = LOCAL_get_modDP_spectra(p, u0, a, oblate, c, group, surface_group_u_index, u_x, v_x, target_group_u_indices, id_int, id_surf, id_ext)
    u_x = u_x(:);
    v_x = real(v_x(:));
    target_group_u_indices = target_group_u_indices(:);
    nt = numel(u_x);
    sp = (p + 1)^2;
    spectra_int = zeros(nt, sp);
    spectra_surf = zeros(nt, sp);
    spectra_ext = zeros(nt, sp);

    dR1_u0 = group.dR1(surface_group_u_index, :);
    dR3_u0 = group.dR3(surface_group_u_index, :);

    if oblate
        anm = -c*(u0^2 + 1)/a .* sqrt((u_x.^2 + 1) ./ (u_x.^2 + v_x.^2));
    else
        anm = 1i*c*(u0^2 - 1)/a .* sqrt((u_x.^2 - 1) ./ (u_x.^2 - v_x.^2));
    end

    if any(id_int)
        spectra_int(id_int, :) = anm(id_int) .* dR3_u0 .* group.dR1(target_group_u_indices(id_int), :);
    end
    if any(id_ext)
        spectra_ext(id_ext, :) = anm(id_ext) .* dR1_u0 .* group.dR3(target_group_u_indices(id_ext), :);
    end
    if any(id_surf)
        spectra_surf(id_surf, :) = anm(id_surf) .* dR1_u0 .* dR3_u0;
    end
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
