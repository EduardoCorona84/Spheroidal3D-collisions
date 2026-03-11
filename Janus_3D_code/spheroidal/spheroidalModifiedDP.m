function modDP = spheroidalModifiedDP(params, lambda, X_trg, mex_opts)
%--------------------------------------------------------------------%
% spheroidalModifiedDP computes the outward-normal derivative of the
% modified Laplace double-layer potential. This is intended for
% self-evaluation.
%    (o) params:
%       (o) sigma = density on the spheroid surface, as a function of (theta,phi).
%           See gl_grid. cos(theta) are Gauss-Legendre nodes and phi are
%           equispaced. sigma is size [np,nf,ns] and each column of sigma
%           corresponds to a different spheroid surface.
%       (o) p = order
%       (o) u0 = 1/eccentricity of the spheroid surface.
%       (o) a = scale of spheroid
%       (o) isReal = Can be set if real output is expected. Then imaginary
%           parts are removed.
%       (o) sigma_coefficients: spherical harmonic coefficients of sigma
%    (o) lambda:
%           Yukawa parameter
%    (o) X = (x,y,z) coordinates of target points. Is a cell of length
%        ns or an nt x 3 x ns matrix.
%    (o) mex_opts (optional): sphwv MEX controls.
%
% Returns modDP as a matrix if X is a matrix or not provided.
%--------------------------------------------------------------------%

p = params.p;
u0 = params.u0;
a = params.a;
isReal = params.isReal;
oblate = params.oblate;
if nargin < 4
    mex_opts = struct();
end
mex_opts = modifiedLaplaceGetMexOptions(mex_opts);

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

% Setup angular/radial separation parameters.
% For oblate:
%   angular side (ASWFnm / transforms) uses c_ang = c_real
%   radial side uses c_rad = c_real on sphwv pure-imag backend.
% For prolate we keep the existing pure-imag convention.
c_ang = zeros(size(a));
c_rad = zeros(size(a));
for k = 1:numel(a)
    if oblate(k)
        c_rad(k) = lambda * a(k);
        c_ang(k) = c_rad(k);
    else
        c_rad(k) = 1j * lambda * a(k);
        c_ang(k) = c_rad(k);
    end
end

if isempty(X_trg)
    [theta2, ~] = gl_grid(p);
    v_surf = cos(theta2);
    nt = numel(theta2);
    modDP = zeros(nt, nf, ns);

    for k = 1:ns
        c_ang_k = c_ang(k);
        c_rad_k = c_rad(k);
        swfc_k = shc_to_swfc(shc(:, :, k), p, c_ang_k);
        Snm_k = LOCAL_compute_Snm(p, c_ang_k);
        u_surf = u0(k) .* ones(nt, 1);

        [~, spectra_surf, ~] = ...
            LOCAL_modDPspectrum(p, u0(k), a(k), oblate(k), c_rad_k, u_surf, v_surf(:), mex_opts);
        modDP(:, :, k) = (Snm_k .* spectra_surf) * swfc_k;
    end

    if isReal
        modDP = real(modDP);
    end
else
    nt = size(X_trg, 1);
    modDP = zeros(nt, nf, ns);

    for k = 1:ns
        c_ang_k = c_ang(k);
        c_rad_k = c_rad(k);
        swfc_k = shc_to_swfc(shc(:, :, k), p, c_ang_k);

        Xtk = X_trg(:, :, k);
        ntk = size(Xtk, 1);
        if ntk == 0
            continue;
        end

        S = cart2spheroidal(Xtk, a(k), oblate(k));
        u_x = S(:, 1);
        v_x = real(S(:, 2));
        phi_x = S(:, 3);

        indices_interior = (u_x < u0(k) - 9e-12);
        indices_surface = (abs(u_x - u0(k)) <= 9e-12);
        indices_exterior = (u_x > u0(k) + 9e-12);

        regions = {indices_interior, indices_surface, indices_exterior};
        [spectra_int_k, spectra_surf_k, spectra_ext_k] = ...
            LOCAL_modDPspectrum(p, u0(k), a(k), oblate(k), c_rad_k, u_x, v_x, mex_opts);
        spectra_regions = {spectra_int_k, spectra_surf_k, spectra_ext_k};

        modDPk = zeros(ntk, nf);
        for r = 1:3
            idx = regions{r};
            if ~any(idx)
                continue;
            end

            u_x_r = u_x(idx);
            v_x_r = v_x(idx);
            phi_x_r = phi_x(idx);
            nt_r = numel(u_x_r);

            spectra_r = spectra_regions{r};
            spectra_r = spectra_r(idx, :);

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
end
end % END MAIN FUNCTION

function [spectra_int, spectra_surf, spectra_ext] = LOCAL_modDPspectrum(p, u0, a, oblate, c, u_x, v_x, mex_opts)

    u_x = u_x(:);
    v_x = real(v_x(:));
    nt = numel(u_x);
    sp = (p + 1)^2;
    spectra_int = zeros(nt, sp);
    spectra_surf = zeros(nt, sp);
    spectra_ext = zeros(nt, sp);

    [~, dR1_u0, ~, dR3_u0] = modifiedLaplaceEvalRadialSphwv(p, u0, c, oblate, mex_opts);

    if oblate
        % In the oblate modified-Laplace convention, c is real while
        % radial functions are evaluated on cc = 1i*c via sphwv pure-imag.
        % The 1i*c produces a -1i phase versus Kernel_Eval.
        % Multiplying by +1i to align conventions gives -c.
        anm = (-c * (u0^2 + 1) / a) .* sqrt((u_x.^2 + 1) ./ (u_x.^2 + v_x.^2));
    else
        anm = (1i * c * (u0^2 - 1) / a) .* sqrt((u_x.^2 - 1) ./ (u_x.^2 - v_x.^2));
    end

    tol = 9e-12;
    id_int = (u_x < u0 - tol);
    id_ext = (u_x > u0 + tol);
    id_surf = ~id_int & ~id_ext;

    if any(id_int)
        [~, dR1_u, ~, ~] = modifiedLaplaceEvalRadialSphwv(p, u_x(id_int), c, oblate, mex_opts);
        spectra_int(id_int, :) = anm(id_int) .* dR3_u0 .* dR1_u;
    end

    if any(id_ext)
        [~, ~, ~, dR3_u] = modifiedLaplaceEvalRadialSphwv(p, u_x(id_ext), c, oblate, mex_opts);
        spectra_ext(id_ext, :) = anm(id_ext) .* dR1_u0 .* dR3_u;
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
