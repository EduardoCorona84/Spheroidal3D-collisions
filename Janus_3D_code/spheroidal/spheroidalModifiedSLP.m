function modSL = spheroidalModifiedSLP(params, lambda, X_trg, mex_opts)
%--------------------------------------------------------------------%
% spheroidalModifiedSLP computes the modified Laplace single-layer potential.
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
% Returns modSL as a matrix if X is a matrix or not provided.
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

% Convert to modified ASWF basis.
Gswfc = zeros(sp, nf, ns);
for k = 1:ns
    c_ang_k = c_ang(k);
    swfc_k = shc_to_swfc(shc(:, :, k), p, c_ang_k);
    Gk = ASWF_Gmatrix(p, u0(k), c_ang_k, 0, oblate(k), 50, 0);
    Gswfc(:, :, k) = Gk \ swfc_k;
end

if isempty(X_trg)
    [theta2, ~] = gl_grid(p);
    nt = numel(theta2);
    modSL = zeros(nt, nf, ns);

    for k = 1:ns
        c_ang_k = c_ang(k);
        c_rad_k = c_rad(k);
        Snm_k = LOCAL_compute_Snm(p, c_ang_k);
        [~, spectra_surf_k, ~] = LOCAL_modSLPspectrum(p, u0(k), a(k), oblate(k), c_rad_k, mex_opts);
        spectra_matrix = repmat(spectra_surf_k(:), 1, nf);
        modSL(:, :, k) = Snm_k * (spectra_matrix .* Gswfc(:, :, k));
    end

    if isReal
        modSL = real(modSL);
    end
else
    nt = size(X_trg, 1);
    modSL = zeros(nt, nf, ns);

    for k = 1:ns
        c_ang_k = c_ang(k);
        c_rad_k = c_rad(k);
        [spectra_int_k, spectra_surf_k, spectra_ext_k] = LOCAL_modSLPspectrum(p, u0(k), a(k), oblate(k), c_rad_k, mex_opts);

        Xtk = X_trg(:, :, k);
        ntk = size(Xtk, 1);
        if ntk == 0
            continue;
        end

        S = cart2spheroidal(Xtk, a(k), oblate(k));
        u_x = S(:, 1);

        indices_interior = (u_x < u0(k) - 9e-12);
        indices_surface = (abs(u_x - u0(k)) <= 9e-12);
        indices_exterior = (u_x > u0(k) + 9e-12);

        regions = {indices_interior, indices_surface, indices_exterior};
        spectra_regions = {spectra_int_k(:), spectra_surf_k(:), spectra_ext_k(:)};

        modSLk = zeros(ntk, nf);
        for r = 1:3
            idx = regions{r};
            if ~any(idx)
                continue;
            end

            Sr = S(idx, :);
            u_x_r = Sr(:, 1);
            v_x_r = Sr(:, 2);
            phi_x_r = Sr(:, 3);

            Fr = LOCAL_solid_swf(p, u0(k), u_x_r, oblate(k), c_rad_k, mex_opts);
            nt_r = length(u_x_r);
            Sr = zeros(nt_r, sp);
            v_row = real(v_x_r);

            for n = 0:p
                An = ASWFnm(n, [], v_row, phi_x_r, c_ang_k, p, 0);
                Sr(:, n^2+1:(n+1)^2) = An;
            end

            SYr = Fr .* Sr;
            spectra_matrix = repmat(spectra_regions{r}, 1, nf);
            modSLcoefs_r = spectra_matrix .* Gswfc(:, :, k);
            modSLk(idx, :) = SYr * modSLcoefs_r;
        end

        if isReal
            modSLk = real(modSLk);
        end
        modSL(:, :, k) = modSLk;
    end
end
end %% END MAIN FUNCTION

function [lambda_int, lambda_surf, lambda_ext] = LOCAL_modSLPspectrum(p, u0, a, oblate, c, mex_opts)

    if oblate
        cnm = -a * c * sqrt(u0^2 + 1);
    else
        cnm = 1j * a * c * sqrt(u0^2 - 1);
    end

    [Rnm1_u0, ~, Rnm3_u0, ~] = modifiedLaplaceEvalRadialSphwv(p, u0, c, oblate, mex_opts);

    lambda_int = cnm .* Rnm3_u0;
    lambda_surf = cnm .* (Rnm1_u0 .* Rnm3_u0);
    lambda_ext = cnm .* Rnm1_u0;
end

function F = LOCAL_solid_swf(p, u0, u_x, oblate, c, mex_opts)
    F = ones(size(u_x, 1), (p + 1)^2);
    idx_int = abs(u_x) - u0 < -1e-14;
    idx_ext = abs(u_x) - u0 > 1e-14;

    if any(idx_int)
        [R1, ~, ~, ~] = modifiedLaplaceEvalRadialSphwv(p, u_x(idx_int), c, oblate, mex_opts);
        F(idx_int, :) = R1;
    end
    if any(idx_ext)
        [~, ~, R3, ~] = modifiedLaplaceEvalRadialSphwv(p, u_x(idx_ext), c, oblate, mex_opts);
        F(idx_ext, :) = R3;
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
