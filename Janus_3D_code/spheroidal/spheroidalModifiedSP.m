function modSP = spheroidalModifiedSP(params, lambda, X_trg)
%--------------------------------------------------------------------%
% spheroidalModifiedSP computes the outward-normal derivative of the
% modified Laplace single-layer potential. Note that this is only intended for
% self-evaluation.
%    (o) params:
%       (o) sigma = density on the spheroid surface, as a function of (theta,phi). 
%           See gl_grid. cos(theta) are gauss-Legendre nodes and phi are 
%           equispaced. sigma is size [np,nf,ns] and each column of sigma 
%           corresponds to a different spheroid surface.
%       (o) p = order
%       (o) u0 = 1/eccentricity of the spheroid surface.
%       (o) a = scale of spheroid
%       (o) isReal = Can be set if real output is expected. then imaginary 
%           parts are removed 
%       (o) sigma_coefficients: spherical harmonic coefficients of sigma
%    (o) lambda:
%           Yukawa parameter
%    (o) X = (x,y,z) coordinates of target points. Is a cell of length
%        ns or an nt x 3 x ns matrix.
%    
% Returns modSP as a matrix if X is a matrix or not provided.
%--------------------------------------------------------------------%


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

if any(oblate)
    error('Not implemented.');
end

% Use same weighted ASWF basis transform as modified SLP.
Gswfc = zeros(sp, nf, ns);
Gcache = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k = 1:ns
    c = 1j * lambda * a(k);
    swfc_k = shc_to_swfc(shc(:, :, k), p, oblate(k), c);
    key = u0(k);
    if ~isKey(Gcache, key)
        Gcache(key) = ASWF_Gmatrix(p, u0(k), c, 0, oblate(k), 50, 0);
    end
    Gk = Gcache(key);
    Gswfc(:, :, k) = Gk \ swfc_k;
end

if isempty(X_trg)
    [theta2, ~] = gl_grid(p);
    v_surf = cos(theta2);
    nt = numel(theta2);
    modSP = zeros(nt, nf, ns);

    for k = 1:ns
        c = 1j * lambda * a(k);
        Snm_k = LOCAL_compute_Snm(p, c);
        u_surf = u0(k) .* ones(nt, 1);
        [~, spectra_surf, ~] = LOCAL_modSPspectrum(p, u0(k), oblate(k), c, u_surf, v_surf(:));
        modSP(:, :, k) = (Snm_k .* spectra_surf) * Gswfc(:, :, k);
    end

    if isReal
        modSP = real(modSP);
    end
else
    nt = size(X_trg, 1);
    modSP = zeros(nt, nf, ns);

    for k = 1:ns
        c = 1j * lambda * a(k);

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
            LOCAL_modSPspectrum(p, u0(k), oblate(k), c, u_x, v_x);
        spectra_regions = {spectra_int_k, spectra_surf_k, spectra_ext_k};

        modSPk = zeros(ntk, nf);
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
                An = ASWFnm(n, [], v_x_r, phi_x_r, c, p, 0);
                Sr(:, n^2+1:(n+1)^2) = An;
            end

            modSPk(idx, :) = (Sr .* spectra_r) * Gswfc(:, :, k);
        end

        if isReal
            modSPk = real(modSPk);
        end
        modSP(:, :, k) = modSPk;
    end
end
end % END MAIN FUNCTION

function [spectra_int, spectra_surf, spectra_ext] = LOCAL_modSPspectrum(p, u0, oblate, c, u_x, v_x)
    if oblate, error('Not implemented.'); end

    u_x = u_x(:);
    v_x = real(v_x(:));
    nt = numel(u_x);
    sp = (p + 1)^2;
    spectra_int = zeros(nt, sp);
    spectra_surf = zeros(nt, sp);
    spectra_ext = zeros(nt, sp);
    R1_u0 = zeros(1, sp);
    R3_u0 = zeros(1, sp);
    dR1_u0 = zeros(1, sp);
    dR3_u0 = zeros(1, sp);

    for n = 0:p
        idx = n^2 + 1:(n + 1)^2;
        [R1_n, dR1_n] = Rnm1(n, [], u0, c);
        [R3_n, dR3_n] = Rnm3(n, [], u0, c);
        R1_u0(idx) = reshape(R1_n, 1, []);
        R3_u0(idx) = reshape(R3_n, 1, []);
        dR1_u0(idx) = reshape(dR1_n, 1, []);
        dR3_u0(idx) = reshape(dR3_n, 1, []);
    end
    anm = (1i * c * sqrt(u0^2 - 1)) .* sqrt((u_x.^2 - 1) ./ (u_x.^2 - v_x.^2));

    tol = 9e-12;
    id_int = (u_x < u0 - tol);
    id_ext = (u_x > u0 + tol);
    id_surf = ~id_int & ~id_ext;

    if any(id_int)
        dR1_u = radial_derivative_block_pswf(p, u_x(id_int), c, 1);
        spectra_int(id_int, :) = anm(id_int) .* R3_u0 .* dR1_u;
    end

    if any(id_ext)
        dR3_u = radial_derivative_block_pswf(p, u_x(id_ext), c, 2);
        spectra_ext(id_ext, :) = anm(id_ext) .* R1_u0 .* dR3_u;
    end

    if any(id_surf)
        spectra_surf(id_surf, :) = anm(id_surf) .* 0.5 * (R3_u0 .* dR1_u0 + R1_u0 .* dR3_u0);
    end
end

function S = LOCAL_compute_Snm(p, c)
    % Compute Ynm blocks used in self-eval branch.

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
