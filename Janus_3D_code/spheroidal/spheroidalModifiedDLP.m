function modDL = spheroidalModifiedDLP(params, lambda, X_trg)
%--------------------------------------------------------------------%
% spheroidalModifiedDLP computes the modified Laplace layer potential (or the double
% layer potential that is associated with Yukawa's equation).
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
%    (o) X = (x,y,z) coordinates of target points. Is a cell of length
%        ns or an nt x 3 x ns matrix.
%    
% Returns modDL as a matrix if X is a matrix or not provided.
%--------------------------------------------------------------------%

p=params.p;
u0=params.u0;
a=params.a;
isReal=params.isReal;
oblate=params.oblate;

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
        gamma_k = 1j * lambda * a(k);
        swfc_k = shc_to_swfc(shc(:, :, k), p, oblate(k), gamma_k);
        Snm_k = LOCAL_compute_Snm(p, gamma_k);
        [~, spectra_surf_k, ~] = LOCAL_modDLPspectrum(p, u0(k), a(k), oblate(k), gamma_k);
        spectra_matrix = repmat(spectra_surf_k(:), 1, nf);
        modDL(:, :, k) = Snm_k * (spectra_matrix .* swfc_k);
    end

    if isReal
        modDL = real(modDL);
    end
else
    nt = size(X_trg, 1);
    modDL = zeros(nt, nf, ns);

    for k = 1:ns
        gamma_k = 1j * lambda * a(k);
        swfc_k = shc_to_swfc(shc(:, :, k), p, oblate(k), gamma_k);
        [spectra_int_k, spectra_surf_k, spectra_ext_k] = ...
            LOCAL_modDLPspectrum(p, u0(k), a(k), oblate(k), gamma_k);

        Xtk = X_trg(:, :, k);
        ntk = size(Xtk, 1);
        if ntk == 0, continue; end

        S = cart2spheroidal(Xtk, a(k), oblate(k));
        u_x = S(:, 1);

        indices_interior = (u_x < u0(k) - 9e-12);
        indices_surface = (abs(u_x - u0(k)) <= 9e-12);
        indices_exterior = (u_x > u0(k) + 9e-12);

        regions = {indices_interior, indices_surface, indices_exterior};
        spectra_regions = {spectra_int_k(:), spectra_surf_k(:), spectra_ext_k(:)};

        modDLk = zeros(ntk, nf);
        for r = 1:3
            idx = regions{r};
            if ~any(idx), continue; end

            Sr = S(idx, :);
            u_x_r = Sr(:, 1);
            v_x_r = Sr(:, 2);
            phi_x_r = Sr(:, 3);

            Fr = LOCAL_solid_swf(p, u0(k), u_x_r, oblate(k), gamma_k);
            nt_r = length(u_x_r);
            Sr = zeros(nt_r, sp);
            v_row = real(v_x_r);

            for n = 0:p
                An = ASWFnm(n, [], v_row, phi_x_r, gamma_k, p, 0);
                Sr(:, n^2+1:(n+1)^2) = An;
            end

            SYr = Fr .* Sr;
            spectra_matrix = repmat(spectra_regions{r}, 1, nf);
            modDLcoefs_r = spectra_matrix .* swfc_k;
            modDLk(idx, :) = SYr * modDLcoefs_r;
        end

        if isReal, modDLk = real(modDLk); end
        modDL(:, :, k) = modDLk;
    end
end
end %% END MAIN FUNCTION

function [lambda_int, lambda_surf, lambda_ext] = LOCAL_modDLPspectrum(p, u0, a, oblate, c)
    sp = (p + 1)^2;

    if oblate
        error('Not implemented.')
    else
        cnm = 1j*c*(u0^2 - 1);

        Rnm1_vec = zeros(1, sp);
        Rnm3_vec = zeros(1, sp);
        dRnm1_vec = zeros(1, sp);
        dRnm3_vec = zeros(1, sp);

        for j = 0:p
            idx = j^2 + 1:(j + 1)^2;
            [R1_j, dR1_j] = Rnm1(j, [], u0, c);
            [R3_j, dR3_j] = Rnm3(j, [], u0, c);
            Rnm1_vec(idx) = R1_j;
            dRnm1_vec(idx) = dR1_j;
            Rnm3_vec(idx) = R3_j;
            dRnm3_vec(idx) = dR3_j;
        end
    end
    if oblate
        error('Not implemented.')
    else
        lambda_int = cnm .* dRnm3_vec;
        lambda_surf = cnm .* (dRnm1_vec .* Rnm3_vec + dRnm3_vec .* Rnm1_vec) / 2;
        lambda_ext = cnm .* dRnm1_vec;
    end
end

function F = LOCAL_solid_swf(p, u0, u_x, oblate, c)
    F = ones(size(u_x, 1), (p+1)^2);
    if oblate, error('Not implemented.'); end
    idx_int = abs(u_x) - u0 < -1e-14;
    idx_ext = abs(u_x) - u0 > 1e-14;

    if any(idx_int)
        F(idx_int, :) = radial_block_pswf(p, u_x(idx_int), c, 1);
    end
    if any(idx_ext)
        F(idx_ext, :) = radial_block_pswf(p, u_x(idx_ext), c, 2);
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
