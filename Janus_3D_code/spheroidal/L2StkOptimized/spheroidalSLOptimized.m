function SL = spheroidalSLOptimized(p, u0, a, oblate, Gshc, isReal, X_trg)
%{
Optimized spheroidal SL evaluation.
%}

if nargin < 7
    X_trg = [];
end

[sp, nf, ns] = size(Gshc);
if sp ~= (p + 1)^2
    error("Gshc first dimension must be (p+1)^2.");
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
    Y = LOCAL_precompute_sl_grid(p);
    nt = size(Y, 1);

    [~, spectra_surf, ~] = LOCAL_SLspectrum(p, u0, a, oblate);
    spectra_surf = reshape(spectra_surf, sp, 1, ns);
    spectra_matrix = repmat(spectra_surf, 1, nf, 1);

    SL_coefs = spectra_matrix .* Gshc;
    SL = Y * reshape(SL_coefs, sp, []);
    SL = reshape(SL, nt, nf, ns);

    if isReal
        SL = real(SL);
    end
else
    if ~isnumeric(X_trg) || ndims(X_trg) ~= 3 || size(X_trg, 2) ~= 3 || size(X_trg, 3) ~= ns
        error("X_trg must be nt x 3 x ns for numeric-only mode.");
    end

    nt = size(X_trg, 1);
    SL = zeros(nt, nf, ns, "like", Gshc);
    [spectra_int, spectra_surf, spectra_ext] = LOCAL_SLspectrum(p, u0, a, oblate);

    for k = 1:ns
        Xtk = X_trg(:, :, k);
        ntk = size(Xtk, 1);
        if ntk == 0, continue; end

        S = cart2spheroidal(Xtk, a(k), oblate(k));
        u_x = S(:, 1);

        indices_interior = (u_x < u0(k) - 9e-12);
        indices_surface = (abs(u_x - u0(k)) <= 9e-12);
        indices_exterior = (u_x > u0(k) + 9e-12);

        regions = {indices_interior, indices_surface, indices_exterior};
        spectra_regions = {spectra_int(:, k), spectra_surf(:, k), spectra_ext(:, k)};

        SLk = zeros(ntk, nf, "like", Gshc);
        for r = 1:3
            idx = regions{r};
            if ~any(idx)
                continue;
            end

            Sr = S(idx, :);
            u_x_r = Sr(:, 1);
            v_x_r = Sr(:, 2);
            phi_x_r = Sr(:, 3);

            Fr = LOCAL_solid_harmonic(p, u0(k), u_x_r, oblate(k));
            nt_r = length(u_x_r);
            Yr = zeros(nt_r, sp, "like", Gshc);
            v_row = real(acos(v_x_r));

            for n = 0:p
                Yn = Ynm(n, [], v_row, phi_x_r);
                Yr(:, n^2+1:(n+1)^2) = Yn;
            end

            FYr = Fr .* Yr;
            spectra_matrix = repmat(spectra_regions{r}, 1, nf);
            SLcoefs_r = spectra_matrix .* Gshc(:, :, k);
            SLk(idx, :) = FYr * SLcoefs_r;
        end

        if isReal
            SLk = real(SLk);
        end
        SL(:, :, k) = SLk;
    end
end
end %% END MAIN FUNCTION

function [lambda_int, lambda_surf, lambda_ext] = LOCAL_SLspectrum(p, u0, a, oblate)
    sp = (p + 1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm = ii - nn.^2 - nn - 1;
    cnm = a .* factorial(nn-mm) ./ factorial(nn+mm) .* (-1).^(mm) .* (oblate .* 1j .* sqrt(u0.^2+1) + ~oblate .* sqrt(u0.^2-1));
    L_obl = legendre_otc(p, 1j .* u0, 1, 1, 1);
    u0_pro = u0 + oblate;
    L_pro = legendre_otc(p, u0_pro, 1, 1, 1);
    Pp = L_pro{1}; Qp = L_pro{2};
    Po = L_obl{1}; Qo = L_obl{2};
    lambda_int = cnm .* (Qo * diag(oblate) + Qp * diag(~oblate));
    lambda_surf = cnm .* ((Po .* Qo) * diag(oblate) + (Pp .* Qp) * diag(~oblate));
    lambda_ext = cnm .* (Po * diag(oblate) + Pp * diag(~oblate));
end

function F = LOCAL_solid_harmonic(p, u0, u_x, oblate)
    F = ones(size(u_x, 1), (p + 1)^2);
    if oblate
        u_x = 1j .* u_x;
    end
    if abs(u_x) - u0 < -1e-14
        PQ = legendre_otc(p, u_x);
        P = PQ{1};
        F = P.';
    elseif abs(u_x) - u0 > 1e-14
        PQ = legendre_otc(p, u_x, 1);
        Q = PQ{2};
        F = Q.';
    end
end

function Y = LOCAL_precompute_sl_grid(p)
    % Precompute Ynm blocks used in self-eval branch (cached by p).
    persistent cache_p cache_Y
    if ~isempty(cache_p) && cache_p == p
        Y = cache_Y;
        return;
    end

    [theta2, phi_k] = gl_grid(p);
    v_k = cos(theta2);
    v_k = v_k(:);
    phi_k = phi_k(:);
    nt_r = length(v_k);
    sp = (p + 1)^2;

    Y = zeros(nt_r, sp);
    v_row = real(acos(v_k)).';
    for n = 0:p
        Yn = Ynm(n, [], v_row, phi_k);
        Y(:, n^2+1:(n+1)^2) = Yn;
    end

    cache_Y = Y;
    cache_p = p;
end

