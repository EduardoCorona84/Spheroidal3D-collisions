function [GDSL_U, GDSL_V, GDSL_PHI] = spheroidalgraddivSLOptimized(p, u0, a, oblate, Gshc_x, Gshc_y, Gshc_z, isReal, X_trg)
%{
Optimized spheroidal graddiv SL evaluation.
%}

if nargin < 9
    X_trg = [];
end

[sp, nf, ns] = size(Gshc_x);
if sp ~= (p + 1)^2
    error("Gshc_x first dimension must be (p+1)^2.");
end
if ~isequal(size(Gshc_x), size(Gshc_y)) || ~isequal(size(Gshc_x), size(Gshc_z))
    error("Gshc_x, Gshc_y, and Gshc_z must have the same size.");
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
    [Yr, v_k, phi_k] = LOCAL_precompute_gd_grid(p);
    nt = length(v_k);

    GDSL_U = zeros(nt, nf, ns, "like", Gshc_x);
    GDSL_V = zeros(nt, nf, ns, "like", Gshc_x);
    GDSL_PHI = zeros(nt, nf, ns, "like", Gshc_x);

    for k = 1:ns
        u_k = u0(k) .* ones(nt, 1, "like", v_k);
        [Ucomp, Vcomp, PHIcomp] = LOCAL_graddivSL_eval( ...
            p, u0(k), a(k), u_k, v_k, phi_k, ...
            Gshc_x(:, :, k), Gshc_y(:, :, k), Gshc_z(:, :, k), oblate(k), Yr);

        if isReal
            Ucomp = real(Ucomp);
            Vcomp = real(Vcomp);
            PHIcomp = real(PHIcomp);
        end

        GDSL_U(:, :, k) = Ucomp;
        GDSL_V(:, :, k) = Vcomp;
        GDSL_PHI(:, :, k) = PHIcomp;
    end
else
    if size(X_trg, 2) ~= 3 || size(X_trg, 3) ~= ns
        error("X_trg must be nt x 3 x ns for numeric-only mode.");
    end

    nt = size(X_trg, 1);
    GDSL_U = zeros(nt, nf, ns, "like", Gshc_x);
    GDSL_V = zeros(nt, nf, ns, "like", Gshc_x);
    GDSL_PHI = zeros(nt, nf, ns, "like", Gshc_x);

    for k = 1:ns
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

        GDSL_k_U = zeros(ntk, nf, "like", Gshc_x);
        GDSL_k_V = zeros(ntk, nf, "like", Gshc_x);
        GDSL_k_PHI = zeros(ntk, nf, "like", Gshc_x);

        for r = 1:3
            idx = regions{r};
            if ~any(idx)
                continue;
            end

            Sr = S(idx, :);
            u_r = Sr(:, 1);
            v_r = Sr(:, 2);
            phi_r = Sr(:, 3);

            v_r_real = real(v_r);
            if norm(abs(v_r_real - v_r)) > 1e-10
                fprintf("v_r imaginary\n ");
            end
            v_r = v_r_real;

            [Ucomp, Vcomp, PHIcomp] = LOCAL_graddivSL_eval( ...
                p, u0(k), a(k), u_r, v_r, phi_r, ...
                Gshc_x(:, :, k), Gshc_y(:, :, k), Gshc_z(:, :, k), oblate(k), []);

            if isReal
                Ucomp = real(Ucomp);
                Vcomp = real(Vcomp);
                PHIcomp = real(PHIcomp);
            end

            GDSL_k_U(idx, :) = Ucomp;
            GDSL_k_V(idx, :) = Vcomp;
            GDSL_k_PHI(idx, :) = PHIcomp;
        end

        GDSL_U(:, :, k) = GDSL_k_U;
        GDSL_V(:, :, k) = GDSL_k_V;
        GDSL_PHI(:, :, k) = GDSL_k_PHI;
    end
end
end %% END MAIN FUNCTION

function [Ucomponent, Vcomponent, PHIcomponent] = LOCAL_graddivSL_eval(p, u0, a, u, v, phi, Gshc_x, Gshc_y, Gshc_z, oblate, Yr)
    nt_r = length(u);

    if isempty(Yr)
        Yr = LOCAL_build_Yr(p, v, phi);
    end

    [nn, mm] = LOCAL_get_nn_mm(p);

    if oblate
        cnm = 1j .* a .* factorial(nn - mm) ./ factorial(nn + mm) .* (-1).^(mm) .* sqrt(u0.^2 + 1);
        L = legendre_otc(p, 1j .* u0, 1, 1, 1);
        if all(abs(u) - u0 < -1e-12)
            gnm = L{2};
        elseif all(abs(u) - u0 > 1e-12)
            gnm = L{1};
        else
            gnm = ones(size(L{1}));
        end
        [Fr, Fp, Fpp] = solid_harmonic_prime(p, u0, 1j .* u);
        common_coeffs = a^(-2) .* cnm .* gnm;
    else
        bnm = a .* factorial(nn - mm) ./ factorial(nn + mm) .* (-1).^(mm) .* sqrt(u0.^2 - 1);
        L = legendre_otc(p, u0, 1, 1, 1);
        if all(abs(u) - u0 < -1e-12)
            gnm = L{2};
        elseif all(abs(u) - u0 > 1e-12)
            gnm = L{1};
        else
            gnm = ones(size(L{1}));
        end
        [Fr, Fp, Fpp] = solid_harmonic_prime(p, u0, u);
        common_coeffs = bnm .* gnm ./ (a.^2);
    end

    coeffs = spheroidalgraddivSLcoefficients(u, v, phi, nn', mm', oblate);
    gshc_types = {'gshcx', 'gshcy', 'gshcz'};
    Gshc_data = {Gshc_x, Gshc_y, Gshc_z};

    Yr0 = Yr(1:nt_r, :);
    Yr1 = Yr(nt_r+1:2*nt_r, :);
    Yr2 = Yr(2*nt_r+1:end, :);

    function gshc_coeff = calculate_gshc_coeff(coeffs_for_gshc)
        f_nm  = coeffs_for_gshc.Ynm;
        f_n1m = coeffs_for_gshc.Yn1m;
        f_n2m = coeffs_for_gshc.Yn2m;

        Ynm_coeff  = f_nm{1}  .* Fr + f_nm{2}  .* Fp + f_nm{3}  .* Fpp;
        Yn1m_coeff = f_n1m{1} .* Fr + f_n1m{2} .* Fp + f_n1m{3} .* Fpp;
        Yn2m_coeff = f_n2m{1} .* Fr + f_n2m{2} .* Fp + f_n2m{3} .* Fpp;

        gshc_coeff = Ynm_coeff  .* Yr0 + ...
                    Yn1m_coeff .* Yr1 + ...
                    Yn2m_coeff .* Yr2;
    end

    Ucomponent = zeros(nt_r, size(Gshc_x, 2), "like", Gshc_x);
    for i = 1:3
        type = gshc_types{i};
        Gshc_coeff = calculate_gshc_coeff(coeffs.U.(type));
        Ucomponent = Ucomponent + (common_coeffs.' .* Gshc_coeff) * Gshc_data{i};
    end

    Vcomponent = zeros(nt_r, size(Gshc_x, 2), "like", Gshc_x);
    for i = 1:3
        type = gshc_types{i};
        Gshc_coeff = calculate_gshc_coeff(coeffs.V.(type));
        Vcomponent = Vcomponent + (common_coeffs.' .* Gshc_coeff) * Gshc_data{i};
    end

    PHIcomponent = zeros(nt_r, size(Gshc_x, 2), "like", Gshc_x);
    for i = 1:3
        type = gshc_types{i};
        Gshc_coeff = calculate_gshc_coeff(coeffs.PHI.(type));
        PHIcomponent = PHIcomponent + (common_coeffs.' .* Gshc_coeff) * Gshc_data{i};
    end
end

function [Yr, v_k, phi_k] = LOCAL_precompute_gd_grid(p)
    persistent cache_p cache_Yr cache_vk cache_phik
    if ~isempty(cache_p) && cache_p == p
        Yr = cache_Yr;
        v_k = cache_vk;
        phi_k = cache_phik;
        return;
    end

    [theta2, phi_k] = gl_grid(p);
    v_k = cos(theta2);
    v_k = v_k(:);
    phi_k = phi_k(:);
    nt = length(v_k);
    sp = (p + 1)^2;

    Yr = zeros(nt * 3, sp);
    v_row = real(acos(v_k)).';
    for n = 0:p
        Yn = Ynm(n, [], v_row, phi_k);
        Yr(1:nt, n^2+1:(n+1)^2) = Yn;
        Yn1 = Ynm(n+1, -n:n, v_row, phi_k);
        Yr(nt+1:2*nt, n^2+1:(n+1)^2) = Yn1;
        Yn2 = Ynm(n+2, -n:n, v_row, phi_k);
        Yr(2*nt+1:end, n^2+1:(n+1)^2) = Yn2;
    end

    cache_Yr = Yr;
    cache_vk = v_k;
    cache_phik = phi_k;
    cache_p = p;
end

function Yr = LOCAL_build_Yr(p, v, phi)
    nt = length(v);
    sp = (p + 1)^2;
    Yr = zeros(nt * 3, sp, "like", v);
    v_row = real(acos(v)).';

    for n = 0:p
        Yn = Ynm(n, [], v_row, phi);
        Yr(1:nt, n^2+1:(n+1)^2) = Yn;
        Yn1 = Ynm(n+1, -n:n, v_row, phi);
        Yr(nt+1:2*nt, n^2+1:(n+1)^2) = Yn1;
        Yn2 = Ynm(n+2, -n:n, v_row, phi);
        Yr(2*nt+1:end, n^2+1:(n+1)^2) = Yn2;
    end
end

function [nn, mm] = LOCAL_get_nn_mm(p)
    persistent cache_p cache_nn cache_mm
    if ~isempty(cache_p) && cache_p == p
        nn = cache_nn;
        mm = cache_mm;
        return;
    end

    sp = (p + 1)^2;
    ii = (1:sp)';
    nn = floor(sqrt(ii - 1));
    mm = ii - nn.^2 - nn - 1;

    cache_p = p;
    cache_nn = nn;
    cache_mm = mm;
end

function [Fr, Fp, Fpp] = solid_harmonic_prime(p, u0, u_x)
    %{
        Solid spheroidal harmonics can be written as f_n^m(u)Y_n^m(v, phi).
        This function returns what Fr = f_n^m is (and its derivative as Fp),
        depending on whether we are in the exterior or interior (or on the surface).
    %}
    if abs(u_x) - u0 < -1e-12 % Interior
        PQ = legendre_otc(p, u_x, 1, 2, 2);
        P = PQ{1}; dP = PQ{3}; ddP = PQ{5};
        Fr = P.'; Fp = dP.'; Fpp = ddP.';
    elseif abs(u_x) - u0 > 1e-12 % Exterior
        PQ = legendre_otc(p, u_x, 1, 2, 2);
        Q = PQ{2}; dQ = PQ{4}; ddQ = PQ{6};
        Fr = Q.'; Fp = dQ.'; Fpp = ddQ.';
    else % On surface: average of interior/exterior
        PQ = legendre_otc(p, u_x, 1, 2, 2);
        P = PQ{1}; dP = PQ{3}; ddP = PQ{5};
        Q = PQ{2}; dQ = PQ{4}; ddQ = PQ{6};
        Fr = P.' .* Q.';
        Fp = (Q.' .* dP.' + P.' .* dQ.') ./ 2;
        Fpp = (Q.' .* ddP.' + P.' .* ddQ.') ./ 2;
    end
end
