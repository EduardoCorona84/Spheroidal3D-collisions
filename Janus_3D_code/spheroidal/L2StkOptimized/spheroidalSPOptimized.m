function [SPx, SPy, SPz] = spheroidalSPOptimized(p, u0, a, oblate, Gshc, isReal, nu_x, nu_y, nu_z, X_trg)
%{
Optimized spheroidal SP evaluation.
%}

[~, nf, ns] = size(Gshc);

if isempty(X_trg)
    [Yr, v_k, phi_k] = LOCAL_precompute_sp_grid(p);
    nt_r = length(v_k);

    SPx = zeros(nt_r, nf, ns);
    SPy = zeros(nt_r, nf, ns);
    SPz = zeros(nt_r, nf, ns);

    for k = 1:ns
        u_k = u0(k) .* ones(nt_r, 1, "like", v_k);
        S = [u_k, v_k, phi_k];

        nu_cart_list = {nu_x, nu_y, nu_z};
        for nu_ind = 1:3
            nu_cart = nu_cart_list{nu_ind}(:, :, k);

            [nu_sph, ~] = cartNu2spheroidal(nu_cart, S, a(k), oblate(k));
            [spectra_nm_prime, spectra_nm, spectra_n1m] = LOCAL_SPspectrum_away( ...
                p, u0(k), u_k, v_k, nu_sph, oblate(k));

            FYr = (spectra_nm_prime + spectra_nm) .* Yr(1:nt_r, :) + spectra_n1m .* Yr(nt_r+1:end, :);
            SP_k = FYr * Gshc(:, :, k);

            if isReal
                SP_k = real(SP_k);
            end

            switch nu_ind
                case 1
                    SPx(:, :, k) = SP_k;
                case 2
                    SPy(:, :, k) = SP_k;
                case 3
                    SPz(:, :, k) = SP_k;
            end
        end
    end
else
    if ~isnumeric(X_trg) || ndims(X_trg) ~= 3 || size(X_trg, 2) ~= 3 || size(X_trg, 3) ~= ns
        error("X_trg must be nt x 3 x ns for numeric-only mode.");
    end

    nt = size(X_trg, 1);
    SPx = zeros(nt, nf, ns);
    SPy = zeros(nt, nf, ns);
    SPz = zeros(nt, nf, ns);

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

        nu_cart_list = {nu_x, nu_y, nu_z};
        for nu_ind = 1:3
            SPk = zeros(ntk, nf, "like", Gshc);
            nu_cart_k = nu_cart_list{nu_ind}(:, :, k);

            for r = 1:3
                idx = regions{r};
                if ~any(idx)
                    continue;
                end

                Sr = S(idx, :);
                u_x_r = Sr(:, 1);
                v_x_r = Sr(:, 2);
                phi_x_r = Sr(:, 3);
                nu_r_cart = nu_cart_k(idx, :);

                [nu_r_sph, ~] = cartNu2spheroidal(nu_r_cart, Sr, a(k), oblate(k));
                [spectra_nm_prime, spectra_nm, spectra_n1m] = LOCAL_SPspectrum_away( ...
                    p, u0(k), u_x_r, v_x_r, nu_r_sph, oblate(k));

                [Fr, Fp] = LOCAL_solid_harmonic_prime(p, u0(k), u_x_r, oblate(k));
                nt_r = length(u_x_r);
                Yr = zeros(nt_r * 2, (p + 1)^2, "like", Gshc);

                for n = 0:p
                    Yn = Ynm(n, [], real(acos(v_x_r))', phi_x_r);
                    Yr(1:nt_r, n^2+1:(n+1)^2) = Yn;
                    Yn1 = Ynm(n+1, -n:n, real(acos(v_x_r))', phi_x_r);
                    yn1_scale = sqrt((2*n+1)/(2*n+3).*(n+(-n:n)+1)./(n-(-n:n)+1));
                    Yr(nt_r+1:end, n^2+1:(n+1)^2) = yn1_scale .* Yn1;
                end

                FYr = (spectra_nm_prime .* Fp + spectra_nm .* Fr) .* Yr(1:nt_r, :) + spectra_n1m .* Fr .* Yr(nt_r+1:end, :);
                SPk(idx, :) = FYr * Gshc(:, :, k);
            end

            if isReal
                SPk = real(SPk);
            end

            switch nu_ind
                case 1
                    SPx(:, :, k) = SPk;
                case 2
                    SPy(:, :, k) = SPk;
                case 3
                    SPz(:, :, k) = SPk;
            end
        end
    end
end
end %% END MAIN FUNCTION

function [F, Fp] = LOCAL_solid_harmonic_prime(p, u0, u_x, oblate)
    F = ones(size(u_x, 1), (p + 1)^2);
    Fp = ones(size(u_x, 1), (p + 1)^2);
    if nargin < 4
        oblate = false;
    end
    if oblate
        u_x = 1j .* u_x;
    end
    if abs(u_x) - u0 < -1e-14
        PQ = legendre_otc(p, u_x, 1, 1);
        P = PQ{1}; dP = PQ{3};
        F = P.'; Fp = dP.';
    elseif abs(u_x) - u0 > 1e-14
        PQ = legendre_otc(p, u_x, 1, 1, 1);
        Q = PQ{2}; dQ = PQ{4};
        F = Q.'; Fp = dQ.';
    end
end

function [Yr, v_k, phi_k] = LOCAL_precompute_sp_grid(p)
    % Precompute Ynm blocks for self-eval branch.
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
    nt_r = length(v_k);
    sp = (p + 1)^2;

    Yr = zeros(nt_r * 2, sp);
    v_row = real(acos(v_k)).';
    for n = 0:p
        Yn = Ynm(n, [], v_row, phi_k);
        Yr(1:nt_r, n^2+1:(n+1)^2) = Yn;
        Yn1 = Ynm(n+1, -n:n, v_row, phi_k);
        yn1_scale = sqrt((2*n+1)/(2*n+3).*(n+(-n:n)+1)./(n-(-n:n)+1));
        Yr(nt_r+1:end, n^2+1:(n+1)^2) = yn1_scale .* Yn1;
    end

    cache_Yr = Yr;
    cache_vk = v_k;
    cache_phik = phi_k;
    cache_p = p;
end

function [lambda_nm_prime,lambda_nm,lambda_n1m]=LOCAL_SPspectrum(p,u0,oblate)
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^(mm).*(oblate.*1j.*sqrt(u0.^2+1)+~oblate.*sqrt(u0.^2-1));
    L_obl=legendre_otc(p,1j.*u0,1,1,1);
    u0_pro=u0+oblate;
    L_pro=legendre_otc(p,u0_pro,1,1,1);
    Pp=L_pro{1}; Qp=L_pro{2}; dPp=L_pro{3}; dQp=L_pro{4};
    Po=L_obl{1}; Qo=L_obl{2}; dPo=L_obl{3}; dQo=L_obl{4};
    lambda_nm_prime=anm.*(oblate.*1j.*sqrt(u0.^2+1)+~oblate.*sqrt(u0.^2-1))./2.*((Po.*dQo+dPo.*Qo)*diag(oblate)+(Pp.*dQp+dPp.*Qp)*diag(~oblate));
    lambda_nm=anm.*(Po.*Qo*diag(oblate)+Pp.*Qp*diag(~oblate));
    lambda_n1m=anm.*(Po.*Qo*diag(oblate)+Pp.*Qp*diag(~oblate));
end

function [lambda_nm_prime, lambda_nm,lambda_n1m] = LOCAL_SPspectrum_away(p,u0,u_x,v_x,nu,oblate)
    % Local copy of spheroidalSP's SPspectrum_away.
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    anm_base=factorial(nn-mm)./factorial(nn+mm).*(-1).^mm;
    nu_u = nu(:,1); nu_v = nu(:,2); nu_phi = nu(:,3);
    if ~oblate
        if norm(abs(u_x)-u0)<9e-12
            [lambda_nm_prime,lambda_nm,lambda_n1m]=LOCAL_SPspectrum(p,u0,oblate);
            lambda_nm_prime = lambda_nm_prime.'./sqrt(u_x.^2-v_x.^2).*nu_u;
            lambda_nm = lambda_nm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2-1).*(1-v_x.^2)));
            lambda_n1m = -lambda_n1m.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2));
        else
            anm=anm_base.*sqrt(u0.^2-1);
            PQ0 = legendre_otc(p,u0,1,1,1);
            if abs(u_x)-u0<1e-14
                gnm=PQ0{2};
            elseif abs(u_x)-u0>1e-14
                gnm=PQ0{1};
            end
            lambda_nm_prime = anm.'.*gnm.'.*sqrt((u_x.^2-1)./(u_x.^2-v_x.^2)).*nu_u;
            lambda_nm = anm.'.*gnm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2-1).*(1-v_x.^2)));
            lambda_n1m = -anm.'.*gnm.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2));
        end
    else
        if norm(abs(u_x)-u0)<9e-12
            [lambda_nm_prime,lambda_nm,lambda_n1m]=LOCAL_SPspectrum(p,u0,oblate);
            lambda_nm_prime = lambda_nm_prime.'./sqrt(u_x.^2+v_x.^2).*nu_u;
            lambda_nm = lambda_nm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2+1).*(1-v_x.^2)));
            lambda_n1m = -lambda_n1m.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2));
        else
            anm=1j.*anm_base.*sqrt(u0.^2+1);
            PQ0 = legendre_otc(p,1j.*u0,1);
            if abs(u_x)-u0<1e-14
                gnm=PQ0{2};
            elseif abs(u_x)-u0>1e-14
                gnm=PQ0{1};
            end
            lambda_nm_prime = 1j.*anm.'.*gnm.'.*sqrt((u_x.^2+1)./(u_x.^2+v_x.^2)).*nu_u;
            lambda_nm = anm.'.*gnm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2+1).*(1-v_x.^2)));
            lambda_n1m = -anm.'.*gnm.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2));
        end
    end
end
