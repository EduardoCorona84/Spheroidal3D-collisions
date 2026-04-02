function target_precomp = spheroidal_cartesian_target_precompute(p, u0, a, oblate, X_trg, include_hessian, allow_surface, surface_is_interior)
%{
Build shared target-side cached data for the Cartesian TSL backend.

This context is intended primarily for off-surface target points. When
allow_surface is true, surface targets are treated using the exterior
branch so that we can add the jump relation separately.
%}
if nargin < 6 || isempty(include_hessian)
    include_hessian = false;
end
if nargin < 7 || isempty(allow_surface)
    allow_surface = false;
end
if nargin < 8 || isempty(surface_is_interior)
    surface_is_interior = false;
end

[~, ~, ns] = size(X_trg);
if isscalar(u0)
    u0 = u0 .* ones(1, ns);
end
if isscalar(a)
    a = a .* ones(1, ns);
end
if isscalar(oblate)
    oblate = oblate .* ones(1, ns);
end

[nn, mm, anm_base] = LOCAL_get_mode_data(p);

target_precomp = struct();
target_precomp.p = p;
target_precomp.nn = nn;
target_precomp.mm = mm;
target_precomp.nrow = (nn + 1).';
target_precomp.mrow = mm.';
target_precomp.n1row = (nn - mm + 1).';
target_precomp.anm_base = anm_base;
target_precomp.include_hessian = include_hessian;

source_template = struct(...
    'nt', 0, ...
    'oblate', false, ...
    'regions', [], ...
    'all_exterior', false ...
);
sources = repmat(source_template, 1, ns);
tol = 9e-12;

for k = 1:ns
    Xtk = X_trg(:, :, k);
    source_precomp = source_template;
    source_precomp.nt = size(Xtk, 1);
    source_precomp.oblate = oblate(k);

    if isempty(Xtk)
        sources(k) = source_precomp;
        continue;
    end

    % Classify the target points and then build one precompute block per region.
    S = cart2spheroidal(Xtk, a(k), oblate(k));
    u_x = S(:, 1);
    indices_interior = find(u_x < u0(k) - tol);
    indices_exterior = find(u_x > u0(k) + tol);
    indices_surface = ~(u_x < u0(k) - tol | u_x > u0(k) + tol);

    regions = repmat(LOCAL_empty_region(), 1, 0);

    % For one-sided self traces, surface nodes can be forced onto the
    % interior branch when building T^-: this is needed to calculate the average.
    if allow_surface && any(indices_surface) && surface_is_interior
        indices_interior = unique([indices_interior; find(indices_surface)]);
        indices_surface = false(size(indices_surface));
    end

    if ~isempty(indices_interior)
        regions(end + 1) = LOCAL_build_region_context( ...
            p, u0(k), a(k), oblate(k), S, indices_interior, true, include_hessian, nn, mm, anm_base);
    end

    % Surface nodes are grouped with the exterior branch when allow_surface
    % is enabled, so we can add any jump term separately.
    surface_or_exterior = unique([indices_exterior; find(indices_surface)]);
    if ~isempty(surface_or_exterior)
        regions(end + 1) = LOCAL_build_region_context( ...
            p, u0(k), a(k), oblate(k), S, surface_or_exterior, false, include_hessian, nn, mm, anm_base);
    end

    source_precomp.regions = regions;
    source_precomp.all_exterior = (numel(regions) == 1) && ~regions(1).is_interior;
    sources(k) = source_precomp;
end

target_precomp.sources = sources;
end %% END MAIN FUNCTION

function region = LOCAL_build_region_context(p, u0, a, oblate, S, idx, is_interior, include_hessian, nn, mm, anm_base)
    region = LOCAL_empty_region();

    u = S(idx, 1);
    v = real(S(idx, 2));
    phi = S(idx, 3);

    [Yr0, Yr1, Yr_nplus1, Yr2] = LOCAL_build_Yr(p, v, phi);
    [Fr, Fp, Fpp] = LOCAL_solid_harmonic(p, u, is_interior, oblate);
    [source_radial_prefactors, source_radial_derivative_prefactors, common_coeffs] = ...
        LOCAL_source_radial_prefactors(p, u0, a, oblate, is_interior, anm_base);

    region.idx = idx;
    region.is_interior = is_interior;
    region.a = a;
    region.oblate = oblate;
    region.u = u;
    region.v = v;
    region.phi = phi;
    region.Yr0 = Yr0;
    region.Yr1 = Yr1;
    region.Yr_nplus1 = Yr_nplus1;
    region.Yr2 = Yr2;
    region.Fr = Fr;
    region.Fp = Fp;
    region.Fpp = Fpp;
    region.source_radial_prefactors = source_radial_prefactors;
    region.source_radial_derivative_prefactors = source_radial_derivative_prefactors;
    region.common_coeffs = common_coeffs;
end

function region = LOCAL_empty_region()
    region = struct( ...
        'idx', [], ...
        'is_interior', false, ...
        'a', [], ...
        'oblate', false, ...
        'u', [], ...
        'v', [], ...
        'phi', [], ...
        'Yr0', [], ...
        'Yr1', [], ...
        'Yr_nplus1', [], ...
        'Yr2', [], ...
        'Fr', [], ...
        'Fp', [], ...
        'Fpp', [], ...
        'source_radial_prefactors', [], ...
        'source_radial_derivative_prefactors', [], ...
        'common_coeffs', []);
end

function [source_radial_prefactors, source_radial_derivative_prefactors, common_coeffs] = ...
        LOCAL_source_radial_prefactors(p, u0, a, oblate, is_interior, anm_base)
    % Retrieve values associated with the eigenvalues of the spheroidal Laplace SLP
    if oblate
        PQ0 = legendre_otc(p, 1i .* u0, 1);
        if is_interior
            gnm = PQ0{2};
        else
            gnm = PQ0{1};
        end

        source_radial_prefactors = (1i .* anm_base .* sqrt(u0.^2 + 1) .* gnm(:)).';
        source_radial_derivative_prefactors = 1i .* source_radial_prefactors;
    else
        PQ0 = legendre_otc(p, u0, 1, 1, 1);
        if is_interior
            gnm = PQ0{2};
        else
            gnm = PQ0{1};
        end

        source_radial_prefactors = (anm_base .* sqrt(u0.^2 - 1) .* gnm(:)).';
        source_radial_derivative_prefactors = source_radial_prefactors;
    end

    common_coeffs = source_radial_prefactors.' ./ a;
end

function [Fr, Fp, Fpp] = LOCAL_solid_harmonic(p, u_x, is_interior, oblate)
    eval_u = u_x;
    if oblate
        eval_u = 1i .* eval_u;
    end

    PQ = legendre_otc(p, eval_u, 1, 2, 2);
    if is_interior
        Fr = PQ{1}.';
        Fp = PQ{3}.';
        Fpp = PQ{5}.';
    else
        Fr = PQ{2}.';
        Fp = PQ{4}.';
        Fpp = PQ{6}.';
    end
end

function [Yr0, Yr1, Yr_nplus1, Yr2] = LOCAL_build_Yr(p, v, phi)
    nt = numel(v);
    sp = (p + 1)^2;
    Yr0 = zeros(nt, sp);
    Yr1 = zeros(nt, sp);
    Yr_nplus1 = zeros(nt, sp);
    Yr2 = zeros(nt, sp);

    theta_row = real(acos(v)).';
    for n = 0:p
        idx = n^2 + 1:(n + 1)^2;
        mvals = -n:n;

        Yr0(:, idx) = Ynm(n, [], theta_row, phi);

        Yn1 = Ynm(n + 1, mvals, theta_row, phi);
        Yr1(:, idx) = Yn1;
        scale1 = sqrt((2 * n + 1) / (2 * n + 3) .* (n + mvals + 1) ./ (n - mvals + 1));
        Yr_nplus1(:, idx) = scale1 .* Yn1;

        Yn2 = Ynm(n + 2, mvals, theta_row, phi);
        Yr2(:, idx) = Yn2;
    end
end

function [nn, mm, anm_base] = LOCAL_get_mode_data(p)
    persistent cache_p cache_nn cache_mm cache_anm_base
    if ~isempty(cache_p) && cache_p == p
        nn = cache_nn;
        mm = cache_mm;
        anm_base = cache_anm_base;
        return;
    end

    sp = (p + 1)^2;
    ii = (1:sp)';
    nn = floor(sqrt(ii - 1));
    mm = ii - nn.^2 - nn - 1;
    anm_base = factorial(nn - mm) ./ factorial(nn + mm) .* (-1).^mm;

    cache_p = p;
    cache_nn = nn;
    cache_mm = mm;
    cache_anm_base = anm_base;
end
