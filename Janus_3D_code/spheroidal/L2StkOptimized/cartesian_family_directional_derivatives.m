function [dplus, dminus, dzero] = cartesian_family_directional_derivatives(p, u0, a, oblate, Gshc, X_trg, target_precomp)
%{
Evaluate the batched Cartesian-family directional derivatives
    dplus  = (d/dx + i d/dy) S^L
    dminus = (d/dx - i d/dy) S^L
    dzero  = d/dz S^L
for scalar spheroidal SLP modal coefficients Gshc.

This is an off-surface evaluator specialized to the directions used by the
Cartesian TSL backend. It avoids the general directional-SP machinery in
spheroidalSPOptimized.
%}
[~, nf, ns] = size(Gshc);
nt = size(X_trg, 1);

dplus = zeros(nt, nf, ns);
dminus = zeros(nt, nf, ns);
dzero = zeros(nt, nf, ns);

for k = 1:ns
    source_precomp = target_precomp.sources(k);
    if source_precomp.nt == 0, continue; end

    for r = 1:numel(source_precomp.regions)
        region = source_precomp.regions(r);
        idx = region.idx;

        [spec_plus_p, spec_plus_0, spec_plus_1, ...
            spec_minus_p, spec_minus_0, spec_minus_1, ...
            spec_zero_p, spec_zero_0, spec_zero_1] = ...
            LOCAL_directional_spectra_from_weights( ...
                region.u, region.v, region.phi, source_precomp.oblate, ...
                region.source_radial_prefactors, region.source_radial_derivative_prefactors, ...
                target_precomp.nrow, target_precomp.mrow, target_precomp.n1row);

        dplus(idx, :, k) = ...
            ((spec_plus_p .* region.Fp + spec_plus_0 .* region.Fr) .* region.Yr0 + ...
                spec_plus_1 .* region.Fr .* region.Yr_nplus1) * Gshc(:, :, k);
        dminus(idx, :, k) = ...
            ((spec_minus_p .* region.Fp + spec_minus_0 .* region.Fr) .* region.Yr0 + ...
                spec_minus_1 .* region.Fr .* region.Yr_nplus1) * Gshc(:, :, k);
        dzero(idx, :, k) = ...
            ((spec_zero_p .* region.Fp + spec_zero_0 .* region.Fr) .* region.Yr0 + ...
                spec_zero_1 .* region.Fr .* region.Yr_nplus1) * Gshc(:, :, k);
    end
end
end %% END MAIN FUNCTION

function [spec_plus_p, spec_plus_0, spec_plus_1, ...
          spec_minus_p, spec_minus_0, spec_minus_1, ...
          spec_zero_p, spec_zero_0, spec_zero_1] = ...
          LOCAL_directional_spectra_from_weights(u_x, v_x, phi_x, oblate, ...
              source_radial_prefactors, source_radial_derivative_prefactors, ...
              nrow, mrow, n1row)
    nt = numel(u_x);

    if oblate
        root_uv = sqrt(u_x.^2 + v_x.^2);
        root_u = sqrt(u_x.^2 + 1);
        root_v = sqrt(1 - v_x.^2);
        phase_plus = exp(1i .* phi_x);
        phase_minus = conj(phase_plus);

        nu_plus_u = (u_x .* root_v ./ root_uv) .* phase_plus;
        nu_plus_v = (-v_x .* root_u ./ root_uv) .* phase_plus;
        nu_plus_phi = 1i .* phase_plus;

        nu_minus_u = (u_x .* root_v ./ root_uv) .* phase_minus;
        nu_minus_v = (-v_x .* root_u ./ root_uv) .* phase_minus;
        nu_minus_phi = -1i .* phase_minus;

        nu_zero_u = v_x .* root_u ./ root_uv;
        nu_zero_v = u_x .* root_v ./ root_uv;
    else
        root_uv = sqrt(u_x.^2 - v_x.^2);
        root_u = sqrt(u_x.^2 - 1);
        root_v = sqrt(1 - v_x.^2);
        phase_plus = exp(1i .* phi_x);
        phase_minus = conj(phase_plus);

        nu_plus_u = (u_x .* root_v ./ root_uv) .* phase_plus;
        nu_plus_v = (-v_x .* root_u ./ root_uv) .* phase_plus;
        nu_plus_phi = 1i .* phase_plus;

        nu_minus_u = (u_x .* root_v ./ root_uv) .* phase_minus;
        nu_minus_v = (-v_x .* root_u ./ root_uv) .* phase_minus;
        nu_minus_phi = -1i .* phase_minus;

        nu_zero_u = v_x .* root_u ./ root_uv;
        nu_zero_v = u_x .* root_v ./ root_uv;
    end

    [spec_plus_p, spec_plus_0, spec_plus_1] = LOCAL_build_directional_spectra( ...
        u_x, v_x, nu_plus_u, nu_plus_v, nu_plus_phi, oblate, ...
        source_radial_prefactors, source_radial_derivative_prefactors, ...
        nrow, mrow, n1row);
    [spec_minus_p, spec_minus_0, spec_minus_1] = LOCAL_build_directional_spectra( ...
        u_x, v_x, nu_minus_u, nu_minus_v, nu_minus_phi, oblate, ...
        source_radial_prefactors, source_radial_derivative_prefactors, ...
        nrow, mrow, n1row);
    [spec_zero_p, spec_zero_0, spec_zero_1] = LOCAL_build_directional_spectra( ...
        u_x, v_x, nu_zero_u, nu_zero_v, zeros(nt, 1), oblate, ...
        source_radial_prefactors, source_radial_derivative_prefactors, ...
        nrow, mrow, n1row);
end

function [spec_prime, spec_zero, spec_n1] = LOCAL_build_directional_spectra( ...
        u_x, v_x, nu_u, nu_v, nu_phi, oblate, ...
        source_radial_prefactors, source_radial_derivative_prefactors, ...
        nrow, mrow, n1row)
    if oblate
        denom_uv = sqrt((u_x.^2 + v_x.^2) .* (1 - v_x.^2));
        denom_phi = sqrt((u_x.^2 + 1) .* (1 - v_x.^2));
        col_prime = sqrt((u_x.^2 + 1) ./ (u_x.^2 + v_x.^2)) .* nu_u;
    else
        denom_uv = sqrt((u_x.^2 - v_x.^2) .* (1 - v_x.^2));
        denom_phi = sqrt((u_x.^2 - 1) .* (1 - v_x.^2));
        col_prime = sqrt((u_x.^2 - 1) ./ (u_x.^2 - v_x.^2)) .* nu_u;
    end

    col_v = v_x .* nu_v ./ denom_uv;
    col_phi = nu_phi ./ denom_phi;
    col_n1 = nu_v ./ denom_uv;

    spec_prime = col_prime * source_radial_derivative_prefactors;
    spec_zero = col_v * (source_radial_prefactors .* nrow) + col_phi * (1i .* source_radial_prefactors .* mrow);
    spec_n1 = -col_n1 * (source_radial_prefactors .* n1row);
end
